use algebra::{Field, PairingEngine};
use algebra::fft::EvaluationDomain;

use crate::groth16::{generator::KeypairAssembly, prover::ProvingAssignment};
use r1cs_core::{Index, SynthesisError};

use rayon::prelude::*;
use std::ops::{AddAssign, SubAssign};

pub(crate) struct R1CStoQAP;

impl R1CStoQAP {
    #[inline]
    pub(crate) fn instance_map_with_evaluation<E: PairingEngine>(
        assembly: &KeypairAssembly<E>,
        t: &E::Fr,
    ) -> Result<(Vec<E::Fr>, Vec<E::Fr>, Vec<E::Fr>, E::Fr, usize, usize), SynthesisError> {
        let domain_size = assembly.num_constraints + (assembly.num_inputs - 1) + 1;
        let domain = EvaluationDomain::<E::Fr>::new(domain_size)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let domain_size = domain.size();

        let zt = domain.evaluate_vanishing_polynomial(*t);

        // Evaluate all Lagrange polynomials
        let coefficients_time = start_timer!(|| "Evaluate Lagrange coefficients");
        let u = domain.evaluate_all_lagrange_coefficients(*t);
        end_timer!(coefficients_time);

        let qap_num_variables = (assembly.num_inputs - 1) + assembly.num_aux;

        let mut a = vec![E::Fr::zero(); qap_num_variables + 1];
        let mut b = vec![E::Fr::zero(); qap_num_variables + 1];
        let mut c = vec![E::Fr::zero(); qap_num_variables + 1];

        for i in 0..assembly.num_inputs {
            a[i] = u[assembly.num_constraints + i];
        }

        for i in 0..assembly.num_constraints {
            for &(ref coeff, index) in assembly.at[i].iter() {
                let index = match index {
                    Index::Input(i) => i,
                    Index::Aux(i) => assembly.num_inputs + i,
                };

                a[index] += &(u[i] * &coeff);
            }
            for &(ref coeff, index) in assembly.bt[i].iter() {
                let index = match index {
                    Index::Input(i) => i,
                    Index::Aux(i) => assembly.num_inputs + i,
                };

                b[index] += &(u[i] * &coeff);
            }
            for &(ref coeff, index) in assembly.ct[i].iter() {
                let index = match index {
                    Index::Input(i) => i,
                    Index::Aux(i) => assembly.num_inputs + i,
                };

                c[index] += &(u[i] * &coeff);
            }
        }

        Ok((a, b, c, zt, qap_num_variables, domain_size))
    }

    #[inline]
    pub(crate) fn witness_map<E: PairingEngine>(
        prover: &ProvingAssignment<E>,
        d1: &E::Fr,
        d2: &E::Fr,
        d3: &E::Fr,
    ) -> Result<(Vec<E::Fr>, Vec<E::Fr>, usize), SynthesisError> {
        #[inline]
        fn evaluate_constraint<E: PairingEngine>(
            terms: &[(E::Fr, Index)],
            assignment: &[E::Fr],
            num_input: usize,
        ) -> E::Fr {
            let mut acc = E::Fr::zero();
            for &(coeff, index) in terms {
                let val = match index {
                    Index::Input(i) => assignment[i],
                    Index::Aux(i) => assignment[num_input + i],
                };
                acc += &(val * &coeff);
            }
            acc
        }

        let zero = E::Fr::zero();
        let one = E::Fr::one();

        let mut full_input_assignment = prover.input_assignment.clone();
        full_input_assignment.extend(prover.aux_assignment.clone());

        let domain =
            EvaluationDomain::<E::Fr>::new(prover.num_constraints + (prover.num_inputs - 1) + 1)
                .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let domain_size = domain.size();

        let mut a = vec![zero; domain_size];
        let mut b = vec![zero; domain_size];
        a[..prover.num_constraints]
            .par_iter_mut()
            .zip(b[..prover.num_constraints].par_iter_mut())
            .zip((&prover.at).par_iter())
            .zip((&prover.bt).par_iter())
            .for_each(|(((a, b), at_i), bt_i)| {
                *a = evaluate_constraint::<E>(&at_i, &full_input_assignment, prover.num_inputs);
                *b = evaluate_constraint::<E>(&bt_i, &full_input_assignment, prover.num_inputs);
            });

        for i in 0..prover.num_inputs {
            a[prover.num_constraints + i] = if i > 0 { full_input_assignment[i] } else { one };
        }

        domain.ifft_in_place(&mut a);
        domain.ifft_in_place(&mut b);

        let mut h: Vec<E::Fr> = vec![zero; domain_size];
        h.par_iter_mut()
            .zip(&a)
            .zip(&b)
            .for_each(|((h_i, a_i), b_i)| *h_i *= &(*d2 * &a_i + &(*d1 * &b_i)));
        h[0].sub_assign(&d3);
        let d1d2 = *d1 * d2;
        h[0].sub_assign(&d1d2);
        h.push(d1d2);

        domain.coset_fft_in_place(&mut a);
        domain.coset_fft_in_place(&mut b);

        let mut ab = domain.mul_polynomials_in_evaluation_domain(&a, &b);
        drop(a);
        drop(b);

        let mut c = vec![zero; domain_size];
        c[..prover.num_constraints]
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, c)| {
                *c = evaluate_constraint::<E>(
                    &prover.ct[i],
                    &full_input_assignment,
                    prover.num_inputs,
                );
            });

        domain.ifft_in_place(&mut c);
        domain.coset_fft_in_place(&mut c);

        ab.par_iter_mut()
            .zip(c)
            .for_each(|(ab_i, c_i)| *ab_i -= &c_i);

        domain.divide_by_vanishing_poly_on_coset_in_place(&mut ab);
        domain.coset_ifft_in_place(&mut ab);

        h[..domain_size - 1]
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, e)| e.add_assign(&ab[i]));

        Ok((full_input_assignment, h, domain_size))
    }
}
