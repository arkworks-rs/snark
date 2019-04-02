use algebra::{
    fft::{
        domain::{EvaluationDomain, Scalar},
        multicore::Worker,
    },
    Field, PairingEngine,
};

use crate::{
    gm17::{generator::KeypairAssembly, prover::ProvingAssignment},
    Index, SynthesisError,
};

use rayon::prelude::*;
use std::ops::{AddAssign, SubAssign};

pub(crate) struct R1CStoSAP;

impl R1CStoSAP {
    #[inline]
    pub(crate) fn instance_map_with_evaluation<E: PairingEngine>(
        assembly: &KeypairAssembly<E>,
        t: &E::Fr,
    ) -> Result<(Vec<E::Fr>, Vec<E::Fr>, E::Fr, usize, usize), SynthesisError> {
        let domain = vec![
            Scalar::<E>(E::Fr::zero());
            2 * assembly.num_constraints + 2 * (assembly.num_inputs - 1) + 1
        ];
        let domain = EvaluationDomain::<E, _>::from_coeffs(domain)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let m = domain.m();

        let zt = domain.z(&t);

        // Evaluate all Lagrange polynomials
        let coefficients_time = timer_start!(|| "Evaluate Lagrange coefficients");
        let u = domain.evaluate_all_lagrange_coefficients(&t);
        timer_end!(coefficients_time);

        let sap_num_variables =
            2 * (assembly.num_inputs - 1) + assembly.num_aux + assembly.num_constraints;
        let extra_var_offset = (assembly.num_inputs - 1) + assembly.num_aux + 1;
        let extra_constr_offset = 2 * assembly.num_constraints;
        let extra_var_offset2 =
            (assembly.num_inputs - 1) + assembly.num_aux + assembly.num_constraints;

        let mut a = vec![E::Fr::zero(); sap_num_variables + 1];
        let mut c = vec![E::Fr::zero(); sap_num_variables + 1];

        for i in 0..assembly.num_constraints {
            let u_2i = u[2 * i];
            let u_2i_plus_1 = u[2 * i + 1];
            let u_add = u_2i + &u_2i_plus_1;
            let u_sub = u_2i - &u_2i_plus_1;

            for &(ref coeff, index) in assembly.at[i].iter() {
                let index = match index {
                    Index::Input(i) => i,
                    Index::Aux(i) => assembly.num_inputs + i,
                };

                a[index] += &(u_add * &coeff);
            }

            for &(ref coeff, index) in assembly.bt[i].iter() {
                let index = match index {
                    Index::Input(i) => i,
                    Index::Aux(i) => assembly.num_inputs + i,
                };

                a[index] += &(u_sub * &coeff);
            }

            for &(ref coeff, index) in assembly.ct[i].iter() {
                let index = match index {
                    Index::Input(i) => i,
                    Index::Aux(i) => assembly.num_inputs + i,
                };

                c[index] += &((u_2i * &coeff).double().double());
            }
            c[extra_var_offset + i].add_assign(&u_add);
        }

        a[0].add_assign(&u[extra_constr_offset]);
        c[0].add_assign(&u[extra_constr_offset]);

        for i in 1..(assembly.num_inputs - 1) + 1 {
            // First extra constraint

            a[i].add_assign(&u[extra_constr_offset + 2 * i - 1]);
            a[0].add_assign(&u[extra_constr_offset + 2 * i - 1]);

            let t_four = u[extra_constr_offset + 2 * i - 1].double().double();

            c[i].add_assign(&t_four);
            c[extra_var_offset2 + i].add_assign(&u[extra_constr_offset + 2 * i - 1]);

            // Second extra constraint

            a[i].add_assign(&u[extra_constr_offset + 2 * i]);
            a[0].sub_assign(&u[extra_constr_offset + 2 * i]);
            c[extra_var_offset2 + i].add_assign(&u[extra_constr_offset + 2 * i]);
        }

        Ok((a, c, zt, sap_num_variables, m))
    }

    #[inline]
    pub(crate) fn witness_map<E: PairingEngine>(
        prover: &ProvingAssignment<E>,
        d1: &E::Fr,
        d2: &E::Fr,
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

        let temp = prover
            .at
            .par_iter()
            .zip(&prover.bt)
            .map(|(a_i, b_i)| {
                let mut extra_var: E::Fr =
                    evaluate_constraint::<E>(&a_i, &full_input_assignment, prover.num_inputs);
                extra_var.sub_assign(&evaluate_constraint::<E>(
                    &b_i,
                    &full_input_assignment,
                    prover.num_inputs,
                ));
                extra_var.square_in_place();
                extra_var
            })
            .collect::<Vec<_>>();
        full_input_assignment.extend(temp);

        for i in 1..prover.num_inputs {
            let mut extra_var = full_input_assignment[i];
            extra_var.sub_assign(&one);
            extra_var.square_in_place();
            full_input_assignment.push(extra_var);
        }

        let m = EvaluationDomain::<E, Scalar<E>>::compute_m_from_num_coeffs(
            2 * prover.num_constraints + 2 * (prover.num_inputs - 1) + 1,
        )
        .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;

        let extra_constr_offset = 2 * prover.num_constraints;
        let extra_var_offset = prover.num_inputs + prover.num_aux;
        let extra_var_offset2 = prover.num_inputs + prover.num_aux + prover.num_constraints - 1;

        let mut a = vec![zero; m];
        a[..2 * prover.num_constraints]
            .par_chunks_mut(2)
            .zip(&prover.at)
            .zip(&prover.bt)
            .for_each(|((chunk, at_i), bt_i)| {
                chunk[0] =
                    evaluate_constraint::<E>(&at_i, &full_input_assignment, prover.num_inputs);
                chunk[0].add_assign(&evaluate_constraint::<E>(
                    &bt_i,
                    &full_input_assignment,
                    prover.num_inputs,
                ));

                chunk[1] =
                    evaluate_constraint::<E>(&at_i, &full_input_assignment, prover.num_inputs);
                chunk[1].sub_assign(&evaluate_constraint::<E>(
                    &bt_i,
                    &full_input_assignment,
                    prover.num_inputs,
                ));
            });
        a[extra_constr_offset] = one;
        for i in 1..prover.num_inputs {
            a[extra_constr_offset + 2 * i - 1] = full_input_assignment[i] + &one;
            a[extra_constr_offset + 2 * i] = full_input_assignment[i] - &one;
        }

        let worker = Worker::new();

        let mut a =
            EvaluationDomain::from_coeffs(a.iter().map(|s| Scalar::<E>(*s)).collect::<Vec<_>>())
                .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        a.ifft(&worker);

        let d1_double = d1.double();
        let mut h = vec![d1_double; m];
        h.par_iter_mut()
            .zip(a.as_ref())
            .for_each(|(h_i, a_i)| *h_i *= &a_i.0);
        h[0].sub_assign(&d2);
        let d1d1 = d1.square();
        h[0].sub_assign(&d1d1);
        h.push(d1d1);

        a.coset_fft(&worker);

        let mut aa = a.clone();
        aa.mul_assign(&worker, &a);

        let mut c = vec![zero; m];
        c[..2 * prover.num_constraints]
            .par_chunks_mut(2)
            .enumerate()
            .for_each(|(i, chunk)| {
                let mut tmp: E::Fr = evaluate_constraint::<E>(
                    &prover.ct[i],
                    &full_input_assignment,
                    prover.num_inputs,
                );
                tmp.double_in_place();
                tmp.double_in_place();

                let assignment = full_input_assignment[extra_var_offset + i];
                chunk[0] = tmp + &assignment;
                chunk[1] = assignment;
            });
        c[extra_constr_offset] = one;
        for i in 1..prover.num_inputs {
            let mut tmp = full_input_assignment[i];
            tmp.double_in_place();
            tmp.double_in_place();

            let assignment = full_input_assignment[extra_var_offset2 + i];
            c[extra_constr_offset + 2 * i - 1] = tmp + &assignment;
            c[extra_constr_offset + 2 * i] = assignment;
        }

        let mut c = EvaluationDomain::<E, _>::from_coeffs(
            c.iter().map(|s| Scalar::<E>(*s)).collect::<Vec<_>>(),
        )
        .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        c.ifft(&worker);
        c.coset_fft(&worker);

        let c = c.into_coeffs();
        aa.as_mut()
            .par_iter_mut()
            .zip(c)
            .for_each(|(aa_i, c_i)| aa_i.0 -= &c_i.0);

        aa.divide_by_z_on_coset(&worker);
        aa.icoset_fft(&worker);

        let aa = aa.into_coeffs();
        h[..m - 1]
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, e)| e.add_assign(&aa[i].0));

        Ok((full_input_assignment, h, m))
    }
}
