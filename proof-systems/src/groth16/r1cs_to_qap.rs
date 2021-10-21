use algebra::fft::domain::get_best_evaluation_domain;
use algebra::{Field, PairingEngine};

use crate::groth16::{generator::KeypairAssembly, prover::ProvingAssignment};
use r1cs_core::{ConstraintSystem, Index, SynthesisError};

use rayon::prelude::*;
use std::ops::Mul;

pub(crate) struct R1CStoQAP;

#[inline]
fn evaluate_constraint<'a, E: PairingEngine>(
    terms: &'a [(E::Fr, Index)],
    assignment: &'a [E::Fr],
    num_inputs: usize,
) -> E::Fr {
    terms
        .par_iter()
        .map(|(coeff, index)| {
            let val = match index {
                Index::Input(i) => assignment[*i],
                Index::Aux(i) => assignment[num_inputs + i],
            };

            if coeff.is_one() {
                val
            } else {
                val.mul(coeff)
            }
        })
        .reduce(|| E::Fr::zero(), |sum, val| sum + &val)
}
/*A R1CS consisting of n constraints in m variables

    (a_{i,0} + Sum_j x_j a_{i,j}) * (b_{i,0} + Sum_{j=1}^m x_j b_{i,j}) =
               c (c_{i,0} + Sum_{j=1}^m x_j c_{i,j}),

i=1..n, translated into the QAP is

    (a_0(Z) + Sum_j x_j a_j(Z)) * (b_0(Z) + Sum_j x_j b_j(Z)) =
                (c_0(Z) + Sum_j x_j c_j(Z))

in F[Z]/(Z^N-1). The polynomials a_j(Z),b_j(Z),c_j(Z) correspond to the column vectors
a_{*,j}, b_{*,j}, c_{*,j} regarded as functions on the FFT domain H = {z^N - 1} = {z_1,..z_N},

    a_j(Z) = Sum_{i=1}^N a_{i,j} L_i (Z),

and similarly to b, c (here, L_i is the Lagrange polynomial at z_i).
*/
impl R1CStoQAP {
    #[inline]
    ///Given a KeypairAssembly, i.e. the constraint-wise description of the R1CS, this function
    ///returns
    ///     - the column polynomials a_j(Z), b_j(Z), c_j(Z) evaluated at the secret point t,
    ///       as vectors
    ///         a = (a_j(t))_{j=0}^m, b=(b_j(t))_{j=0}^m, c=(c_j(t))_{j=0}^m,
    ///       as well as
    ///     - the vanishing polynomial of the FFT domain H evaluated at t,
    ///         zt= v_H(t),
    ///       and
    ///     - the number qap_num_ variables = m of QAP variables, as well as
    ///     - the domain size |H|.
    pub(crate) fn instance_map_with_evaluation<E: PairingEngine>(
        assembly: &KeypairAssembly<E>,
        t: &E::Fr,
    ) -> Result<(Vec<E::Fr>, Vec<E::Fr>, Vec<E::Fr>, E::Fr, usize, usize), SynthesisError> {
        let domain_size = assembly.num_constraints + (assembly.num_inputs - 1) + 1;
        let domain = get_best_evaluation_domain::<E::Fr>(domain_size)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let domain_size = domain.size();

        let zt = domain.evaluate_vanishing_polynomial(*t);

        //Evaluate all Lagrange polynomials L_i(t) for "i in H".
        let coefficients_time = start_timer!(|| "Evaluate Lagrange coefficients");
        let u = domain.evaluate_all_lagrange_coefficients(*t);
        end_timer!(coefficients_time);

        //one variable in the R1CS is always for the constants
        let qap_num_variables = (assembly.num_inputs - 1) + assembly.num_aux;

        let mut a = vec![E::Fr::zero(); qap_num_variables + 1];
        let mut b = vec![E::Fr::zero(); qap_num_variables + 1];
        let mut c = vec![E::Fr::zero(); qap_num_variables + 1];

        //The points i= n, n+1, .. ,n +l-1 correspond to the copy-paste constraints to
        //load x_1,..,x_l with the public input
        for i in 0..assembly.num_inputs {
            a[i] = u[assembly.num_constraints + i];
        }

        //constraint-wise, i.e. row-wise scanning of m_{i,j}!=0
        //and incrementing column m_j(t), m=a,b,c.
        for i in 0..assembly.num_constraints {
            for &(ref coeff, index) in assembly.at[i].iter() {
                //convert R1CS index into the corresponding number from [0,...,m]
                let index = match index {
                    Index::Input(j) => j, //if index is an input variable, return it as it is
                    Index::Aux(j) => assembly.num_inputs + j, //if index is a private variable, shift it
                };
                //update the column sum corresponding to the variable.
                a[index] += &(u[i] * coeff);
            }
            for &(ref coeff, index) in assembly.bt[i].iter() {
                let index = match index {
                    Index::Input(i) => i,
                    Index::Aux(i) => assembly.num_inputs + i,
                };

                b[index] += &(u[i] * coeff);
            }
            for &(ref coeff, index) in assembly.ct[i].iter() {
                let index = match index {
                    Index::Input(i) => i,
                    Index::Aux(i) => assembly.num_inputs + i,
                };

                c[index] += &(u[i] * coeff);
            }
        }

        Ok((a, b, c, zt, qap_num_variables, domain_size))
    }

    //computes the coefficients of the quotient polynomial
    // h(Z) = (a(Z)*b(Z)-c(Z))/v_H(Z)
    //from the witness assignments of the circuit.
    //We have deg(h(Z))= deg(a(Z))+deg(b(Z)) - deg v_H(Z) <= |H|-1 + |H|-1 - |H|
    //  = |H|-2.
    #[inline]
    pub(crate) fn witness_map<E: PairingEngine>(
        prover: &ProvingAssignment<E>,
    ) -> Result<Vec<E::Fr>, SynthesisError> {
        let zero = E::Fr::zero();

        let num_inputs = prover.input_assignment.len();
        let num_constraints = prover.num_constraints();

        let full_input_assignment =
            [&prover.input_assignment[..], &prover.aux_assignment[..]].concat();

        // including the copy-paste constraints for the public inputs, the full
        // number of constraints equals 'num_constraints + num_inputs'.
        let domain = get_best_evaluation_domain::<E::Fr>(num_constraints + num_inputs)
            .ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let domain_size = domain.size();
        /*
                println!("num_constraints: {}", num_constraints);
                println!("num_inputs: {}", num_inputs);
                println!("Domain H size: {}", domain_size);
        */
        let mut a = vec![zero; domain_size];
        let mut b = vec![zero; domain_size];
        //compute the evaluations of a(Z), b(Z) on H
        //for each constraint, i=1..n, compute m_i= sum_j xj mij, m=a,b,c.
        a[..num_constraints]
            .par_iter_mut()
            .zip(b[..num_constraints].par_iter_mut())
            .zip((&prover.at).par_iter())
            .zip((&prover.bt).par_iter())
            .for_each(|(((a, b), at_i), bt_i)| {
                *a = evaluate_constraint::<E>(&at_i, &full_input_assignment, num_inputs);
                *b = evaluate_constraint::<E>(&bt_i, &full_input_assignment, num_inputs);
            });
        //the further a_i, i=n+1,..,n+l are for the public inputs
        for i in 0..num_inputs {
            a[num_constraints + i] = full_input_assignment[i];
        }
        //compute the coefficients of a(Z) and b(Z)
        domain.ifft_in_place(&mut a);
        domain.ifft_in_place(&mut b);
        //and their values on the double-sized FFT domain
        //K = H v coset(H)
        //actually, here we compute ONLY their evals over the
        //coset of H (because we now that a(Z)*b(Z)-c(Z) will be zero on H).
        domain.coset_fft_in_place(&mut a);
        domain.coset_fft_in_place(&mut b);

        //compute the product of the evals of a(Z)*b(Z) on the coset of H
        let mut ab = domain.mul_polynomials_in_evaluation_domain(&a, &b)?;
        drop(a);
        drop(b);

        //compute values of c(Z) on H
        let mut c = vec![zero; domain_size];
        c[..num_constraints]
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, c)| {
                *c = evaluate_constraint::<E>(&prover.ct[i], &full_input_assignment, num_inputs);
            });

        //extrapolate c(Z) from H to the coset of H and
        //compute a(Z)*b(Z)-c(Z) on this coset of H
        domain.ifft_in_place(&mut c);
        domain.coset_fft_in_place(&mut c);

        ab.par_iter_mut()
            .zip(c)
            .for_each(|(ab_i, c_i)| *ab_i -= &c_i);

        // compute quotient polynomial (a(Z)*b(Z)-c(Z))/v_H(Z)
        // from the coset evaluations
        domain.divide_by_vanishing_poly_on_coset_in_place(&mut ab);
        domain.coset_ifft_in_place(&mut ab);

        // we drop the leading coefficient, as deg(h(Z)) = n-2.
        assert!(ab.pop().unwrap() == zero);
        Ok(ab)
    }
}
