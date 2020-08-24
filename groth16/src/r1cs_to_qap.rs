use algebra_core::{One, PairingEngine, Zero};
use ff_fft::{cfg_iter, cfg_iter_mut, EvaluationDomain};

use crate::Vec;
use core::ops::{AddAssign, Deref};
use r1cs_core::{ConstraintSystemRef, SynthesisError};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[inline]
fn evaluate_constraint<'a, LHS, RHS, R>(terms: &'a [(LHS, usize)], assignment: &'a [RHS]) -> R
where
    LHS: One + Send + Sync + PartialEq,
    RHS: Send + Sync + core::ops::Mul<&'a LHS, Output = RHS> + Copy,
    R: Zero + Send + Sync + AddAssign<RHS> + core::iter::Sum,
{
    // Need to wrap in a closure when using Rayon
    #[cfg(feature = "parallel")]
    let zero = || R::zero();
    #[cfg(not(feature = "parallel"))]
    let zero = R::zero();

    let res = cfg_iter!(terms).fold(zero, |mut sum, (coeff, index)| {
        let val = &assignment[*index];

        if coeff.is_one() {
            sum += *val;
        } else {
            sum += val.mul(coeff);
        }

        sum
    });

    // Need to explicitly call `.sum()` when using Rayon
    #[cfg(feature = "parallel")]
    return res.sum();
    #[cfg(not(feature = "parallel"))]
    return res;
}

pub(crate) struct R1CStoQAP;

impl R1CStoQAP {
    #[inline]
    pub(crate) fn instance_map_with_evaluation<E: PairingEngine, D: EvaluationDomain<E::Fr>>(
        cs: ConstraintSystemRef<E::Fr>,
        t: &E::Fr,
    ) -> Result<(Vec<E::Fr>, Vec<E::Fr>, Vec<E::Fr>, E::Fr, usize, usize), SynthesisError> {
        let matrices = cs.to_matrices().unwrap();
        let domain_size = cs.num_constraints() + cs.num_instance_variables();
        let domain = D::new(domain_size).ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let domain_size = domain.size();

        let zt = domain.evaluate_vanishing_polynomial(*t);

        // Evaluate all Lagrange polynomials
        let coefficients_time = start_timer!(|| "Evaluate Lagrange coefficients");
        let u = domain.evaluate_all_lagrange_coefficients(*t);
        end_timer!(coefficients_time);

        let qap_num_variables = (cs.num_instance_variables() - 1) + cs.num_witness_variables();

        let mut a = vec![E::Fr::zero(); qap_num_variables + 1];
        let mut b = vec![E::Fr::zero(); qap_num_variables + 1];
        let mut c = vec![E::Fr::zero(); qap_num_variables + 1];

        for i in 0..cs.num_instance_variables() {
            a[i] = u[cs.num_constraints() + i];
        }

        for i in 0..cs.num_constraints() {
            for &(ref coeff, index) in &matrices.a[i] {
                a[index] += &(u[i] * coeff);
            }
            for &(ref coeff, index) in &matrices.b[i] {
                b[index] += &(u[i] * coeff);
            }
            for &(ref coeff, index) in &matrices.c[i] {
                c[index] += &(u[i] * coeff);
            }
        }

        Ok((a, b, c, zt, qap_num_variables, domain_size))
    }

    #[inline]
    pub(crate) fn witness_map<E: PairingEngine, D: EvaluationDomain<E::Fr>>(
        prover: ConstraintSystemRef<E::Fr>,
    ) -> Result<Vec<E::Fr>, SynthesisError> {
        let matrices = prover.to_matrices().unwrap();
        let zero = E::Fr::zero();
        let num_inputs = prover.num_instance_variables();
        let num_constraints = prover.num_constraints();
        let cs = prover.borrow().unwrap();
        let prover = cs.deref();

        let full_assignment = [
            prover.instance_assignment.as_slice(),
            prover.witness_assignment.as_slice(),
        ]
        .concat();

        let domain =
            D::new(num_constraints + num_inputs).ok_or(SynthesisError::PolynomialDegreeTooLarge)?;
        let domain_size = domain.size();

        let mut a = vec![zero; domain_size];
        let mut b = vec![zero; domain_size];

        cfg_iter_mut!(a[..num_constraints])
            .zip(cfg_iter_mut!(b[..num_constraints]))
            .zip(cfg_iter!(&matrices.a))
            .zip(cfg_iter!(&matrices.b))
            .for_each(|(((a, b), at_i), bt_i)| {
                *a = evaluate_constraint(&at_i, &full_assignment);
                *b = evaluate_constraint(&bt_i, &full_assignment);
            });

        for i in 0..num_inputs {
            a[num_constraints + i] = full_assignment[i];
        }

        domain.ifft_in_place(&mut a);
        domain.ifft_in_place(&mut b);

        domain.coset_fft_in_place(&mut a);
        domain.coset_fft_in_place(&mut b);

        let mut ab = domain.mul_polynomials_in_evaluation_domain(&a, &b);
        drop(a);
        drop(b);

        let mut c = vec![zero; domain_size];
        cfg_iter_mut!(c[..prover.num_constraints])
            .enumerate()
            .for_each(|(i, c)| {
                *c = evaluate_constraint(&matrices.c[i], &full_assignment);
            });

        domain.ifft_in_place(&mut c);
        domain.coset_fft_in_place(&mut c);

        cfg_iter_mut!(ab)
            .zip(c)
            .for_each(|(ab_i, c_i)| *ab_i -= &c_i);

        domain.divide_by_vanishing_poly_on_coset_in_place(&mut ab);
        domain.coset_ifft_in_place(&mut ab);

        Ok(ab)
    }
}
