//! This module contains an `EvaluationDomain` abstraction for
//! performing various kinds of polynomial arithmetic on top of
//! the scalar field.
//!
//! In pairing-based SNARKs we need to calculate a quotient
//! polynomial over a target polynomial with roots at distinct
//! points associated with each constraint of the constraint
//! system.

#[doc(hidden)]
pub mod domain_selector;
pub use self::domain_selector::*;

#[doc(hidden)]
pub mod basic_radix_2_domain;
pub use self::basic_radix_2_domain::*;

#[doc(hidden)]
pub mod mixed_radix_2_domain;
pub use self::mixed_radix_2_domain::*;

#[cfg(test)]
mod test;

use crate::{
    SparsePolynomial,
    multicore::Worker,
};
use crate::PrimeField;
use rayon::prelude::*;
use std::hash::Hash;
use std::fmt::Debug;
use rand::Rng;

/// Defines a domain over which finite field (I)FFTs can be performed.
pub trait EvaluationDomainImpl<F: PrimeField>: Sized + Copy + Clone + Eq + PartialEq + Debug + Hash + Default
{
    /// Creates an EvaluationDomain of size at least `num_coeffs`.
    fn new(num_coeffs: usize) -> Option<Self>;

    /// Return the size of a domain that is large enough for evaluations of a polynomial
    /// having `num_coeffs` coefficients.
    fn compute_size_of_domain(num_coeffs: usize) -> Option<usize>;

    /// Returns the size of the domain
    fn size(&self) -> usize;

    /// Interprets size as a field element F and returns its inverse in the field
    fn size_inv(&self) -> F;

    /// Returns the generator of the multiplicative subgroup over which FFT is performed
    fn group_gen(&self) -> F;

    /// Sample an element that is *not* in the domain.
    fn sample_element_outside_domain<R: Rng>(&self, rng: &mut R) -> F {
        let mut t = F::rand(rng);
        while self.evaluate_vanishing_polynomial(t).is_zero() {
            t = F::rand(rng);
        }
        t
    }

    /// Compute a FFT.
    fn fft(&self, coeffs: &[F]) -> Vec<F> {
        let mut coeffs = coeffs.to_vec();
        self.fft_in_place(&mut coeffs);
        coeffs
    }

    /// Compute a FFT, modifying the vector in place.
    fn fft_in_place(&self, coeffs: &mut Vec<F>);

    /// Compute a FFT over a coset of the domain.
    fn coset_fft(&self, coeffs: &[F]) -> Vec<F> {
        let mut coeffs = coeffs.to_vec();
        self.coset_fft_in_place(&mut coeffs);
        coeffs
    }

    /// Compute a FFT over a coset of the domain, modifying the input vector
    /// in place.
    fn coset_fft_in_place(&self, coeffs: &mut Vec<F>);

    /// Compute a IFFT.
    fn ifft(&self, evals: &[F]) -> Vec<F> {
        let mut evals = evals.to_vec();
        self.ifft_in_place(&mut evals);
        evals
    }

    /// Compute a IFFT, modifying the vector in place.
    fn ifft_in_place(&self, evals: &mut Vec<F>);

    /// Compute a IFFT over a coset of the domain.
    fn coset_ifft(&self, evals: &[F]) -> Vec<F> {
        let mut evals = evals.to_vec();
        self.coset_ifft_in_place(&mut evals);
        evals
    }

    /// Compute a IFFT over a coset of the domain, modifying the input vector in place.
    fn coset_ifft_in_place(&self, evals: &mut Vec<F>);

    /// Return the sparse vanishing polynomial.
    fn vanishing_polynomial(&self) -> SparsePolynomial<F> {
        let coeffs = vec![(0, -F::one()), (self.size(), F::one())];
        SparsePolynomial::from_coefficients_vec(coeffs)
    }

    /// This evaluates the vanishing polynomial for this domain at tau.
    /// For multiplicative subgroups, this polynomial is `z(X) = X^self.size - 1`.
    fn evaluate_vanishing_polynomial(&self, tau: F) -> F {
        tau.pow(&[self.size() as u64]) - &F::one()
    }

    /// The target polynomial is the zero polynomial in our
    /// evaluation domain, so we must perform division over
    /// a coset.
    fn divide_by_vanishing_poly_on_coset_in_place(&self, evals: &mut [F]) {
        let i = self
            .evaluate_vanishing_polynomial(F::multiplicative_generator())
            .inverse()
            .unwrap();

        Worker::new().scope(evals.len(), |scope, chunk| {
            for evals in evals.chunks_mut(chunk) {
                scope.spawn(move |_| evals.iter_mut().for_each(|eval| *eval *= &i));
            }
        });
    }

    /// Perform O(n) multiplication of two polynomials that are presented by their
    /// evaluations in the domain.
    /// Returns the evaluations of the product over the domain.
    #[must_use]
    fn mul_polynomials_in_evaluation_domain(&self, self_evals: &[F], other_evals: &[F]) -> Vec<F> {
        assert_eq!(self_evals.len(), other_evals.len());
        let mut result = self_evals.to_vec();
        result.par_iter_mut().zip(other_evals).for_each(|(a,b)| *a *= b);
        result
    }

    /// Evaluate all the lagrange polynomials defined by this domain at the point
    /// `tau`.
    fn evaluate_all_lagrange_coefficients(&self, tau: F) -> Vec<F> {
        // Evaluate all Lagrange polynomials
        let size = self.size();
        let t_size = tau.pow(&[size as u64]);
        let one = F::one();
        if t_size.is_one() {
            let mut u = vec![F::zero(); size];
            let mut omega_i = one;
            for i in 0..size {
                if omega_i == tau {
                    u[i] = one;
                    break;
                }
                omega_i *= &self.group_gen();
            }
            u
        } else {
            use crate::fields::batch_inversion;

            let mut l = (t_size - &one) * &self.size_inv();
            let mut r = one;
            let mut u = vec![F::zero(); size];
            let mut ls = vec![F::zero(); size];
            for i in 0..size {
                u[i] = tau - &r;
                ls[i] = l;
                l *= &self.group_gen();
                r *= &self.group_gen();
            }

            batch_inversion(u.as_mut_slice());
            u.par_iter_mut().zip(ls).for_each(|(tau_minus_r, l)| {
                *tau_minus_r = l * tau_minus_r;
            });
            u
        }
    }

    #[doc(hidden)]
    fn distribute_powers(coeffs: &mut Vec<F>, g: F) {
        Worker::new().scope(coeffs.len(), |scope, chunk| {
            for (i, v) in coeffs.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut u = g.pow(&[(i * chunk) as u64]);
                    for v in v.iter_mut() {
                        *v *= &u;
                        u *= &g;
                    }
                });
            }
        });
    }

    /// Return an iterator over the elements of the domain.
    fn elements(&self) -> Elements<F> {
        Elements {
            cur_elem: F::one(),
            cur_pow: 0,
            size: self.size() as u64,
            group_gen: self.group_gen()
        }
    }
}

/// An iterator over the elements of the domain.
pub struct Elements<F: PrimeField> {
    cur_elem:   F,
    cur_pow:    u64,
    size:       u64,
    group_gen:  F,
}

impl<F: PrimeField> Iterator for Elements<F> {
    type Item = F;
    fn next(&mut self) -> Option<F> {
        if self.cur_pow == self.size {
            None
        } else {
            let cur_elem = self.cur_elem;
            self.cur_elem *= &self.group_gen;
            self.cur_pow += 1;
            Some(cur_elem)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{EvaluationDomain, EvaluationDomainImpl};
    use crate::Field;
    use crate::fields::mnt6753::fr::Fr;
    use rand::{Rng, thread_rng};

    #[test]
    fn vanishing_polynomial_evaluation() {
        let rng = &mut thread_rng();
        for coeffs in 0..18 {
            let size = 1 << coeffs;
            let domain = EvaluationDomain::<Fr>::new(size).unwrap();
            let z = domain.vanishing_polynomial();
            for _ in 0..100 {
                let point = rng.gen();
                assert_eq!(z.evaluate(point), domain.evaluate_vanishing_polynomial(point))
            }
        }
    }

    #[test]
    fn vanishing_polynomial_vanishes_on_domain() {
        for coeffs in 0..18 {
            let size = 1 << coeffs;
            let domain = EvaluationDomain::<Fr>::new(size).unwrap();
            let z = domain.vanishing_polynomial();
            for point in domain.elements() {
                assert!(z.evaluate(point).is_zero())
            }
        }
    }

    #[test]
    fn size_of_elements() {
        for coeffs in 1..18 {
            let size = 1 << coeffs;
            let domain = EvaluationDomain::<Fr>::new(size).unwrap();
            let domain_size = domain.size();
            assert_eq!(domain_size, domain.elements().collect::<Vec<_>>().len());
        }
    }

    #[test]
    fn elements_contents() {
        for coeffs in 1..18 {
            let size = 1 << coeffs;
            let domain = EvaluationDomain::<Fr>::new(size).unwrap();
            for (i, element) in domain.elements().enumerate() {
                assert_eq!(element, domain.group_gen().pow([i as u64]));
            }
        }
    }
}