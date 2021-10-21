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
    Error,
};
use crate::PrimeField;
use rayon::prelude::*;
//use std::hash::Hash;
use std::fmt::Debug;
use std::any::Any;
use rand::Rng;

/// Defines a domain over which finite field (I)FFTs can be performed.
pub trait EvaluationDomain<F: PrimeField>: Debug + Send + Sync
{
    /// Returns the size of the domain
    fn size(&self) -> usize;

    fn size_as_field_element(&self) -> F {
        F::from_repr(F::BigInt::from(self.size() as u64))
    }

    /// Interprets size as a field element F and returns its inverse in the field
    fn size_inv(&self) -> F;

    /// Returns the generator of the multiplicative subgroup over which FFT is performed
    fn group_gen(&self) -> F;

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
    fn mul_polynomials_in_evaluation_domain(&self, self_evals: &[F], other_evals: &[F]) -> Result<Vec<F>, Error> {
        if self_evals.len() != other_evals.len() {
            Err(format!("Evals sizes are not same"))?
        }
        let mut result = self_evals.to_vec();
        result.par_iter_mut().zip(other_evals).for_each(|(a,b)| *a *= b);
        Ok(result)
    }

    /// Given an arbitrary field element `tau`, compute the Lagrange kernel 
    ///     L(z,tau) = 1/n * z * (1 - tau^n)/(z -tau).
    /// The Lagrange kernel is useful when one needs to evaluate many polynomials given in 
    /// Lagrange representation at that given point.
    /// This implementation works also if `tau` is selected from the domain.
    fn evaluate_all_lagrange_coefficients(&self, tau: F) -> Vec<F> {
        // Evaluate all Lagrange polynomials
        let size = self.size();
        let t_size = tau.pow(&[size as u64]);
        let one = F::one();
        if t_size.is_one() {
            // if tau is from the domain itself, then L(z,tau) is "trivial"
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
            // we compute L(z,tau) = 1/n * z * (tau^n - 1)/(tau - z)
            // using batch inversion for (tau - z), z over H.
            use crate::fields::batch_inversion;

            let mut l = (t_size - &one) * &self.size_inv();
            let mut r = one;
            let mut u = vec![F::zero(); size];
            let mut ls = vec![F::zero(); size];
            // u[i] = tau - z  at z = g^i,
            // ls[i] = (tau^n - 1)/n * g^i, 
            for i in 0..size {
                u[i] = tau - &r;
                ls[i] = l;
                l *= &self.group_gen();
                r *= &self.group_gen();
            }

            batch_inversion(u.as_mut_slice());
            // We compute L(z,tau) = u[i]*ls[i]. 
            u.par_iter_mut().zip(ls).for_each(|(tau_minus_r, l)| {
                *tau_minus_r *= l;
            });
            u
        }
    }

    /// Given an index which assumes the first elements of this domain are the elements of
    /// another (sub)domain with size size_s, this returns the actual index into this domain.
    fn reindex_by_subdomain(&self, other_size: usize, index: usize) -> Result<usize, Error> {
        if self.size() < other_size {
            Err(format!(
                "'self' domain size must be bigger than 'other' domain size. 'self' size: {}, 'other' size: {}",
                self.size(),
                other_size
            ))?
        }
        if other_size == 0 {
            Err(format!("'other' size must be bigger than 0"))?
        }
        // Let this subgroup be G, and the subgroup we're re-indexing by be S.
        // Since its a subgroup, the 0th element of S is at index 0 in G, the first element of S is at
        // index |G|/|S|, the second at 2*|G|/|S|, etc.
        // Thus for an index i that corresponds to S, the index in G is i*|G|/|S|
        let period = self.size() / other_size;
        if index < other_size {
            Ok(index * period)
        } else {
            // Let i now be the index of this element in G \ S
            // Let x be the number of elements in G \ S, for every element in S. Then x = (|G|/|S| - 1).
            // At index i in G \ S, the number of elements in S that appear before the index in G to which
            // i corresponds to, is floor(i / x) + 1.
            // The +1 is because index 0 of G is S_0, so the position is offset by at least one.
            // The floor(i / x) term is because after x elements in G \ S, there is one more element from S
            // that will have appeared in G.
            if period <= 1 {
                Err(format!("'period' must be bigger than 1"))?
            }
            let i = index - other_size;
            let x = period - 1;
            Ok(i + (i / x) + 1)
        }
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

    // Support to PartialEq to make this trait a trait object
    fn eq(&self, other: & dyn EvaluationDomain<F>) -> bool;

    fn as_any(&self) -> & dyn Any;

    // Support to Clone to make this trait a trait object
    fn clone_and_box(&self) -> Box<dyn EvaluationDomain<F>>;
}

impl<'a, 'b, F: PrimeField> PartialEq<dyn EvaluationDomain<F>+'b> for dyn EvaluationDomain<F>+'a {
    fn eq(&self, other: &(dyn EvaluationDomain<F>+'b)) -> bool {
        EvaluationDomain::<F>::eq(self, other)
    }
}

impl<F: PrimeField> Clone for Box<dyn EvaluationDomain<F>>
{
    fn clone(&self) -> Box<dyn EvaluationDomain<F>> {
        self.clone_and_box()
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

pub fn sample_element_outside_domain<
    F: PrimeField,
    R: Rng
>(domain: &Box<dyn EvaluationDomain<F>>, rng: &mut R) -> F {
    let mut t = F::rand(rng);
    while domain.evaluate_vanishing_polynomial(t).is_zero() {
        t = F::rand(rng);
    }
    t
}

#[cfg(test)]
mod tests {
    use crate::get_best_evaluation_domain;
    use crate::Field;
    use crate::fields::bls12_381::fr::Fr;
    use rand::{Rng, thread_rng};

    #[test]
    fn vanishing_polynomial_evaluation() {
        let rng = &mut thread_rng();
        for coeffs in 0..18 {
            let size = 1 << coeffs;
            let domain = get_best_evaluation_domain::<Fr>(size).unwrap();
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
            let domain = get_best_evaluation_domain::<Fr>(size).unwrap();
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
            let domain = get_best_evaluation_domain::<Fr>(size).unwrap();
            let domain_size = domain.size();
            assert_eq!(domain_size, domain.elements().collect::<Vec<_>>().len());
        }
    }

    #[test]
    fn elements_contents() {
        for coeffs in 1..18 {
            let size = 1 << coeffs;
            let domain = get_best_evaluation_domain::<Fr>(size).unwrap();
            for (i, element) in domain.elements().enumerate() {
                assert_eq!(element, domain.group_gen().pow([i as u64]));
            }
        }
    }
}
