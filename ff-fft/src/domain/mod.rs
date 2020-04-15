//! This module contains an `EvaluationDomain` abstraction for
//! performing various kinds of polynomial arithmetic on top of
//! the scalar field.
//!
//! In pairing-based SNARKs like GM17, we need to calculate
//! a quotient polynomial over a target polynomial with roots
//! at distinct points associated with each constraint of the
//! constraint system. In order to be efficient, we try to choose these
//! roots to be the powers of a 2^n root of unity in the field.
//! This allows us to perform polynomial operations in O(n)
//! by performing an O(n log n) FFT over such a domain.
//!
//! If the 2-adicity of the field is too small, but a small subgroup
//! over the field is defined, we can try to build a mixed-radix evaluation
//! domain.

use crate::Vec;
use algebra_core::FftField;
use core::{fmt, hash};
use rand::Rng;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub mod general;
pub mod mixed_radix;
pub mod radix2;
pub(crate) mod utils;

pub use general::GeneralEvaluationDomain;
pub use mixed_radix::MixedRadixEvaluationDomain;
pub use radix2::Radix2EvaluationDomain;

/// Defines a domain over which finite field (I)FFTs can be performed. Works
/// only for fields that have a large multiplicative subgroup of size that is
/// a power-of-2.
pub trait EvaluationDomain<F: FftField>:
    Copy + Clone + hash::Hash + Eq + PartialEq + fmt::Debug
{
    /// The type of the elements iterator.
    type Elements: Iterator<Item = F> + Sized;

    /// Sample an element that is *not* in the domain.
    fn sample_element_outside_domain<R: Rng>(&self, rng: &mut R) -> F {
        let mut t = F::rand(rng);
        while self.evaluate_vanishing_polynomial(t).is_zero() {
            t = F::rand(rng);
        }
        t
    }

    /// Construct a domain that is large enough for evaluations of a polynomial
    /// having `num_coeffs` coefficients.
    fn new(num_coeffs: usize) -> Option<Self>;

    /// Return the size of a domain that is large enough for evaluations of a
    /// polynomial having `num_coeffs` coefficients.
    fn compute_size_of_domain(num_coeffs: usize) -> Option<usize>;

    /// Return the size of `self`.
    fn size(&self) -> usize;

    /// Compute a FFT.
    #[inline]
    fn fft<T: DomainCoeff<F>>(&self, coeffs: &[T]) -> Vec<T> {
        let mut coeffs = coeffs.to_vec();
        self.fft_in_place(&mut coeffs);
        coeffs
    }

    /// Compute a FFT, modifying the vector in place.
    fn fft_in_place<T: DomainCoeff<F>>(&self, coeffs: &mut Vec<T>);

    /// Compute a IFFT.
    #[inline]
    fn ifft<T: DomainCoeff<F>>(&self, evals: &[T]) -> Vec<T> {
        let mut evals = evals.to_vec();
        self.ifft_in_place(&mut evals);
        evals
    }

    /// Compute a IFFT, modifying the vector in place.
    fn ifft_in_place<T: DomainCoeff<F>>(&self, evals: &mut Vec<T>);

    /// Distribute the powers of `g` among the `coeffs`.
    fn distribute_powers<T: DomainCoeff<F>>(coeffs: &mut [T], g: F) {
        let mut pow = F::one();
        coeffs.iter_mut().for_each(|c| {
            *c *= pow;
            pow *= &g
        })
    }

    /// Compute a FFT over a coset of the domain.
    #[inline]
    fn coset_fft<T: DomainCoeff<F>>(&self, coeffs: &[T]) -> Vec<T> {
        let mut coeffs = coeffs.to_vec();
        self.coset_fft_in_place(&mut coeffs);
        coeffs
    }

    /// Compute a FFT over a coset of the domain, modifying the input vector
    /// in place.
    #[inline]
    fn coset_fft_in_place<T: DomainCoeff<F>>(&self, coeffs: &mut Vec<T>) {
        Self::distribute_powers(coeffs, F::multiplicative_generator());
        self.fft_in_place(coeffs);
    }

    /// Compute a IFFT over a coset of the domain.
    #[inline]
    fn coset_ifft<T: DomainCoeff<F>>(&self, evals: &[T]) -> Vec<T> {
        let mut evals = evals.to_vec();
        self.coset_ifft_in_place(&mut evals);
        evals
    }

    /// Compute a IFFT over a coset of the domain, modifying the input vector in
    /// place.
    #[inline]
    fn coset_ifft_in_place<T: DomainCoeff<F>>(&self, evals: &mut Vec<T>) {
        self.ifft_in_place(evals);
        Self::distribute_powers(evals, F::multiplicative_generator().inverse().unwrap());
    }

    /// Evaluate all the lagrange polynomials defined by this domain at the
    /// point `tau`.
    fn evaluate_all_lagrange_coefficients(&self, tau: F) -> Vec<F>;

    /// Return the sparse vanishing polynomial.
    fn vanishing_polynomial(&self) -> crate::SparsePolynomial<F>;

    /// This evaluates the vanishing polynomial for this domain at tau.
    fn evaluate_vanishing_polynomial(&self, tau: F) -> F;

    /// Return an iterator over the elements of the domain.
    fn elements(&self) -> Self::Elements;

    /// The target polynomial is the zero polynomial in our
    /// evaluation domain, so we must perform division over
    /// a coset.
    fn divide_by_vanishing_poly_on_coset_in_place(&self, evals: &mut [F]) {
        let i = self
            .evaluate_vanishing_polynomial(F::multiplicative_generator())
            .inverse()
            .unwrap();

        cfg_iter_mut!(evals).for_each(|eval| *eval *= &i);
    }

    /// Given an index which assumes the first elements of this domain are the
    /// elements of another (sub)domain with size size_s,
    /// this returns the actual index into this domain.
    fn reindex_by_subdomain(&self, other: Self, index: usize) -> usize {
        assert!(self.size() >= other.size());
        // Let this subgroup be G, and the subgroup we're re-indexing by be S.
        // Since its a subgroup, the 0th element of S is at index 0 in G, the first
        // element of S is at index |G|/|S|, the second at 2*|G|/|S|, etc.
        // Thus for an index i that corresponds to S, the index in G is i*|G|/|S|
        let period = self.size() / other.size();
        if index < other.size() {
            index * period
        } else {
            // Let i now be the index of this element in G \ S
            // Let x be the number of elements in G \ S, for every element in S. Then x =
            // (|G|/|S| - 1). At index i in G \ S, the number of elements in S
            // that appear before the index in G to which i corresponds to, is
            // floor(i / x) + 1. The +1 is because index 0 of G is S_0, so the
            // position is offset by at least one. The floor(i / x) term is
            // because after x elements in G \ S, there is one more element from S
            // that will have appeared in G.
            let i = index - other.size();
            let x = period - 1;
            i + (i / x) + 1
        }
    }

    /// Perform O(n) multiplication of two polynomials that are presented by
    /// their evaluations in the domain.
    /// Returns the evaluations of the product over the domain.
    ///
    /// Assumes that the domain is large enough to allow for successful
    /// interpolation after multiplication.
    #[must_use]
    fn mul_polynomials_in_evaluation_domain(&self, self_evals: &[F], other_evals: &[F]) -> Vec<F> {
        assert_eq!(self_evals.len(), other_evals.len());
        let mut result = self_evals.to_vec();

        cfg_iter_mut!(result)
            .zip(other_evals)
            .for_each(|(a, b)| *a *= b);

        result
    }
}

/// Types that can be FFT-ed must implement this trait.
pub trait DomainCoeff<F: FftField>:
    Copy
    + Send
    + Sync
    + core::ops::AddAssign
    + core::ops::SubAssign
    + algebra_core::Zero
    + core::ops::MulAssign<F>
{
}

impl<T, F> DomainCoeff<F> for T
where
    F: FftField,
    T: Copy
        + Send
        + Sync
        + core::ops::AddAssign
        + core::ops::SubAssign
        + algebra_core::Zero
        + core::ops::MulAssign<F>,
{
}
