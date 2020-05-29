//! This module contains a `GeneralEvaluationDomain` for
//! performing various kinds of polynomial arithmetic on top of
//! a FFT-friendly finite field.
//!
//! It is a wrapper around specific implementations of `EvaluationDomain` that
//! automatically chooses the most efficient implementation
//! depending on the number of coefficients and the two-adicity of the prime.

pub use crate::domain::utils::Elements;
use crate::domain::{
    DomainCoeff, EvaluationDomain, MixedRadixEvaluationDomain, Radix2EvaluationDomain,
};
use crate::Vec;
use algebra_core::{FftField, FftParameters};

/// Defines a domain over which finite field (I)FFTs can be performed.
/// Generally tries to build a radix-2 domain and falls back to a mixed-radix
/// domain if the radix-2 multiplicative subgroup is too small.
#[derive(Copy, Clone, Hash, Eq, PartialEq, Debug)]
pub enum GeneralEvaluationDomain<F: FftField> {
    /// Radix-2 domain
    Radix2(Radix2EvaluationDomain<F>),
    /// Mixed-radix domain
    MixedRadix(MixedRadixEvaluationDomain<F>),
}

impl<F: FftField> EvaluationDomain<F> for GeneralEvaluationDomain<F> {
    type Elements = GeneralElements<F>;

    /// Construct a domain that is large enough for evaluations of a polynomial
    /// having `num_coeffs` coefficients.
    ///
    /// If the field specifies a small subgroup for a mixed-radix FFT and
    /// the radix-2 FFT cannot be constructed, this method tries
    /// constructing a mixed-radix FFT instead.
    fn new(num_coeffs: usize) -> Option<Self> {
        let domain = Radix2EvaluationDomain::new(num_coeffs);
        if let Some(domain) = domain {
            return Some(GeneralEvaluationDomain::Radix2(domain));
        }

        if F::FftParams::SMALL_SUBGROUP_BASE.is_some() {
            return Some(GeneralEvaluationDomain::MixedRadix(
                MixedRadixEvaluationDomain::new(num_coeffs)?,
            ));
        }

        None
    }

    fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        let domain_size = Radix2EvaluationDomain::<F>::compute_size_of_domain(num_coeffs);
        if let Some(domain_size) = domain_size {
            return Some(domain_size);
        }

        if F::FftParams::SMALL_SUBGROUP_BASE.is_some() {
            return Some(MixedRadixEvaluationDomain::<F>::compute_size_of_domain(
                num_coeffs,
            )?);
        }

        None
    }

    #[inline]
    fn size(&self) -> usize {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.size(),
            GeneralEvaluationDomain::MixedRadix(domain) => domain.size(),
        }
    }

    #[inline]
    fn fft_in_place<T: DomainCoeff<F>>(&self, coeffs: &mut Vec<T>) {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.fft_in_place(coeffs),
            GeneralEvaluationDomain::MixedRadix(domain) => domain.fft_in_place(coeffs),
        }
    }

    #[inline]
    fn ifft_in_place<T: DomainCoeff<F>>(&self, evals: &mut Vec<T>) {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.ifft_in_place(evals),
            GeneralEvaluationDomain::MixedRadix(domain) => domain.ifft_in_place(evals),
        }
    }

    #[inline]
    fn coset_fft_in_place<T: DomainCoeff<F>>(&self, coeffs: &mut Vec<T>) {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.coset_fft_in_place(coeffs),
            GeneralEvaluationDomain::MixedRadix(domain) => domain.coset_fft_in_place(coeffs),
        }
    }

    #[inline]
    fn coset_ifft_in_place<T: DomainCoeff<F>>(&self, evals: &mut Vec<T>) {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.coset_ifft_in_place(evals),
            GeneralEvaluationDomain::MixedRadix(domain) => domain.coset_ifft_in_place(evals),
        }
    }

    #[inline]
    fn evaluate_all_lagrange_coefficients(&self, tau: F) -> Vec<F> {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => {
                domain.evaluate_all_lagrange_coefficients(tau)
            }
            GeneralEvaluationDomain::MixedRadix(domain) => {
                domain.evaluate_all_lagrange_coefficients(tau)
            }
        }
    }

    #[inline]
    fn vanishing_polynomial(&self) -> crate::SparsePolynomial<F> {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.vanishing_polynomial(),
            GeneralEvaluationDomain::MixedRadix(domain) => domain.vanishing_polynomial(),
        }
    }

    #[inline]
    fn evaluate_vanishing_polynomial(&self, tau: F) -> F {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.evaluate_vanishing_polynomial(tau),
            GeneralEvaluationDomain::MixedRadix(domain) => {
                domain.evaluate_vanishing_polynomial(tau)
            }
        }
    }

    /// Return an iterator over the elements of the domain.
    fn elements(&self) -> GeneralElements<F> {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => {
                GeneralElements::BasicElements(domain.elements())
            }
            GeneralEvaluationDomain::MixedRadix(domain) => {
                GeneralElements::BasicElements(domain.elements())
            }
        }
    }
}

/// A generalized version of an iterator over the elements of a domain.
pub enum GeneralElements<F: FftField> {
    /// A basic iterator over the elements of a domain (currently, the only one in use).
    BasicElements(Elements<F>),
}

impl<F: FftField> Iterator for GeneralElements<F> {
    type Item = F;

    #[inline]
    fn next(&mut self) -> Option<F> {
        match self {
            GeneralElements::BasicElements(it) => it.next(),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{EvaluationDomain, GeneralEvaluationDomain};
    use algebra::{bls12_381::Fr, mnt6_753::Fr as MNT6Fr};
    use algebra_core::{test_rng, Zero};
    use rand::Rng;

    #[test]
    fn vanishing_polynomial_evaluation() {
        let rng = &mut test_rng();
        for coeffs in 0..10 {
            let domain = GeneralEvaluationDomain::<Fr>::new(coeffs).unwrap();
            let z = domain.vanishing_polynomial();
            for _ in 0..100 {
                let point = rng.gen();
                assert_eq!(
                    z.evaluate(point),
                    domain.evaluate_vanishing_polynomial(point)
                )
            }
        }

        for coeffs in 15..17 {
            let domain = GeneralEvaluationDomain::<MNT6Fr>::new(coeffs).unwrap();
            let z = domain.vanishing_polynomial();
            for _ in 0..100 {
                let point = rng.gen();
                assert_eq!(
                    z.evaluate(point),
                    domain.evaluate_vanishing_polynomial(point)
                )
            }
        }
    }

    #[test]
    fn vanishing_polynomial_vanishes_on_domain() {
        for coeffs in 0..1000 {
            let domain = GeneralEvaluationDomain::<Fr>::new(coeffs).unwrap();
            let z = domain.vanishing_polynomial();
            for point in domain.elements() {
                assert!(z.evaluate(point).is_zero())
            }
        }
    }

    #[test]
    fn size_of_elements() {
        for coeffs in 1..10 {
            let size = 1 << coeffs;
            let domain = GeneralEvaluationDomain::<Fr>::new(size).unwrap();
            let domain_size = domain.size();
            assert_eq!(domain_size, domain.elements().count());
        }
    }
}
