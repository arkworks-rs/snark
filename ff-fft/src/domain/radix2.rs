//! This module defines `Radix2EvaluationDomain`, an `EvaluationDomain`
//! for performing various kinds of polynomial arithmetic on top of
//! fields that are FFT-friendly. `Radix2EvaluationDomain` supports
//! FFTs of size at most `2^F::TWO_ADICITY`.

pub use crate::domain::utils::Elements;
use crate::domain::{
    utils::{best_fft, bitreverse},
    DomainCoeff, EvaluationDomain,
};
use crate::Vec;
use algebra_core::{FftField, FftParameters};
use core::convert::TryFrom;
use core::fmt;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Defines a domain over which finite field (I)FFTs can be performed. Works
/// only for fields that have a large multiplicative subgroup of size that is
/// a power-of-2.
#[derive(Copy, Clone, Hash, Eq, PartialEq)]
pub struct Radix2EvaluationDomain<F: FftField> {
    /// The size of the domain.
    pub size: u64,
    /// `log_2(self.size)`.
    pub log_size_of_group: u32,
    /// Size of the domain as a field element.
    pub size_as_field_element: F,
    /// Inverse of the size in the field.
    pub size_inv: F,
    /// A generator of the subgroup.
    pub group_gen: F,
    /// Inverse of the generator of the subgroup.
    pub group_gen_inv: F,
    /// Multiplicative generator of the finite field.
    pub generator_inv: F,
}

impl<F: FftField> fmt::Debug for Radix2EvaluationDomain<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Radix-2 multiplicative subgroup of size {}", self.size)
    }
}

impl<F: FftField> EvaluationDomain<F> for Radix2EvaluationDomain<F> {
    type Elements = Elements<F>;

    /// Construct a domain that is large enough for evaluations of a polynomial
    /// having `num_coeffs` coefficients.
    fn new(num_coeffs: usize) -> Option<Self> {
        // Compute the size of our evaluation domain
        let size = num_coeffs.next_power_of_two() as u64;
        let log_size_of_group = size.trailing_zeros();

        // libfqfft uses > https://github.com/scipr-lab/libfqfft/blob/e0183b2cef7d4c5deb21a6eaf3fe3b586d738fe0/libfqfft/evaluation_domain/domains/basic_radix2_domain.tcc#L33
        if log_size_of_group > F::FftParams::TWO_ADICITY {
            return None;
        }

        // Compute the generator for the multiplicative subgroup.
        // It should be the 2^(log_size_of_group) root of unity.
        let group_gen = F::get_root_of_unity(usize::try_from(size).unwrap())?;
        // Check that it is indeed the 2^(log_size_of_group) root of unity.
        debug_assert_eq!(group_gen.pow([size]), F::one());
        let size_as_field_element = F::from(size);
        let size_inv = size_as_field_element.inverse()?;

        Some(Radix2EvaluationDomain {
            size,
            log_size_of_group,
            size_as_field_element,
            size_inv,
            group_gen,
            group_gen_inv: group_gen.inverse()?,
            generator_inv: F::multiplicative_generator().inverse()?,
        })
    }

    fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        let size = num_coeffs.next_power_of_two();
        if size.trailing_zeros() > F::FftParams::TWO_ADICITY {
            None
        } else {
            Some(size)
        }
    }

    #[inline]
    fn size(&self) -> usize {
        usize::try_from(self.size).unwrap()
    }

    #[inline]
    fn fft_in_place<T: DomainCoeff<F>>(&self, coeffs: &mut Vec<T>) {
        coeffs.resize(self.size(), T::zero());
        best_fft(
            coeffs,
            self.group_gen,
            self.log_size_of_group,
            serial_radix2_fft::<T, F>,
        )
    }

    #[inline]
    fn ifft_in_place<T: DomainCoeff<F>>(&self, evals: &mut Vec<T>) {
        evals.resize(self.size(), T::zero());
        best_fft(
            evals,
            self.group_gen_inv,
            self.log_size_of_group,
            serial_radix2_fft::<T, F>,
        );
        cfg_iter_mut!(evals).for_each(|val| *val *= self.size_inv);
    }

    #[inline]
    fn coset_ifft_in_place<T: DomainCoeff<F>>(&self, evals: &mut Vec<T>) {
        self.ifft_in_place(evals);
        Self::distribute_powers(evals, self.generator_inv);
    }

    fn evaluate_all_lagrange_coefficients(&self, tau: F) -> Vec<F> {
        // Evaluate all Lagrange polynomials
        let size = self.size();
        let t_size = tau.pow(&[self.size]);
        let one = F::one();
        if t_size.is_one() {
            let mut u = vec![F::zero(); size];
            let mut omega_i = one;
            for i in 0..size {
                if omega_i == tau {
                    u[i] = one;
                    break;
                }
                omega_i *= &self.group_gen;
            }
            u
        } else {
            use algebra_core::fields::batch_inversion;

            let mut l = (t_size - &one) * &self.size_inv;
            let mut r = one;
            let mut u = vec![F::zero(); size];
            let mut ls = vec![F::zero(); size];
            for i in 0..size {
                u[i] = tau - &r;
                ls[i] = l;
                l *= &self.group_gen;
                r *= &self.group_gen;
            }

            batch_inversion(u.as_mut_slice());

            cfg_iter_mut!(u).zip(ls).for_each(|(tau_minus_r, l)| {
                *tau_minus_r = l * *tau_minus_r;
            });

            u
        }
    }

    fn vanishing_polynomial(&self) -> crate::SparsePolynomial<F> {
        let coeffs = vec![(0, -F::one()), (self.size(), F::one())];
        crate::SparsePolynomial::from_coefficients_vec(coeffs)
    }

    /// This evaluates the vanishing polynomial for this domain at tau.
    /// For multiplicative subgroups, this polynomial is `z(X) = X^self.size -
    /// 1`.
    fn evaluate_vanishing_polynomial(&self, tau: F) -> F {
        tau.pow(&[self.size]) - &F::one()
    }

    /// Return an iterator over the elements of the domain.
    fn elements(&self) -> Elements<F> {
        Elements {
            cur_elem: F::one(),
            cur_pow: 0,
            size: self.size,
            group_gen: self.group_gen,
        }
    }
}

pub(crate) fn serial_radix2_fft<T: DomainCoeff<F>, F: FftField>(a: &mut [T], omega: F, log_n: u32) {
    let n =
        u32::try_from(a.len()).expect("cannot perform FFTs larger on vectors of len > (1 << 32)");
    assert_eq!(n, 1 << log_n);

    for k in 0..n {
        let rk = bitreverse(k, log_n);
        if k < rk {
            a.swap(rk as usize, k as usize);
        }
    }

    let mut m = 1;
    for _ in 0..log_n {
        let w_m = omega.pow(&[(n / (2 * m)) as u64]);

        let mut k = 0;
        while k < n {
            let mut w = F::one();
            for j in 0..m {
                let mut t = a[(k + j + m) as usize];
                t *= w;
                let mut tmp = a[(k + j) as usize];
                tmp -= t;
                a[(k + j + m) as usize] = tmp;
                a[(k + j) as usize] += t;
                w.mul_assign(&w_m);
            }

            k += 2 * m;
        }

        m *= 2;
    }
}

#[cfg(test)]
mod tests {
    use crate::{EvaluationDomain, Radix2EvaluationDomain};
    use algebra::bls12_381::Fr;
    use algebra_core::{test_rng, Field, Zero};
    use rand::Rng;

    #[test]
    fn vanishing_polynomial_evaluation() {
        let rng = &mut test_rng();
        for coeffs in 0..10 {
            let domain = Radix2EvaluationDomain::<Fr>::new(coeffs).unwrap();
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
            let domain = Radix2EvaluationDomain::<Fr>::new(coeffs).unwrap();
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
            let domain = Radix2EvaluationDomain::<Fr>::new(size).unwrap();
            let domain_size = domain.size();
            assert_eq!(domain_size, domain.elements().count());
        }
    }

    #[test]
    fn elements_contents() {
        for coeffs in 1..10 {
            let size = 1 << coeffs;
            let domain = Radix2EvaluationDomain::<Fr>::new(size).unwrap();
            for (i, element) in domain.elements().enumerate() {
                assert_eq!(element, domain.group_gen.pow([i as u64]));
            }
        }
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn parallel_fft_consistency() {
        use super::serial_radix2_fft;
        use crate::domain::utils::parallel_fft;
        use crate::Vec;
        use algebra::bls12_381::Bls12_381;
        use algebra_core::{test_rng, PairingEngine, UniformRand};
        use core::cmp::min;

        fn test_consistency<E: PairingEngine, R: Rng>(rng: &mut R, max_coeffs: u32) {
            for _ in 0..5 {
                for log_d in 0..max_coeffs {
                    let d = 1 << log_d;

                    let mut v1 = (0..d).map(|_| E::Fr::rand(rng)).collect::<Vec<_>>();
                    let mut v2 = v1.clone();

                    let domain = Radix2EvaluationDomain::new(v1.len()).unwrap();

                    for log_cpus in log_d..min(log_d + 1, 3) {
                        parallel_fft::<E::Fr, E::Fr>(
                            &mut v1,
                            domain.group_gen,
                            log_d,
                            log_cpus,
                            serial_radix2_fft::<E::Fr, E::Fr>,
                        );
                        serial_radix2_fft::<E::Fr, E::Fr>(&mut v2, domain.group_gen, log_d);

                        assert_eq!(v1, v2);
                    }
                }
            }
        }

        let rng = &mut test_rng();

        test_consistency::<Bls12_381, _>(rng, 10);
    }
}
