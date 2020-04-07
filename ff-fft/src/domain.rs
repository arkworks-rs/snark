//! This module contains an `EvaluationDomain` abstraction for
//! performing various kinds of polynomial arithmetic on top of
//! the scalar field.
//!
//! In pairing-based SNARKs like GM17, we need to calculate
//! a quotient polynomial over a target polynomial with roots
//! at distinct points associated with each constraint of the
//! constraint system. In order to be efficient, we choose these
//! roots to be the powers of a 2^n root of unity in the field.
//! This allows us to perform polynomial operations in O(n)
//! by performing an O(n log n) FFT over such a domain.

use crate::Vec;
use algebra_core::{FpParameters, PrimeField};
use core::fmt;
use rand::Rng;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Defines a domain over which finite field (I)FFTs can be performed. Works
/// only for fields that have a large multiplicative subgroup of size that is
/// a power-of-2.
#[derive(Copy, Clone, Hash, Eq, PartialEq)]
pub struct EvaluationDomain<F: PrimeField> {
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

impl<F: PrimeField> fmt::Debug for EvaluationDomain<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Multiplicative subgroup of size {}", self.size)
    }
}

/// Types that can be FFT-ed must implement this trait.
pub trait DomainCoeff<F: PrimeField>:
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
    F: PrimeField,
    T: Copy
        + Send
        + Sync
        + core::ops::AddAssign
        + core::ops::SubAssign
        + algebra_core::Zero
        + core::ops::MulAssign<F>,
{
}

impl<F: PrimeField> EvaluationDomain<F> {
    /// Sample an element that is *not* in the domain.
    pub fn sample_element_outside_domain<R: Rng>(&self, rng: &mut R) -> F {
        let mut t = F::rand(rng);
        while self.evaluate_vanishing_polynomial(t).is_zero() {
            t = F::rand(rng);
        }
        t
    }

    /// Construct a domain that is large enough for evaluations of a polynomial
    /// having `num_coeffs` coefficients.
    pub fn new(num_coeffs: usize) -> Option<Self> {
        // Compute the size of our evaluation domain
        let size = num_coeffs.next_power_of_two() as u64;
        let log_size_of_group = size.trailing_zeros();

        if log_size_of_group >= F::Params::TWO_ADICITY {
            return None;
        }

        // Compute the generator for the multiplicative subgroup.
        // It should be 2^(log_size_of_group) root of unity.
        let mut group_gen = F::root_of_unity();
        for _ in log_size_of_group..F::Params::TWO_ADICITY {
            group_gen.square_in_place();
        }

        let size_as_bigint = F::BigInt::from(size);
        let size_as_field_element = F::from_repr(size_as_bigint);
        let size_inv = size_as_field_element.inverse()?;

        Some(EvaluationDomain {
            size,
            log_size_of_group,
            size_as_field_element,
            size_inv,
            group_gen,
            group_gen_inv: group_gen.inverse()?,
            generator_inv: F::multiplicative_generator().inverse()?,
        })
    }

    /// Return the size of a domain that is large enough for evaluations of a
    /// polynomial having `num_coeffs` coefficients.
    pub fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        let size = num_coeffs.next_power_of_two();
        if size.trailing_zeros() < F::Params::TWO_ADICITY {
            Some(size)
        } else {
            None
        }
    }

    /// Return the size of `self`.
    pub fn size(&self) -> usize {
        self.size as usize
    }

    /// Compute a FFT.
    pub fn fft<T: DomainCoeff<F>>(&self, coeffs: &[T]) -> Vec<T> {
        let mut coeffs = coeffs.to_vec();
        self.fft_in_place(&mut coeffs);
        coeffs
    }

    /// Compute a FFT, modifying the vector in place.
    pub fn fft_in_place<T: DomainCoeff<F>>(&self, coeffs: &mut Vec<T>) {
        coeffs.resize(self.size(), T::zero());
        best_fft(coeffs, self.group_gen, self.log_size_of_group)
    }

    /// Compute a IFFT.
    pub fn ifft<T: DomainCoeff<F>>(&self, evals: &[T]) -> Vec<T> {
        let mut evals = evals.to_vec();
        self.ifft_in_place(&mut evals);
        evals
    }

    /// Compute a IFFT, modifying the vector in place.
    #[inline]
    pub fn ifft_in_place<T: DomainCoeff<F>>(&self, evals: &mut Vec<T>) {
        evals.resize(self.size(), T::zero());
        best_fft(evals, self.group_gen_inv, self.log_size_of_group);
        cfg_iter_mut!(evals).for_each(|val| *val *= self.size_inv);
    }

    fn distribute_powers<T: DomainCoeff<F>>(coeffs: &mut [T], g: F) {
        let mut pow = F::one();
        coeffs.iter_mut().for_each(|c| {
            *c *= pow;
            pow *= &g
        })
    }

    /// Compute a FFT over a coset of the domain.
    pub fn coset_fft<T: DomainCoeff<F>>(&self, coeffs: &[T]) -> Vec<T> {
        let mut coeffs = coeffs.to_vec();
        self.coset_fft_in_place(&mut coeffs);
        coeffs
    }

    /// Compute a FFT over a coset of the domain, modifying the input vector
    /// in place.
    pub fn coset_fft_in_place<T: DomainCoeff<F>>(&self, coeffs: &mut Vec<T>) {
        Self::distribute_powers(coeffs, F::multiplicative_generator());
        self.fft_in_place(coeffs);
    }

    /// Compute a IFFT over a coset of the domain.
    pub fn coset_ifft<T: DomainCoeff<F>>(&self, evals: &[T]) -> Vec<T> {
        let mut evals = evals.to_vec();
        self.coset_ifft_in_place(&mut evals);
        evals
    }

    /// Compute a IFFT over a coset of the domain, modifying the input vector in
    /// place.
    pub fn coset_ifft_in_place<T: DomainCoeff<F>>(&self, evals: &mut Vec<T>) {
        self.ifft_in_place(evals);
        Self::distribute_powers(evals, self.generator_inv);
    }

    /// Evaluate all the lagrange polynomials defined by this domain at the
    /// point `tau`.
    pub fn evaluate_all_lagrange_coefficients(&self, tau: F) -> Vec<F> {
        // Evaluate all Lagrange polynomials
        let size = self.size as usize;
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

    /// Return the sparse vanishing polynomial.
    pub fn vanishing_polynomial(&self) -> crate::SparsePolynomial<F> {
        let coeffs = vec![(0, -F::one()), (self.size(), F::one())];
        crate::SparsePolynomial::from_coefficients_vec(coeffs)
    }

    /// This evaluates the vanishing polynomial for this domain at tau.
    /// For multiplicative subgroups, this polynomial is `z(X) = X^self.size -
    /// 1`.
    pub fn evaluate_vanishing_polynomial(&self, tau: F) -> F {
        tau.pow(&[self.size]) - &F::one()
    }

    /// Return an iterator over the elements of the domain.
    pub fn elements(&self) -> Elements<F> {
        Elements {
            cur_elem: F::one(),
            cur_pow: 0,
            domain: *self,
        }
    }

    /// The target polynomial is the zero polynomial in our
    /// evaluation domain, so we must perform division over
    /// a coset.
    pub fn divide_by_vanishing_poly_on_coset_in_place(&self, evals: &mut [F]) {
        let i = self
            .evaluate_vanishing_polynomial(F::multiplicative_generator())
            .inverse()
            .unwrap();

        cfg_iter_mut!(evals).for_each(|eval| *eval *= &i);
    }

    /// Given an index which assumes the first elements of this domain are the
    /// elements of another (sub)domain with size size_s,
    /// this returns the actual index into this domain.
    pub fn reindex_by_subdomain(&self, other: Self, index: usize) -> usize {
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
    pub fn mul_polynomials_in_evaluation_domain(
        &self,
        self_evals: &[F],
        other_evals: &[F],
    ) -> Vec<F> {
        assert_eq!(self_evals.len(), other_evals.len());
        let mut result = self_evals.to_vec();

        cfg_iter_mut!(result)
            .zip(other_evals)
            .for_each(|(a, b)| *a *= b);

        result
    }
}

#[cfg(feature = "parallel")]
fn best_fft<T: DomainCoeff<F>, F: PrimeField>(a: &mut [T], omega: F, log_n: u32) {
    fn log2_floor(num: usize) -> u32 {
        assert!(num > 0);
        let mut pow = 0;
        while (1 << (pow + 1)) <= num {
            pow += 1;
        }
        pow
    }

    let num_cpus = rayon::current_num_threads();
    let log_cpus = log2_floor(num_cpus);
    if log_n <= log_cpus {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, omega, log_n, log_cpus);
    }
}

#[cfg(not(feature = "parallel"))]
fn best_fft<T: DomainCoeff<F>, F: PrimeField>(a: &mut [T], omega: F, log_n: u32) {
    serial_fft(a, omega, log_n)
}

#[inline]
fn bitreverse(mut n: u32, l: u32) -> u32 {
    let mut r = 0;
    for _ in 0..l {
        r = (r << 1) | (n & 1);
        n >>= 1;
    }
    r
}

pub(crate) fn serial_fft<T: DomainCoeff<F>, F: PrimeField>(a: &mut [T], omega: F, log_n: u32) {
    let n = a.len() as u32;
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

#[cfg(feature = "parallel")]
pub(crate) fn parallel_fft<T: DomainCoeff<F>, F: PrimeField>(
    a: &mut [T],
    omega: F,
    log_n: u32,
    log_cpus: u32,
) {
    assert!(log_n >= log_cpus);

    let num_cpus = 1 << log_cpus;
    let log_new_n = log_n - log_cpus;
    let mut tmp = vec![vec![T::zero(); 1 << log_new_n]; num_cpus];
    let new_omega = omega.pow(&[num_cpus as u64]);

    tmp.par_iter_mut().enumerate().for_each(|(j, tmp)| {
        // Shuffle into a sub-FFT
        let omega_j = omega.pow(&[j as u64]);
        let omega_step = omega.pow(&[(j as u64) << log_new_n]);

        let mut elt = F::one();
        for i in 0..(1 << log_new_n) {
            for s in 0..num_cpus {
                let idx = (i + (s << log_new_n)) % (1 << log_n);
                let mut t = a[idx];
                t *= elt;
                tmp[i] += t;
                elt *= &omega_step;
            }
            elt *= &omega_j;
        }

        // Perform sub-FFT
        serial_fft(tmp, new_omega, log_new_n);
    });

    let mask = (1 << log_cpus) - 1;
    a.iter_mut()
        .enumerate()
        .for_each(|(i, a)| *a = tmp[i & mask][i >> log_cpus]);
}

/// An iterator over the elements of the domain.
pub struct Elements<F: PrimeField> {
    cur_elem: F,
    cur_pow: u64,
    domain: EvaluationDomain<F>,
}

impl<F: PrimeField> Iterator for Elements<F> {
    type Item = F;
    fn next(&mut self) -> Option<F> {
        if self.cur_pow == self.domain.size {
            None
        } else {
            let cur_elem = self.cur_elem;
            self.cur_elem *= &self.domain.group_gen;
            self.cur_pow += 1;
            Some(cur_elem)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::EvaluationDomain;
    use algebra::bls12_381::Fr;
    use algebra_core::{test_rng, Field, Zero};
    use rand::Rng;

    #[test]
    fn vanishing_polynomial_evaluation() {
        let rng = &mut test_rng();
        for coeffs in 0..10 {
            let domain = EvaluationDomain::<Fr>::new(coeffs).unwrap();
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
            let domain = EvaluationDomain::<Fr>::new(coeffs).unwrap();
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
            let domain = EvaluationDomain::<Fr>::new(size).unwrap();
            let domain_size = domain.size();
            assert_eq!(domain_size, domain.elements().count());
        }
    }

    #[test]
    fn elements_contents() {
        for coeffs in 1..10 {
            let size = 1 << coeffs;
            let domain = EvaluationDomain::<Fr>::new(size).unwrap();
            for (i, element) in domain.elements().enumerate() {
                assert_eq!(element, domain.group_gen.pow([i as u64]));
            }
        }
    }
}
