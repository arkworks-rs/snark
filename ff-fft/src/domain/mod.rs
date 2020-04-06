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
//! over the field is defined, we will try to build a mixed-radix evaluation
//! domain.

use self::utils::mixed_radix_fft_permute;
use crate::{domain::utils::best_mixed_domain_size, Vec};
use algebra_core::{fields::utils::k_adicity, FpParameters, PrimeField};
use core::fmt;
use rand::Rng;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub(crate) mod utils;

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
    ///
    /// If the field specifies a small subgroup for a mixed-radix FFT and
    /// the radix-2 FFT cannot be constructed, this method tries
    /// constructing a mixed-radix FFT instead.
    pub fn new(num_coeffs: usize) -> Option<Self> {
        let domain = Self::construct(num_coeffs);
        if domain.is_some() {
            return domain;
        }

        if F::Params::SMALL_SUBGROUP_BASE.is_some() {
            let s = best_mixed_domain_size::<F>(num_coeffs);
            return Self::construct(s);
        }

        None
    }

    /// Construct a domain that is large enough for evaluations of a polynomial
    /// having `num_coeffs` coefficients.
    pub fn construct(num_coeffs: usize) -> Option<Self> {
        // Compute the size of our evaluation domain
        let size;
        let log_size_of_group;
        let group_gen;

        if let Some(small_subgroup_base) = F::Params::SMALL_SUBGROUP_BASE {
            let q = small_subgroup_base as usize;
            let q_adicity = k_adicity(q, num_coeffs);
            let q_part = q.pow(q_adicity);

            let two_adicity = k_adicity(2, num_coeffs);
            let two_part = 1 << two_adicity;

            size = num_coeffs as u64;
            log_size_of_group = two_adicity;

            if num_coeffs != q_part * two_part {
                return None;
            }

            group_gen = F::get_root_of_unity(num_coeffs)?;
        } else {
            // Compute the size of our evaluation domain
            size = num_coeffs.next_power_of_two() as u64;
            log_size_of_group = size.trailing_zeros();

            // TODO: Check > vs. >= here.
            // libfqfft uses > https://github.com/scipr-lab/libfqfft/blob/e0183b2cef7d4c5deb21a6eaf3fe3b586d738fe0/libfqfft/evaluation_domain/domains/basic_radix2_domain.tcc#L33
            if log_size_of_group > F::Params::TWO_ADICITY {
                return None;
            }

            // Compute the generator for the multiplicative subgroup.
            // It should be 2^(log_size_of_group) root of unity.
            group_gen = F::get_root_of_unity(num_coeffs)?;
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
    fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        if let Some(small_subgroup_base) = F::Params::SMALL_SUBGROUP_BASE {
            let q = small_subgroup_base as usize;
            let q_adicity = k_adicity(q, num_coeffs);
            let q_part = q.pow(q_adicity);

            let two_adicity = k_adicity(2, num_coeffs);
            let two_part = 1 << two_adicity;

            if num_coeffs == q_part * two_part {
                Some(num_coeffs)
            } else {
                None
            }
        } else {
            let size = num_coeffs.next_power_of_two();
            // TODO: Check > vs. >= here.
            if size.trailing_zeros() <= F::Params::TWO_ADICITY {
                Some(size)
            } else {
                None
            }
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
#[inline]
pub(crate) fn best_fft<T: DomainCoeff<F>, F: PrimeField>(a: &mut [T], omega: F, log_n: u32) {
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
    if F::Params::SMALL_SUBGROUP_BASE.is_some() {
        serial_mixed_radix_fft(a, omega, log_n)
    } else {
        serial_radix2_fft(a, omega, log_n)
    }
}

pub(crate) fn serial_radix2_fft<T: DomainCoeff<F>, F: PrimeField>(
    a: &mut [T],
    omega: F,
    log_n: u32,
) {
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

pub(crate) fn serial_mixed_radix_fft<T: DomainCoeff<F>, F: PrimeField>(
    a: &mut [T],
    omega: F,
    two_adicity: u32,
) {
    // Conceptually, this FFT first splits into 2 sub-arrays two_adicity many times,
    // and then splits into q sub-arrays q_adicity many times.

    let n = a.len();
    let q = F::Params::SMALL_SUBGROUP_BASE.unwrap() as usize;

    let q_adicity = k_adicity(q, n);
    let q_part = q.pow(q_adicity);
    let two_part = 1 << two_adicity;

    assert_eq!(n, q_part * two_part);

    let mut m = 1; // invariant: m = 2^{s-1}

    if q_adicity > 0 {
        // If we're using the other radix, we have to do two things differently than in
        // the radix 2 case. 1. Applying the index permutation is a bit more
        // complicated. It isn't an involution (like it is in the radix 2 case)
        // so we need to remember which elements we've moved as we go along
        // and can't use the trick of just swapping when processing the first element of
        // a 2-cycle.
        //
        // 2. We need to do q_adicity many merge passes, each of which is a bit more
        // complicated than the specialized q=2 case.

        // Applying the permutation
        let mut seen = vec![false; n];
        for k in 0..n {
            let mut i = k;
            let mut a_i = a[i];
            while !seen[i] {
                let dest = mixed_radix_fft_permute(two_adicity, q_adicity, q, n, i);

                let a_dest = a[dest];
                a[dest] = a_i;

                seen[i] = true;

                a_i = a_dest;
                i = dest;
            }
        }

        let omega_q = omega.pow(&[(n / q) as u64]);
        let mut qth_roots = Vec::with_capacity(q);
        qth_roots.push(F::one());
        for i in 1..q {
            qth_roots.push(qth_roots[i - 1] * omega_q);
        }

        let mut terms = vec![T::zero(); q - 1];

        // Doing the q_adicity passes.
        for _ in 0..q_adicity {
            let w_m = omega.pow(&[(n / (q * m)) as u64]);
            let mut k = 0;
            while k < n {
                let mut w_j = F::one(); // w_j is omega_m ^ j
                for j in 0..m {
                    let base_term = a[k + j];
                    let mut w_j_i = w_j;
                    for i in 1..q {
                        terms[i - 1] = a[k + j + i * m];
                        terms[i - 1] *= w_j_i;
                        w_j_i *= w_j;
                    }

                    for i in 0..q {
                        a[k + j + i * m] = base_term;
                        for l in 1..q {
                            let mut tmp = terms[l - 1];
                            tmp *= qth_roots[(i * l) % q];
                            a[k + j + i * m] += tmp;
                        }
                    }

                    w_j *= w_m;
                }

                k += q * m;
            }
            m *= q;
        }
    } else {
        // swapping in place (from Storer's book)
        for k in 0..n {
            let rk = bitreverse(k as u32, two_adicity) as usize;
            if k < rk {
                a.swap(k, rk);
            }
        }
    }

    for _ in 0..two_adicity {
        // w_m is 2^s-th root of unity now
        let w_m = omega.pow(&[(n / (2 * m)) as u64]);

        let mut k = 0;
        while k < n {
            let mut w = F::one();
            for j in 0..m {
                let mut t = a[(k + m) + j];
                t *= w;
                a[(k + m) + j] = a[k + j];
                a[(k + m) + j] -= t;
                a[k + j] += t;
                w *= w_m;
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

    let m = a.len();
    let num_chunks = 1 << (log_cpus as usize);
    assert_eq!(m % num_chunks, 0);
    let m_div_num_chunks = m / num_chunks;

    let mut tmp = vec![vec![T::zero(); m_div_num_chunks]; num_chunks];
    let new_omega = omega.pow(&[num_chunks as u64]);
    let new_two_adicity = k_adicity(2, m_div_num_chunks);

    tmp.par_iter_mut().enumerate().for_each(|(j, tmp)| {
        // Shuffle into a sub-FFT
        let omega_j = omega.pow(&[j as u64]);
        let omega_step = omega.pow(&[(j * m_div_num_chunks) as u64]);

        let mut elt = F::one();
        for i in 0..m_div_num_chunks {
            for s in 0..num_chunks {
                let idx = (i + (s * m_div_num_chunks)) % m;
                let mut t = a[idx];
                t *= elt;
                tmp[i] += t;
                elt *= &omega_step;
            }
            elt *= &omega_j;
        }

        // Perform sub-FFT
        serial_fft(tmp, new_omega, new_two_adicity);
    });

    a.iter_mut()
        .enumerate()
        .for_each(|(i, a)| *a = tmp[i % num_chunks][i / num_chunks]);
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
