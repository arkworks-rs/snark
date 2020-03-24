use crate::{FpParameters, PrimeField};
use crate::{multicore::Worker, EvaluationDomainImpl};
use std::fmt;
use rayon::prelude::*;

/// Defines a domain over which finite field (I)FFTs can be performed. Works
/// only for fields that have a large multiplicative subgroup of size that is
/// a power-of-2. The roots of the target polynomial will be the powers of a
/// 2^n root of unity in the field. This allows us to perform polynomial operations
/// in O(n) by performing an O(n log n) FFT over such a domain.
#[derive(Copy, Clone, Hash, Eq, PartialEq, Default)]
pub struct BasicRadix2Domain<F: PrimeField> {
    /// The size of the domain.
    pub size:                  u64,
    /// `log_2(self.size)`.
    pub log_size_of_group:     u32,
    /// Size of the domain as a field element.
    pub size_as_field_element: F,
    /// Inverse of the size in the field.
    pub size_inv:              F,
    /// A generator of the subgroup.
    pub group_gen:             F,
    /// Inverse of the generator of the subgroup.
    pub group_gen_inv:         F,
    /// Multiplicative generator of the finite field.
    pub generator_inv:         F,
}

impl<F: PrimeField> fmt::Debug for BasicRadix2Domain<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Multiplicative subgroup of size {}", self.size)
    }
}

impl<F: PrimeField> BasicRadix2Domain<F> {

    /// Given an index which assumes the first elements of this domain are the elements of
    /// another (sub)domain with size size_s,
    /// this returns the actual index into this domain.
    pub fn reindex_by_subdomain(&self, other: Self, index: usize) -> usize {
        assert!(self.size() >= other.size());
        // Let this subgroup be G, and the subgroup we're re-indexing by be S.
        // Since its a subgroup, the 0th element of S is at index 0 in G, the first element of S is at
        // index |G|/|S|, the second at 2*|G|/|S|, etc.
        // Thus for an index i that corresponds to S, the index in G is i*|G|/|S|
        let period = self.size() / other.size();
        if index < other.size() {
            index * period
        } else {
            // Let i now be the index of this element in G \ S
            // Let x be the number of elements in G \ S, for every element in S. Then x = (|G|/|S| - 1).
            // At index i in G \ S, the number of elements in S that appear before the index in G to which
            // i corresponds to, is floor(i / x) + 1.
            // The +1 is because index 0 of G is S_0, so the position is offset by at least one.
            // The floor(i / x) term is because after x elements in G \ S, there is one more element from S
            // that will have appeared in G.
            let i = index - other.size();
            let x = period - 1;
            i + (i / x) + 1
        }
    }

    fn best_fft(a: &mut [F], worker: &Worker, omega: F, log_n: u32) {
        let log_cpus = worker.log_num_cpus();

        if log_n <= log_cpus {
            Self::serial_fft(a, omega, log_n);
        } else {
            Self::parallel_fft(a, worker, omega, log_n, log_cpus);
        }
    }

    pub(crate) fn serial_fft(a: &mut [F], omega: F, log_n: u32) {
        #[inline]
        fn bitreverse(mut n: u32, l: u32) -> u32 {
            let mut r = 0;
            for _ in 0..l {
                r = (r << 1) | (n & 1);
                n >>= 1;
            }
            r
        }
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
                    t *= &w;
                    let mut tmp = a[(k + j) as usize];
                    tmp -= &t;
                    a[(k + j + m) as usize] = tmp;
                    a[(k + j) as usize] += &t;
                    w.mul_assign(&w_m);
                }

                k += 2 * m;
            }

            m *= 2;
        }
    }

    pub(crate) fn parallel_fft(
        a: &mut [F],
        worker: &Worker,
        omega: F,
        log_n: u32,
        log_cpus: u32,
    ) {
        assert!(log_n >= log_cpus);

        let num_cpus = 1 << log_cpus;
        let log_new_n = log_n - log_cpus;
        let mut tmp = vec![vec![F::zero(); 1 << log_new_n]; num_cpus];
        let new_omega = omega.pow(&[num_cpus as u64]);

        worker.scope(0, |scope, _| {
            let a = &*a;

            for (j, tmp) in tmp.iter_mut().enumerate() {
                scope.spawn(move |_| {
                    // Shuffle into a sub-FFT
                    let omega_j = omega.pow(&[j as u64]);
                    let omega_step = omega.pow(&[(j as u64) << log_new_n]);

                    let mut elt = F::one();
                    for i in 0..(1 << log_new_n) {
                        for s in 0..num_cpus {
                            let idx = (i + (s << log_new_n)) % (1 << log_n);
                            let mut t = a[idx];
                            t *= &elt;
                            tmp[i] += &t;
                            elt *= &omega_step;
                        }
                        elt *= &omega_j;
                    }

                    // Perform sub-FFT
                    Self::serial_fft(tmp, new_omega, log_new_n);
                });
            }
        });

        worker.scope(a.len(), |scope, chunk| {
            let tmp = &tmp;

            for (idx, a) in a.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut idx = idx * chunk;
                    let mask = (1 << log_cpus) - 1;
                    for a in a {
                        *a = tmp[idx & mask][idx >> log_cpus];
                        idx += 1;
                    }
                });
            }
        });
    }
}

impl<F: PrimeField> EvaluationDomainImpl<F> for BasicRadix2Domain<F> {

    fn new(num_coeffs: usize) -> Option<Self>
    {
        // Compute the size of our evaluation domain
        let (size, log_size_of_group) = match Self::compute_size_of_domain(num_coeffs) {
            Some(size) => (size, size.trailing_zeros()),
            _ => return None,
        };

        // Compute the generator for the multiplicative subgroup.
        // It should be 2^(log_size_of_group) root of unity.
        let mut group_gen = F::root_of_unity();
        for _ in log_size_of_group..F::Params::TWO_ADICITY {
            group_gen.square_in_place();
        }

        let size_as_bigint = F::BigInt::from(size as u64);
        let size_as_field_element = F::from_repr(size_as_bigint);
        let size_inv = size_as_field_element.inverse()?;

        Some(Self{
            size: size as u64,
            log_size_of_group,
            size_as_field_element,
            size_inv,
            group_gen,
            group_gen_inv: group_gen.inverse()?,
            generator_inv: F::multiplicative_generator().inverse()?
        })
    }

    fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        let size = num_coeffs.next_power_of_two();
        if size.trailing_zeros() <= F::Params::TWO_ADICITY {
            Some(size)
        } else {
            None
        }
    }

    fn size(&self) -> usize {
        self.size.clone() as usize
    }

    fn size_inv(&self) -> F {
        self.size_inv.clone()
    }

    fn group_gen(&self) -> F {
        self.group_gen.clone()
    }

    fn fft_in_place(&self, coeffs: &mut Vec<F>) {
        coeffs.resize(self.size(), F::zero());
        Self::best_fft(coeffs, &Worker::new(), self.group_gen, self.log_size_of_group)
    }

    fn coset_fft_in_place(&self, coeffs: &mut Vec<F>) {
        Self::distribute_powers(coeffs, F::multiplicative_generator());
        self.fft_in_place(coeffs);
    }

    #[inline]
    fn ifft_in_place(&self, evals: &mut Vec<F>) {
        evals.resize(self.size(), F::zero());
        Self::best_fft(evals, &Worker::new(), self.group_gen_inv, self.log_size_of_group);
        evals.par_iter_mut().for_each(|val| *val *= &self.size_inv);
    }

    fn coset_ifft_in_place(&self, evals: &mut Vec<F>) {
        self.ifft_in_place(evals);
        Self::distribute_powers(evals, self.generator_inv);
    }
}