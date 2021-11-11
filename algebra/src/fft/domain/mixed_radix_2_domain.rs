use crate::{multicore::Worker, EvaluationDomain};
use crate::{FpParameters, PrimeField};
use rayon::prelude::*;
use std::any::Any;
use std::fmt;

/// Defines a domain over which finite field (I)FFTs can be performed. Works
/// only for fields that have a large multiplicative subgroup of size that is
/// a power-of-2 and a small multiplicative subgroup of size that is a power of
/// q, with q being a prime number. The evaluation domain will have size 2^k * q^s.
/// In this case, k steps of radix-2 algorithm are applied at the beginning of
/// the algorithm, while the rest of the transform is carried out by s steps of the
/// radix-q algorithm.
#[derive(Copy, Clone, Hash, Eq, PartialEq, Default)]
pub struct MixedRadix2Domain<F: PrimeField> {
    /// The size of the domain.
    pub size: u64,
    /// `log_2(self.size)`.
    pub log_size_of_group: u32,
    /// Inverse of the size in the field.
    pub size_inv: F,
    /// A generator of the subgroup.
    pub group_gen: F,
    /// Inverse of the generator of the subgroup.
    pub group_gen_inv: F,
    /// Multiplicative generator of the finite field.
    pub generator_inv: F,
}

impl<F: PrimeField> fmt::Debug for MixedRadix2Domain<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Multiplicative subgroup of size {}", self.size)
    }
}

impl<F: PrimeField> MixedRadix2Domain<F> {
    pub fn new(num_coeffs: usize) -> Option<Self> {
        let q = F::Params::SMALL_SUBGROUP_BASE.unwrap();

        // Compute the size of our evaluation domain
        let (size, log_size_of_group, q_adicity, two_adicity) =
            match Self::compute_size_of_domain(num_coeffs) {
                Some(size) => {
                    let q_adicity = Self::k_adicity(q, size as u64);
                    let two_adicity = Self::k_adicity(2, size as u64);
                    (size as u64, size.trailing_zeros(), q_adicity, two_adicity)
                }
                _ => return None,
            };

        let q_as_bigint = F::BigInt::from(q);
        let mut group_gen = match F::full_root_of_unity() {
            Some(v) => v,
            None => return None,
        };
        for _ in q_adicity..F::Params::SMALL_SUBGROUP_POWER.unwrap() {
            group_gen = group_gen.pow(&q_as_bigint);
        }
        for _ in two_adicity..F::Params::TWO_ADICITY as u64 {
            group_gen.square_in_place();
        }
        let size_as_bigint = F::BigInt::from(size);
        let size_as_field_element = F::from_repr(size_as_bigint);
        let size_inv = size_as_field_element.inverse()?;

        let group_gen_inv = group_gen.inverse()?;
        let generator_inv = F::multiplicative_generator().inverse()?;
        Some(Self {
            size,
            log_size_of_group,
            size_inv,
            group_gen,
            group_gen_inv,
            generator_inv,
        })
    }

    //Returns: min { n : N | n = 2^k * q^s, n >= num_coeffs, s <= small_subgroup_power, k <= TWO_ADICITY }
    pub fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        let mut best = std::u64::MAX;
        for b in 0..F::Params::SMALL_SUBGROUP_POWER.unwrap() + 1 {
            let mut r = F::Params::SMALL_SUBGROUP_BASE.unwrap().pow(b as u32);
            let mut two_adicity = 0;
            while r < num_coeffs as u64 {
                r *= 2;
                two_adicity += 1;
            }
            if two_adicity <= F::Params::TWO_ADICITY {
                best = best.min(r);
            }
        }

        let q = F::Params::SMALL_SUBGROUP_BASE.unwrap();
        let q_adicity = Self::k_adicity(q, best);
        let q_part = q.pow(q_adicity as u32);

        let two_adicity = Self::k_adicity(2, best);
        let two_part = 1 << two_adicity;

        if best != (two_part * q_part)
            || two_adicity > F::Params::TWO_ADICITY as u64
            || q_adicity > F::Params::SMALL_SUBGROUP_POWER.unwrap()
        {
            return None;
        }

        return Some(best as usize);
    }

    fn k_adicity(k: u64, n: u64) -> u64 {
        let mut r = 0;
        let mut ctr = n.clone();
        while ctr > 1 {
            if ctr % k == 0 {
                r += 1;
                ctr /= k;
            } else {
                return r;
            }
        }
        return r;
    }

    /// Algorithm from [CLRS 2n Ed, pp. 864]
    /// This FFT first splits into 2 sub-arrays two_adicity many times,
    /// And then splits into q sub-arrays q_adicity many times.
    pub(crate) fn mixed_serial_fft(a: &mut [F], omega: F, log_n: u32) {
        let n = a.len() as u64;
        let q = F::Params::SMALL_SUBGROUP_BASE.unwrap() as usize;
        let n_over_q = F::BigInt::from(n / q as u64);

        let q_adicity = Self::k_adicity(q as u64, n);
        let two_adicity = Self::k_adicity(2, n);
        let mut m = 1; // invariant: m = 2^{s-1}

        if q_adicity > 0 {
            // If we're using the other radix, we have to do two things differently than in the radix 2 case.
            // 1. Applying the index permutation is a bit more complicated. It isn't an involution
            // (like it is in the radix 2 case) so we need to remember which elements we've moved as we go along
            // and can't use the trick of just swapping when processing the first element of a 2-cycle.
            //
            // 2. We need to do q_adicity many merge passes, each of which is a bit more complicated than the
            // specialized q=2 case.

            // The algorithm reindexes the FFT domain C by reversing the digits of the mixed-radix representation
            // x = (b_0 + b_1*2 + ... + b_{two_adicity-1}* 2^{two_adicity-1}
            //        +  2^{two_adicity} (x_0  + x_1*q+..+ x_{q_adicity -1}*q^{q_adicity - 1}),
            // see the mixed_radix_fft_permute() below.

            // Applying the permutation
            let mut seen = vec![false; n as usize];
            for k in 0..n {
                let mut i = k as usize;
                let mut a_i = a[i];
                while !seen[i] {
                    let dest =
                        Self::mixed_radix_fft_permute(two_adicity, q_adicity, q as u64, n, i);
                    let a_dest = a[dest];
                    a[dest] = a_i;

                    seen[i] = true;

                    a_i = a_dest;
                    i = dest;
                }
            }

            // We recursively compute the FFTs over the cosets of C_{q*m} from the
            // FFTs over the cosets of C_m, starting with m=1.
            // With this convention,
            //          new_a[k' || i || j ] =  Sum_{l=0..q}  w^{ (i*m + j) * l} * a[k' || l || j]
            // where w is a generator of C_{q*m} (and hence w^{m*i*l} is a q-th
            // unit root
            //          qth_roots[i*l mod q] = g^{n/q* i*l}).
            let omega_q = omega.pow(n_over_q);
            let mut qth_roots = vec![F::one(); q];
            for i in 1..q {
                qth_roots[i] = qth_roots[i - 1] * &omega_q;
            }

            let mut terms = vec![F::one(); q - 1];

            for _ in 0..q_adicity {
                let n_over_q_times_m = F::BigInt::from(n / ((q * m) as u64));
                // w_m is the generator of the cyclic subgroup C_{q*m}
                let w_m = omega.pow(n_over_q_times_m);
                let mut k = 0;
                // k enumerates the partition of C_n into cosets of C_{q*m}
                while k < (n as usize) {
                    let mut w_j = F::one(); // w_j keeps track of omega_m ^ j
                                            // compute the FFT for the coset C_{q*m} at k.
                    for j in 0..m {
                        //  terms[i-1] = w^{i*j} * a[k || i || j], i= 1..q
                        let base_term = a[k + j];
                        let mut w_j_i = w_j.clone(); // w_j_i keeps track of the powers w^{j*i}
                        for i in 1..q {
                            terms[i - 1] = w_j_i * &a[k + j + (i * m)];
                            w_j_i *= &w_j;
                        }

                        for i in 0..q {
                            //  a[k || i || j] <-  Sum_{l=0..q}  w^{ (i*m + j) * l} * a[k' || l || j] =
                            //                  =  Sum_{l=0..q} qth_roots[(i*l)%q] * w^{ (l * j} * a[k' || l || j]
                            a[k + j + (i * m)] = base_term;
                            for l in 1..q {
                                a[k + j + (i * m)] += &(qth_roots[(i * l) % q] * &terms[l - 1]);
                            }
                        }
                        w_j *= &w_m;
                    }
                    // choose next coset of C_{q*m}
                    k += q * m;
                }
                m *= q;
            }
        } else {
            #[inline]
            fn bitreverse(mut n: u32, l: u32) -> u32 {
                let mut r = 0;
                for _ in 0..l {
                    r = (r << 1) | (n & 1);
                    n >>= 1;
                }
                r
            }

            // Swapping in place (from Storer's book)
            for k in 0..(n as u32) {
                let rk = bitreverse(k, log_n);
                if k < rk {
                    a.swap(rk as usize, k as usize);
                }
            }
        }

        //2-adic part
        for _ in 0..two_adicity {
            let w_m = omega.pow(&[(n / (2 * m as u64))]); // w_m is 2^s-th root of unity now
            let mut k = 0;
            while k < n as usize {
                let mut w = F::one();
                for j in 0..m {
                    let mut t = a[k + j + m];
                    t *= &w;
                    let mut tmp = a[k + j];
                    tmp -= &t;
                    a[k + j + m] = tmp;
                    a[k + j] += &t;
                    w *= &w_m;
                }

                k += 2 * m;
            }

            m *= 2;
        }
    }

    /// This is the permutation obtained by splitting into
    /// 2 groups two_adicity times and then
    /// q groups q_adicity many times
    /// It can be efficiently described as follows:
    /// write
    /// i = 2^0 b_0 + 2^1 b_1 + ... + 2^{two_adicity - 1} b_{two_adicity - 1}
    ///     + 2^two_adicity ( x_0 + q^1 x_1 + .. + q^{q_adicity-1} x_{q_adicity-1})
    /// We want to return
    /// j = b_0 (N/2) + b_1 (N/ 2^2) + ... + b_{two_adicity-1} (N/ 2^two_adicity)
    ///     + x_0 (N / 2^two_adicity / q) + .. + x_{q_adicity-1} (N / 2^two_adicity / q^q_adicity)
    fn mixed_radix_fft_permute(
        two_adicity: u64,
        q_adicity: u64,
        q: u64,
        n: u64,
        idx: usize,
    ) -> usize {
        let mut res = 0;
        let mut shift = n;
        let mut i = idx as u64;

        for _ in 0..two_adicity {
            shift = shift / 2;
            res += (i % 2) * shift;
            i = i / 2;
        }

        for _ in 0..q_adicity {
            shift = shift / q;
            res += (i % q) * shift;
            i = i / q;
        }

        return res as usize;
    }

    pub(crate) fn mixed_parallel_fft(
        a: &mut [F],
        worker: &Worker,
        omega: F,
        log_n: u32,
        log_cpus: u32,
    ) {
        let num_cpus = 1 << log_cpus;
        let m = a.len();
        let two_adicity = Self::k_adicity(2, m as u64);
        let two_part = 1 << two_adicity;

        if two_part < num_cpus as u64 {
            Self::mixed_serial_fft(a, omega, log_n);
            return;
        }

        let log_new_n = m / num_cpus;

        let mut tmp = vec![vec![F::zero(); log_new_n]; num_cpus];
        let new_omega = omega.pow(&[num_cpus as u64]);

        worker.scope(0, |scope, _| {
            let a = &*a;

            for (j, tmp) in tmp.iter_mut().enumerate() {
                scope.spawn(move |_| {
                    // Shuffle into a sub-FFT
                    let omega_j = omega.pow(&[j as u64]);
                    let omega_step = omega.pow(&[(j * log_new_n) as u64]);

                    let mut elt = F::one();
                    for i in 0..log_new_n {
                        for s in 0..num_cpus {
                            let idx = (i + (s * log_new_n)) % m;
                            let mut t = a[idx];
                            t *= &elt;
                            tmp[i] += &t;
                            elt *= &omega_step;
                        }
                        elt *= &omega_j;
                    }
                    // Perform sub-FFT
                    Self::mixed_serial_fft(tmp, new_omega, log_new_n as u32);
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

    fn best_fft(a: &mut [F], _worker: &Worker, omega: F, log_n: u32) {
        let log_cpus = _worker.log_num_cpus();

        if log_n <= log_cpus {
            return Self::mixed_serial_fft(a, omega, log_n);
        } else {
            return Self::mixed_parallel_fft(a, _worker, omega, log_n, log_cpus);
        }
    }
}

impl<F: PrimeField> EvaluationDomain<F> for MixedRadix2Domain<F> {
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
        Self::best_fft(
            coeffs,
            &Worker::new(),
            self.group_gen,
            self.log_size_of_group,
        )
    }

    fn coset_fft_in_place(&self, coeffs: &mut Vec<F>) {
        Self::distribute_powers(coeffs, F::multiplicative_generator());
        self.fft_in_place(coeffs);
    }

    #[inline]
    fn ifft_in_place(&self, evals: &mut Vec<F>) {
        evals.resize(self.size(), F::zero());
        Self::best_fft(
            evals,
            &Worker::new(),
            self.group_gen_inv,
            self.log_size_of_group,
        );
        evals.par_iter_mut().for_each(|val| *val *= &self.size_inv);
    }

    fn coset_ifft_in_place(&self, evals: &mut Vec<F>) {
        self.ifft_in_place(evals);
        Self::distribute_powers(evals, self.generator_inv);
    }

    fn eq(&self, other: &dyn EvaluationDomain<F>) -> bool {
        other
            .as_any()
            .downcast_ref::<Self>()
            .map_or(false, |x| x == self)
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn clone_and_box(&self) -> Box<dyn EvaluationDomain<F>> {
        Box::new((*self).clone())
    }
}
