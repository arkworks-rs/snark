use crate::{FpParameters, PrimeField};
use crate::{multicore::Worker, EvaluationDomain};
use std::fmt;
use rayon::prelude::*;
use std::any::Any;

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

    pub fn new(num_coeffs: usize) -> Option<Self>
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

    pub fn compute_size_of_domain(num_coeffs: usize) -> Option<usize> {
        let size = num_coeffs.next_power_of_two();
        if size.trailing_zeros() <= F::Params::TWO_ADICITY {
            Some(size)
        } else {
            None
        }
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

    fn best_fft(a: &mut [F], worker: &Worker, omega: F, log_n: u32) {
        let log_cpus = worker.log_num_cpus();

        if log_n <= log_cpus {
            Self::serial_fft(a, omega, log_n);
        } else {
            Self::parallel_fft(a, worker, omega, log_n, log_cpus);
        }
    }

    /// Computes the radix-2 FFT of a[0..n] over an FFT domain {z: z^n - 1 = 0} of size n=2^log_n, 
    /// given a generator  omega for this domain.
    ///
    /// The algorithm reindexes the FFT domain C_n={0,1}^log(n) by reversing the bit order. 
    /// This makes the enumeration of the FFT domain C_n into the cosets of C_m with m|n slightly
    /// more intuitive: For fixed k' from 0..n/m, 
    ///              [k' || j] = k' * m + j, 
    /// where j=0..m, goes through the coset of C_m at k', and varying k' enumerates the partition
    /// of C_n.
    pub(crate) fn serial_fft(a: &mut [F], omega: F, log_n: u32) {
        // inverts the bit order of an l bit integer n
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
        debug_assert_eq!(n, 1 << log_n);

        // reindex a[0..n] as described above
        for k in 0..n {
            let rk = bitreverse(k, log_n);
            if k < rk {
                a.swap(rk as usize, k as usize);
            }
        }

        // We recursively compute the FFTs over the cosets of C_{2*m} from the 
        // FFTs over the cosets of C_m, starting with m=1. 
        // With this convention, 
        //          new_a[k' || 0 || j] = a[k' || 0 || j] +   w^{j} * a[k' || 1 || j]
        //          new_a[k' || 1 || j] = a[k' || 0 || j] + w^{m+j} * a[k' || 1 || j],
        // where w is a generator of C_{2m} (and hence w^{m+j} = -w^j).
        let mut m = 1;
        for _ in 0..log_n {
            // w_m is the generator of the cyclic subgroup C_{2m}
            let w_m = omega.pow(&[(n / (2 * m)) as u64]);

            let mut k = 0;
            // k enumerates the partition of C_n into cosets of C_{2m}
            while k < n {
                let mut w = F::one();
                // compute the FFT for the coset C_{2m} at k. 
                for j in 0..m {
                    //     a[k + m + j] <- a[k+j] - w_m * a[k + j + m] 
                    //     a[k + j]     <- a[k+j] + w_m * a[k + j + m]
                    let mut t = a[(k + j + m) as usize];
                    t *= &w;
                    let mut tmp = a[(k + j) as usize];
                    tmp -= &t;
                    a[(k + j + m) as usize] = tmp;
                    a[(k + j) as usize] += &t;
                    w.mul_assign(&w_m);
                }
                // choose next coset of C_{2m}
                k += 2 * m;
            }

            m *= 2;
        }
    }
    
    /// To parallelize over cpu=2^log_cpu cores, we split the computation of the FFT over 
    /// C_n = C_cpu x C_new in the following manner:
    /// FFT(f(x), k) = Sum_{(g,h) in C_cpus x C_new} f(gh x)* omega^{k * gh x} =
    ///             = Sum_{h in C_new} [ Sum_{g in C_cpus} f(g hx) * omega^{j *g hx} ] * omega^{i * cpus * hx}
    /// where k = i*cpus + j, with j in 0..cpus. 
    ///
    /// The inner sums 
    ///          f_j(x)=  Sum_{g in C_cpus} f(g x) * omega^{j *g x},
    /// are computed in a preparation step, and the "big" ones
    ///          phi_j(i) = Sum_{h in C_new} f_j(h x) omega^{i*cpus*hx},
    /// i=0..n/cpus, via a call of serial_fft.
    pub(crate) fn parallel_fft(
        a: &mut [F],
        worker: &Worker,
        omega: F,
        log_n: u32,
        log_cpus: u32,
    ) {
        debug_assert!(log_n >= log_cpus);

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
                    
                    // Compute f_j(x) for x in C_new.
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
                    // phi_j = FFT(f_j) over C_new
                    Self::serial_fft(tmp, new_omega, log_new_n);
                });
            }
        });

        // j=0..cpus, i=0..n/cpus,
        // FFT(f, idx = i*cpus + j) = phi_j(i) = FFT(f_j)[i] 
        worker.scope(a.len(), |scope, chunk| {
            let tmp = &tmp;

            for (idx, a) in a.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut idx = idx * chunk; // compute index from chunk index
                    let mask = (1 << log_cpus) - 1;
                    for a in a {
                        *a = tmp[idx & mask][idx >> log_cpus]; // idx & mask = idx mod cpus = j, idx>>log_cpus = i. 
                        idx += 1;
                    }
                });
            }
        });
    }
}

impl<F: PrimeField> EvaluationDomain<F> for BasicRadix2Domain<F> {

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

    fn eq(&self, other: & dyn EvaluationDomain<F>) -> bool {
        other.as_any().downcast_ref::<Self>().map_or(false, |x| x == self)
    }

    fn as_any(&self) -> & dyn Any {
        self
    }

    fn clone_and_box(&self) -> Box<dyn EvaluationDomain<F>> {
        Box::new((*self).clone())
    }
}