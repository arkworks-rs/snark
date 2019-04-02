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

use crate::{Field, FpParameters, PairingEngine, PrimeField, ProjectiveCurve, SquareRootField};

use super::multicore::Worker;
use std::{
    marker::PhantomData,
    ops::{AddAssign, MulAssign, SubAssign},
};

pub struct EvaluationDomain<E: PairingEngine<Fr = G::ScalarField>, G: DomainGroup> {
    pub(crate) coeffs:    Vec<G>,
    pub(crate) exp:       u32,
    pub(crate) omega:     G::ScalarField,
    pub(crate) omega_inv: G::ScalarField,
    pub(crate) geninv:    G::ScalarField,
    pub(crate) m_inv:     G::ScalarField,
    pub(crate) m:         u64,
    _engine:              PhantomData<E>,
}

impl<E: PairingEngine<Fr = G::ScalarField>, G: DomainGroup> Clone for EvaluationDomain<E, G> {
    fn clone(&self) -> Self {
        Self {
            coeffs: self.coeffs.clone(),
            ..*self
        }
    }
}

impl<E: PairingEngine<Fr = G::ScalarField>, G: DomainGroup> EvaluationDomain<E, G> {
    pub fn compute_m_from_num_coeffs(num_coeffs: usize) -> Option<usize> {
        let m = num_coeffs.next_power_of_two();
        let exp = m.trailing_zeros(); // exp = log_2(m)
        if exp < <G::ScalarField as PrimeField>::Params::TWO_ADICITY {
            Some(m)
        } else {
            None
        }
    }

    pub fn m(&self) -> usize {
        self.m as usize
    }

    pub fn as_ref(&self) -> &[G] {
        &self.coeffs
    }

    pub fn as_mut(&mut self) -> &mut [G] {
        &mut self.coeffs
    }

    pub fn into_coeffs(self) -> Vec<G> {
        self.coeffs
    }

    pub fn from_coeffs(mut coeffs: Vec<G>) -> Option<EvaluationDomain<E, G>> {
        // Compute the size of our evaluation domain
        let m = coeffs.len().next_power_of_two() as u64;
        let exp = m.trailing_zeros(); // exp = log_2(m)

        if exp >= <G::ScalarField as PrimeField>::Params::TWO_ADICITY {
            return None;
        }

        // Compute omega, the 2^exp primitive root of unity
        let mut omega = G::ScalarField::root_of_unity();
        for _ in exp..<G::ScalarField as PrimeField>::Params::TWO_ADICITY {
            omega.square_in_place();
        }

        // Extend the coeffs vector with zeroes if necessary
        coeffs.resize(m as usize, G::group_zero());

        let m_bigint = <<E as PairingEngine>::Fr as PrimeField>::BigInt::from(m);
        let m_inv = G::ScalarField::from_repr(m_bigint).inverse().unwrap();

        Some(EvaluationDomain {
            coeffs,
            exp,
            omega,
            omega_inv: omega.inverse().unwrap(),
            geninv: G::ScalarField::multiplicative_generator()
                .inverse()
                .unwrap(),
            m_inv,
            m,
            _engine: PhantomData,
        })
    }

    pub fn fft(&mut self, worker: &Worker) {
        best_fft(&mut self.coeffs, worker, &self.omega, self.exp);
    }

    pub fn ifft(&mut self, worker: &Worker) {
        best_fft(&mut self.coeffs, worker, &self.omega_inv, self.exp);

        worker.scope(self.coeffs.len(), |scope, chunk| {
            let m_inv = self.m_inv;

            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v {
                        v.group_mul_assign(&m_inv);
                    }
                });
            }
        });
    }

    pub fn distribute_powers(&mut self, worker: &Worker, g: G::ScalarField) {
        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (i, v) in self.coeffs.chunks_mut(chunk).enumerate() {
                scope.spawn(move |_| {
                    let mut u = g.pow(&[(i * chunk) as u64]);
                    for v in v.iter_mut() {
                        v.group_mul_assign(&u);
                        u.mul_assign(&g);
                    }
                });
            }
        });
    }

    pub fn coset_fft(&mut self, worker: &Worker) {
        self.distribute_powers(worker, G::ScalarField::multiplicative_generator());
        self.fft(worker);
    }

    pub fn icoset_fft(&mut self, worker: &Worker) {
        let geninv = self.geninv;

        self.ifft(worker);
        self.distribute_powers(worker, geninv);
    }

    pub fn evaluate_all_lagrange_coefficients(&self, tau: &G::ScalarField) -> Vec<G::ScalarField> {
        // Evaluate all Lagrange polynomials
        let m = self.m as usize;
        let t_m = tau.pow(&[self.m]);
        let one = G::ScalarField::one();
        if t_m == one {
            let mut u = vec![G::ScalarField::zero(); m];
            let mut omega_i = one;
            for i in 0..m {
                if omega_i == *tau {
                    u[i] = one;
                    break;
                }
                omega_i.mul_assign(&self.omega);
            }
            u
        } else {
            use crate::fields::batch_inversion;
            use rayon::prelude::*;

            let mut l = (t_m - &one) * &self.m_inv;
            let mut r = one;
            let mut u = vec![G::ScalarField::zero(); m];
            let mut ls = vec![G::ScalarField::zero(); m];
            for i in 0..m {
                u[i] = *tau - &r;
                ls[i] = l;
                l.mul_assign(&self.omega);
                r.mul_assign(&self.omega);
            }

            batch_inversion(u.as_mut_slice());
            u.par_iter_mut().zip(ls).for_each(|(tau_minus_r, l)| {
                *tau_minus_r = l * tau_minus_r;
            });
            u
        }
    }

    /// This evaluates t(tau) for this domain, which is
    /// tau^m - 1 for these radix-2 domains.
    pub fn z(&self, tau: &G::ScalarField) -> G::ScalarField {
        let mut tmp = tau.pow(&[self.coeffs.len() as u64]);
        tmp.sub_assign(&G::ScalarField::one());

        tmp
    }

    /// The target polynomial is the zero polynomial in our
    /// evaluation domain, so we must perform division over
    /// a coset.
    pub fn divide_by_z_on_coset(&mut self, worker: &Worker) {
        let i = self
            .z(&G::ScalarField::multiplicative_generator())
            .inverse()
            .unwrap();

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for v in self.coeffs.chunks_mut(chunk) {
                scope.spawn(move |_| {
                    for v in v {
                        v.group_mul_assign(&i);
                    }
                });
            }
        });
    }

    /// Perform O(n) multiplication of two polynomials in the domain.
    pub fn mul_assign(&mut self, worker: &Worker, other: &EvaluationDomain<E, Scalar<E>>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.group_mul_assign(&b.0);
                    }
                });
            }
        });
    }

    /// Perform O(n) subtraction of one polynomial from another in the domain.
    pub fn sub_assign(&mut self, worker: &Worker, other: &EvaluationDomain<E, G>) {
        assert_eq!(self.coeffs.len(), other.coeffs.len());

        worker.scope(self.coeffs.len(), |scope, chunk| {
            for (a, b) in self
                .coeffs
                .chunks_mut(chunk)
                .zip(other.coeffs.chunks(chunk))
            {
                scope.spawn(move |_| {
                    for (a, b) in a.iter_mut().zip(b.iter()) {
                        a.group_sub_assign(&b);
                    }
                });
            }
        });
    }
}

pub trait DomainGroup: Sized + Copy + Clone + Send + Sync {
    type ScalarField: PrimeField + SquareRootField + Into<<Self::ScalarField as PrimeField>::BigInt>;
    fn group_zero() -> Self;
    fn group_mul_assign(&mut self, by: &Self::ScalarField);
    fn group_add_assign(&mut self, other: &Self);
    fn group_sub_assign(&mut self, other: &Self);
}

pub struct Point<G: ProjectiveCurve>(pub G);

impl<G: ProjectiveCurve> PartialEq for Point<G> {
    fn eq(&self, other: &Point<G>) -> bool {
        self.0 == other.0
    }
}

impl<G: ProjectiveCurve> Copy for Point<G> {}

impl<G: ProjectiveCurve> Clone for Point<G> {
    fn clone(&self) -> Point<G> {
        *self
    }
}

impl<G: ProjectiveCurve> DomainGroup for Point<G> {
    type ScalarField = G::ScalarField;

    fn group_zero() -> Self {
        Point(G::zero())
    }

    fn group_mul_assign(&mut self, by: &Self::ScalarField) {
        self.0.mul_assign(by.into_repr());
    }

    fn group_add_assign(&mut self, other: &Self) {
        self.0.add_assign(&other.0);
    }

    fn group_sub_assign(&mut self, other: &Self) {
        self.0.sub_assign(&other.0);
    }
}

pub struct Scalar<E: PairingEngine>(pub E::Fr);

impl<E: PairingEngine> PartialEq for Scalar<E> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl<E: PairingEngine> Copy for Scalar<E> {}

impl<E: PairingEngine> Clone for Scalar<E> {
    fn clone(&self) -> Scalar<E> {
        *self
    }
}

impl<E: PairingEngine> DomainGroup for Scalar<E> {
    type ScalarField = E::Fr;

    fn group_zero() -> Self {
        Scalar(E::Fr::zero())
    }
    fn group_mul_assign(&mut self, by: &Self::ScalarField) {
        self.0.mul_assign(by);
    }
    fn group_add_assign(&mut self, other: &Self) {
        self.0.add_assign(&other.0);
    }
    fn group_sub_assign(&mut self, other: &Self) {
        self.0.sub_assign(&other.0);
    }
}

fn best_fft<G: DomainGroup>(a: &mut [G], worker: &Worker, omega: &G::ScalarField, log_n: u32) {
    let log_cpus = worker.log_num_cpus();

    if log_n <= log_cpus {
        serial_fft(a, omega, log_n);
    } else {
        parallel_fft(a, worker, omega, log_n, log_cpus);
    }
}

pub(crate) fn serial_fft<G: DomainGroup>(a: &mut [G], omega: &G::ScalarField, log_n: u32) {
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
            let mut w = G::ScalarField::one();
            for j in 0..m {
                let mut t = a[(k + j + m) as usize];
                t.group_mul_assign(&w);
                let mut tmp = a[(k + j) as usize];
                tmp.group_sub_assign(&t);
                a[(k + j + m) as usize] = tmp;
                a[(k + j) as usize].group_add_assign(&t);
                w.mul_assign(&w_m);
            }

            k += 2 * m;
        }

        m *= 2;
    }
}

pub(crate) fn parallel_fft<G: DomainGroup>(
    a: &mut [G],
    worker: &Worker,
    omega: &G::ScalarField,
    log_n: u32,
    log_cpus: u32,
) {
    assert!(log_n >= log_cpus);

    let num_cpus = 1 << log_cpus;
    let log_new_n = log_n - log_cpus;
    let mut tmp = vec![vec![G::group_zero(); 1 << log_new_n]; num_cpus];
    let new_omega = omega.pow(&[num_cpus as u64]);

    worker.scope(0, |scope, _| {
        let a = &*a;

        for (j, tmp) in tmp.iter_mut().enumerate() {
            scope.spawn(move |_| {
                // Shuffle into a sub-FFT
                let omega_j = omega.pow(&[j as u64]);
                let omega_step = omega.pow(&[(j as u64) << log_new_n]);

                let mut elt = G::ScalarField::one();
                for i in 0..(1 << log_new_n) {
                    for s in 0..num_cpus {
                        let idx = (i + (s << log_new_n)) % (1 << log_n);
                        let mut t = a[idx];
                        t.group_mul_assign(&elt);
                        tmp[i].group_add_assign(&t);
                        elt.mul_assign(&omega_step);
                    }
                    elt.mul_assign(&omega_j);
                }

                // Perform sub-FFT
                serial_fft(tmp, &new_omega, log_new_n);
            });
        }
    });

    // TODO: does this hurt or help?
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
