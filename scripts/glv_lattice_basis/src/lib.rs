extern crate algebra;
extern crate algebra_core;
extern crate num_traits;

mod arithmetic;

use algebra::BigInteger;
use algebra_core::fields::PrimeField;
pub use arithmetic::*;

// We work on arrays of size 3
// We assume that |E(F_q)| < R = 2^{ceil(limbs/2) * 64}
pub fn get_lattice_basis<F: PrimeField>(
    n: F::BigInt,
    lambda: F::BigInt,
) -> (
    (F::BigInt, (bool, F::BigInt)),
    (F::BigInt, (bool, F::BigInt)),
) {
    let mut r = [n, lambda, n];
    let one = F::one();
    let zero = F::zero();
    let mut t: [F; 3] = [zero, one, zero];
    let max_num_bits_lattice = (F::BigInt::from_slice(F::characteristic()).num_bits() - 1) / 2 + 1;

    let sqrt_n = as_f64(n.as_ref()).sqrt();

    println!("Log sqrtn: {}", sqrt_n.log2());

    let mut i = 0;
    // While r_i >= sqrt(n), we perform the extended euclidean algorithm so that
    // si*n + ti*lambda = ri then return the vectors (r_i, (sign(t_i), |t_i|)),
    // (r_i+1, (sign(t_i+1), |t_i+1|)) Notice this makes ri + (-ti)*lambda = 0
    // mod n, which is what we desire for our short lattice basis
    while as_f64(r[i % 3].as_ref()) >= sqrt_n {
        // while i < 20 {
        let (q, rem): (F::BigInt, F::BigInt) =
            div_with_remainder::<F::BigInt>(r[i % 3], r[(i + 1) % 3]);
        r[(i + 2) % 3] = rem;
        let int_q = F::from_repr(q).unwrap();
        t[(i + 2) % 3] = t[i % 3] - int_q * (t[(i + 1) % 3]);

        i += 1;
    }
    let just_computed = (i + 1) % 3;
    let (neg_flag1, t1) = if t[just_computed].into_repr().num_bits() <= max_num_bits_lattice {
        (false, t[just_computed].into_repr())
    } else {
        (true, t[just_computed].neg().into_repr())
    };
    let vec_1 = (r[just_computed], (neg_flag1, t1));

    let prev = i % 3;
    let (neg_flag2, t2) = if t[prev].into_repr().num_bits() <= max_num_bits_lattice {
        (false, t[prev].into_repr())
    } else {
        (true, t[prev].neg().into_repr())
    };
    let vec_2 = (r[prev], (neg_flag2, t2));

    (vec_1, vec_2)
}

pub fn recompose_integer<F: PrimeField>(k1: F, k2: F, lambda: F) -> F {
    k1 - &(k2 * &lambda)
}

fn as_f64(bigint_ref: &[u64]) -> f64 {
    let mut n_float: f64 = 0.0;
    for (i, limb) in bigint_ref.iter().enumerate() {
        n_float += (*limb as f64) * 2f64.powf((i as f64) * 64f64)
    }
    n_float
}
