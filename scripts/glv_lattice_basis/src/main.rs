extern crate algebra_core;
extern crate algebra;
extern crate num_traits;

use algebra::bw6_761::Fr;
use algebra_core::{biginteger::{BigInteger384, BigInteger}, fields::PrimeField};
mod arithmetic;
use crate::arithmetic::div_with_remainder;
use std::ops::Neg;
use num_traits::Zero;

fn main() {
    let n = BigInteger384([
        0x8508c00000000001,
        0x170b5d4430000000,
        0x1ef3622fba094800,
        0x1a22d9f300f5138f,
        0xc63b05c06ca1493b,
        0x1ae3a4617c510ea,
    ]);
    let lambda = BigInteger384([
        0x8508c00000000001,
        0x452217cc90000000,
        0xc5ed1347970dec00,
        0x619aaf7d34594aab,
        0x9b3af05dd14f6ec,
        0x0
    ]);
    // println!("{:?}",);

    let vecs = get_lattice_basis::<Fr>(n, lambda);
    for vec in [vecs.0, vecs.1].iter() {
        println!("vec: {:?}", vec);
        let (s1, (flag, t1)) = vec;
        debug_assert_eq!(recompose_integer(
            Fr::from_repr(*s1).unwrap(),
            if !flag {
                Fr::from_repr(*t1).unwrap()
            } else {
                Fr::from_repr(*t1).unwrap().neg()
            },
            Fr::from_repr(lambda).unwrap()
        ),
        Fr::zero());
    }
}

// We work on arrays of size 3
// We assume that |E(F_q)| < R = 2^{ceil(limbs/2) * 64}
fn get_lattice_basis<F: PrimeField>(n: F::BigInt, lambda: F::BigInt) -> ((F::BigInt, (bool, F::BigInt)), (F::BigInt, (bool, F::BigInt)))
{
    let mut r: Vec<F::BigInt> = vec![n, lambda, n];
    let one = F::from(F::BigInt::from(1));
    let zero = F::from(F::BigInt::from(0));
    let mut t = [zero, one, zero];
    let max_num_bits_lattice = F::BigInt::from_slice(F::characteristic()).num_bits() / 2 + 1;

    let sqrt_n = as_f64(n.as_ref()).sqrt();

    let mut i = 0;
    // While r_i >= sqrt(n), we then return the vectors (r_i, t_i), (r_i+1, t_i+1)
    while as_f64(r[(i + 1) % 3].as_ref()) >= sqrt_n {
        let (q, rem): (F::BigInt, F::BigInt) = div_with_remainder::<F::BigInt>(r[i % 3], r[(i + 1) % 3]);
        r[(i + 2) % 3] = rem;
        let int_q = F::from(q);
        t[(i + 2) % 3] = t[i % 3] - int_q * (t[(i + 1) % 3]);
        i += 1;
    }

    // we do a conversion from the fields into
    let (neg_flag1, t1) = if t[(i + 1) % 3].into_repr().num_bits() <= max_num_bits_lattice {
        (false, t[(i + 1) % 3].into_repr())
    } else {
        (true, t[(i + 1) % 3].neg().into_repr())
    };
    let (neg_flag2, t2) = if t[(i + 2) % 3].into_repr().num_bits() <= max_num_bits_lattice{
        (false, t[(i + 2) % 3].into_repr())
    } else {
        (true, t[(i + 2) % 3].neg().into_repr())
    };
    let vec_1 = (r[(i + 1) % 3], (neg_flag1, t1));
    let vec_2 = (r[(i + 2) % 3], (neg_flag2, t2));

    (vec_1, vec_2)
}

fn recompose_integer<F: PrimeField>(k1: F, k2: F, lambda: F) -> F {
    k1 + &(k2 * &lambda)
}

fn as_f64(bigint_ref: &[u64]) -> f64 {
    let mut n_float: f64 = 0.0;
    for (i, limb) in bigint_ref.iter().enumerate() {
        n_float += (*limb as f64) * 2f64.powf((i as f64) * 64f64)
    }
    n_float
}
//
// struct iBigInteger<BigInt: BigInteger> {
//     value: BigInt,
//     neg: bool,
// }
//
// impl iBigInteger {}
//
// impl<BigInt: BigInteger> Mul for iBigInteger<BigInt> {
//     fn mul_assign(&mut self, other: &Self) {
//         self.value *= other.value;
//         match (self.neg, other.neg) {
//             (true, true) => self.neg(),
//             (false, true) => self.neg(),
//             _ => (),
//         }
//     }
// }
//
// impl<BigInt: BigInteger> Neg for iBigInteger<BigInt> {
//     fn neg(&mut self) {
//         if self.neg {
//             self.neg = false;
//         } else {
//             self.neg = true;
//         }
//     }
// }
//
// impl<BigInt: BigInteger> Sub for iBigInteger<BigInt> {
//     fn sub_assign(&mut self, other: &Self) {
//         self.add_nocarry(other.neg());
//     }
// }
//
// impl<BigInt: BigInteger> Add for iBigInteger<BigInt> {
//     fn add_assign(&mut self, other: &Self) {
//         // If operators have the same sign, just add the values
//         if self.neg + other.neg == false {
//             self.value += other.value;
//         } else {
//             if self.value > other.value {
//                 self.sub_noborrow(other);
//             } else {
//                 let mut tmp = other.clone();
//                 tmp.sub_noborrow(self.value);
//                 self.value = tmp;
//                 self.neg();
//             }
//         }
//     }
// }
//
// impl<BigInt: BigInteger> From<BigInt> for iBigInteger<BigInt> {
//     #[inline]
//     fn from(val: BigInt) -> iBigInteger<BigInt> {
//         iBigInteger::<BigInt>{
//             value: val,
//             neg: false,
//         }
//     }
// }
