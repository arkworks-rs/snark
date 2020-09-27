extern crate algebra;
extern crate algebra_core;
extern crate num_traits;

use algebra::bw6_761::{Fq, Fr};
use algebra_core::{
    biginteger::{BigInteger, BigInteger384, BigInteger768},
    fields::PrimeField,
};
use glv_lattice_basis::*;
use num_traits::Zero;
use std::ops::Neg;

fn main() {
    let _omega_g1 = BigInteger768([
        0x962140000000002a,
        0xc547ba8a4000002f,
        0xb6290012d96f8819,
        0xf2f082d4dcb5e37c,
        0xc65759fc45183151,
        0x8e0a235a0a398300,
        0xab5e57926fa70184,
        0xee4a737f73b6f952,
        0x2d17be416c5e4426,
        0x6c1f31e53bd9603c,
        0xaa846c61024e4cca,
        0x531dc16c6ecd27,
    ]);
    let _omega_g2 = BigInteger768([
        0x5e7bc00000000060,
        0x214983de30000053,
        0x5fe3f89c11811c1e,
        0xa5b093ed79b1c57b,
        0xab8579e02ed3cddc,
        0xf87fa59308c07a8f,
        0x5870636cb60d217f,
        0x823132b971cdefc6,
        0x256ab7ae14297a1a,
        0x4d06e68545f7e64c,
        0x27035cdf02acb274,
        0xcfca638f1500e3,
    ]);
    println!(
        "const OMEGA: Self::BaseField = {:?};",
        Fq::from_repr(_omega_g2).unwrap()
    );
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
        0x0,
    ]);
    println!(
        "const LAMBDA: Self::ScalarField = {:?};",
        Fr::from_repr(lambda).unwrap()
    );

    let vecs = get_lattice_basis::<Fr>(n, lambda);

    for (i, vec) in [vecs.0, vecs.1].iter().enumerate() {
        // println!("vec: {:?}", vec);
        let (s1, (flag, t1)) = vec;

        let mut t1_big = BigInteger768::from_slice(t1.as_ref());
        let n_big = BigInteger768::from_slice(n.as_ref());
        t1_big.muln(BigInteger384::NUM_LIMBS as u32 * 64);
        let (g1_big, _) = div_with_remainder::<BigInteger768>(t1_big, n_big);
        let g1 = BigInteger384::from_slice(g1_big.as_ref());

        println!("/// |round(B{} * R / n)|", i + 1);
        println!(
            "const Q{}: <Self::ScalarField as PrimeField>::BigInt = {:?};",
            ((i + 1) % 2) + 1,
            g1
        );
        println!(
            "const B{}: <Self::ScalarField as PrimeField>::BigInt = {:?};",
            i + 1,
            t1
        );
        println!("const B{}_IS_NEG: bool = {:?};", i + 1, flag);

        debug_assert_eq!(
            recompose_integer(
                Fr::from_repr(*s1).unwrap(),
                if !flag {
                    Fr::from_repr(*t1).unwrap()
                } else {
                    Fr::from_repr(*t1).unwrap().neg()
                },
                Fr::from_repr(lambda).unwrap()
            ),
            Fr::zero()
        );
    }
    println!("const R_BITS: u32 = {:?};", BigInteger384::NUM_LIMBS * 64);
}
