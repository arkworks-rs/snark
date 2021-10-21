use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::ops::{AddAssign, MulAssign, SubAssign};

use algebra::{
    fields::mnt4753::{
        fq::Fq, fq2::Fq2, fr::Fr, Fq4
    },
    curves::{
        mnt4753::{
            MNT4, G1Affine, G1Projective as G1, G2Affine,
            G2Projective as G2, MNT4_753Parameters as Parameters,
        },
        models::mnt4::{G1Prepared, G2Prepared},
    },
    biginteger::BigInteger768 as FqRepr,
    BigInteger, Field, PairingEngine, PrimeField, ProjectiveCurve, SquareRootField, UniformRand,
};

ec_bench!();
f_bench!(1, Fq2, Fq2, fq2);
f_bench!(2, Fq4, Fq4, fq4);
f_bench!(Fq, Fq, FqRepr, FqRepr, fq);
pairing_bench!(MNT4, Fq4, prepared_v);
