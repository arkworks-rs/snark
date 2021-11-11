use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::ops::{AddAssign, MulAssign, SubAssign};

use algebra::{
    biginteger::BigInteger768 as FqRepr,
    curves::{
        mnt6753::{
            G1Affine, G1Projective as G1, G2Affine, G2Projective as G2,
            MNT6_753Parameters as Parameters, MNT6,
        },
        models::mnt6::{G1Prepared, G2Prepared},
    },
    fields::mnt6753::{fq::Fq, fq3::Fq3, fr::Fr, Fq6},
    BigInteger, Field, PairingEngine, PrimeField, ProjectiveCurve, SquareRootField, UniformRand,
};

ec_bench!();
f_bench!(1, Fq3, Fq3, fq3);
f_bench!(2, Fq6, Fq6, fq6);
f_bench!(Fq, Fq, FqRepr, FqRepr, fq);
pairing_bench!(MNT6, Fq6, prepared_v);
