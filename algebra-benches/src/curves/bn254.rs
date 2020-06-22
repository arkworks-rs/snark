use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::ops::{AddAssign, MulAssign, SubAssign};

use algebra::{
    biginteger::{BigInteger256 as FrRepr, BigInteger256 as FqRepr},
    bn::{G1Prepared, G2Prepared},
    bn254::{
        fq::Fq, fq2::Fq2, fr::Fr, Bn254, Fq12, G1Affine, G1Projective as G1, G2Affine,
        G2Projective as G2, Parameters,
    },
    BigInteger, Field, PairingEngine, PrimeField, ProjectiveCurve, SquareRootField, UniformRand,
};

ec_bench!();
f_bench!(1, Fq2, Fq2, fq2);
f_bench!(2, Fq12, Fq12, fq12);
f_bench!(Fq, Fq, FqRepr, FqRepr, fq);
f_bench!(Fr, Fr, FrRepr, FrRepr, fr);
pairing_bench!(Bn254, Fq12, prepared_v);
