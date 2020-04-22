use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::ops::{AddAssign, MulAssign, SubAssign};

use algebra::{bls12_377::{Fq12, fq2::Fq2, fq::Fq, fr::Fr, G1Affine, G1Projective as G1,
    Bls12_377, G2Affine, G2Projective as G2, Parameters},
    PairingEngine,
    BigInteger, Field, SquareRootField, PrimeField, ProjectiveCurve, UniformRand,
    biginteger::{BigInteger384 as FqRepr, BigInteger256 as FrRepr},
    bls12::{G1Prepared, G2Prepared},};

ec_bench!();
f_bench!(1, Fq2, Fq2, fq2);
f_bench!(2, Fq12, Fq12, fq12);
f_bench!(Fq, Fq, FqRepr, FqRepr, fq);
f_bench!(Fr, Fr, FrRepr, FrRepr, fr);
pairing_bench!(Bls12_377, Fq12, prepared_v);
