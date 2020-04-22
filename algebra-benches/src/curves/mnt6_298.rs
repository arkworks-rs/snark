use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::ops::{AddAssign, MulAssign, SubAssign};

use algebra::{mnt6_298::{Fq6, fq3::Fq3, fq::Fq, fr::Fr, G1Affine, G1Projective as G1,
    MNT6_298, G2Affine, G2Projective as G2, Parameters},
    PairingEngine,
    Field, SquareRootField, PrimeField, ProjectiveCurve, UniformRand,
    BigInteger, biginteger::{BigInteger320 as FqRepr},
    mnt6::{G1Prepared, G2Prepared}};

ec_bench!();
f_bench!(1, Fq3, Fq3, fq3);
f_bench!(2, Fq6, Fq6, fq6);
f_bench!(Fq, Fq, FqRepr, FqRepr, fq);
pairing_bench!(MNT6_298, Fq6, prepared_v);
