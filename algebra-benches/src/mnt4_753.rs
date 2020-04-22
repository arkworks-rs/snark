use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::ops::{AddAssign, MulAssign, SubAssign};

use algebra::{mnt4_753::{Fq4, fq2::Fq2, fq::Fq, fr::Fr, G1Affine, G1Projective as G1,
    MNT4_753, G2Affine, G2Projective as G2, Parameters},
    PairingEngine,
    Field, SquareRootField, PrimeField, ProjectiveCurve, UniformRand,
    BigInteger, biginteger::{BigInteger768 as FqRepr},
    mnt4::{G1Prepared, G2Prepared}};

ec_bench!();
f_bench!(1, Fq2, Fq2, fq2);
f_bench!(2, Fq4, Fq4, fq4);
f_bench!(Fq, Fq, FqRepr, FqRepr, fq);
pairing_bench!(MNT4_753, Fq4, prepared_v);
