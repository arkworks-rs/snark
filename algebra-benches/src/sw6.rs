use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::ops::{AddAssign, MulAssign, SubAssign};

use algebra::{sw6::{Fq6, fq3::Fq3, fq::Fq, fr::Fr, G1Affine, G1Projective as G1,
    SW6, G2Affine, G2Projective as G2},
    PairingEngine,
    Field, SquareRootField, PrimeField, ProjectiveCurve, UniformRand,
    BigInteger, biginteger::{BigInteger832 as FqRepr, BigInteger384 as FrRepr}
    };

ec_bench!();
f_bench!(1, Fq3, Fq3, fq3);
f_bench!(2, Fq6, Fq6, fq6);
f_bench!(Fq, Fq, FqRepr, FqRepr, fq);
f_bench!(Fr, Fr, FrRepr, FrRepr, fr);
pairing_bench!(SW6, Fq6, affine_v);
