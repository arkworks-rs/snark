use crate::mnt6_298::{self, Fq, Fr};
use algebra_core::{
    biginteger::BigInteger320,
    curves::{
        mnt6,
        models::{ModelParameters, SWModelParameters},
    },
    field_new,
};

pub type G1Affine = mnt6::G1Affine<mnt6_298::Parameters>;
pub type G1Projective = mnt6::G1Projective<mnt6_298::Parameters>;
pub type G1Prepared = mnt6::G1Prepared<mnt6_298::Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl ModelParameters for Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl SWModelParameters for Parameters {
    /// COEFF_A =
    #[rustfmt::skip]
    const COEFF_A: Fq = field_new!(Fq, BigInteger320([
        0xb9b2411bfd0eafef,
        0xc61a10fadd9fecbd,
        0x89f128e59811f3fb,
        0x980c0f780adadabb,
        0x9ba1f11320,
    ]));

    /// COEFF_B =
    #[rustfmt::skip]
    const COEFF_B: Fq = field_new!(Fq, BigInteger320([
        0xa94cb16ed8e733b,
        0xe1ed15e8119bae6,
        0xae927592157c8121,
        0x990dbcbc6661cf95,
        0xecff0892ef,
    ]));

    /// COFACTOR = 1
    const COFACTOR: &'static [u64] = &[1];

    /// COFACTOR^(-1) mod r =
    /// 1
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger320([
        1784298994435064924,
        16852041090100268533,
        14258261760832875328,
        2961187778261111191,
        1929014752195,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G1_GENERATOR_X, G1_GENERATOR_Y);
}

/// G1_GENERATOR_X =
#[rustfmt::skip]
pub const G1_GENERATOR_X: Fq = field_new!(Fq, BigInteger320([
    0x1a663562f74e1d24,
    0xc1d1d583fccd1b79,
    0xda077538a9763df2,
    0x70c4a4ea36aa01d9,
    0x86537578a8,
]));

/// G1_GENERATOR_Y =
#[rustfmt::skip]
pub const G1_GENERATOR_Y: Fq = field_new!(Fq, BigInteger320([
    0x7ad5bfd16dcfffb2,
    0x88dd739252215070,
    0x43f137a8b517b339,
    0x9a7fac709a8c463c,
    0x3140fbc3593,
]));
