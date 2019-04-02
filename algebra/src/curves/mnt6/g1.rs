use crate::{
    biginteger::BigInteger320,
    bytes::ToBytes,
    curves::{
        mnt6::{G2Affine, MNT6},
        models::{ModelParameters, SWModelParameters},
        short_weierstrass_projective::{GroupAffine, GroupProjective},
        AffineCurve, PairingCurve, PairingEngine,
    },
    fields::mnt6::{Fq, Fq3, Fq6, Fr},
};
use std::io::{Result as IoResult, Write};

pub type G1Affine = GroupAffine<MNT6G1Parameters>;
pub type G1Projective = GroupProjective<MNT6G1Parameters>;

impl PairingCurve for G1Affine {
    type Engine = MNT6;
    type Prepared = G1Prepared;
    type PairWith = G2Affine;
    type PairingResult = Fq6;

    fn prepare(&self) -> Self::Prepared {
        Self::Prepared::from_affine(self)
    }

    fn pairing_with(&self, other: &Self::PairWith) -> Self::PairingResult {
        MNT6::pairing(*self, *other)
    }
}

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct MNT6G1Parameters;

impl ModelParameters for MNT6G1Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl SWModelParameters for MNT6G1Parameters {
    /// COEFF_A =
    const COEFF_A: Fq = Fq::new(BigInteger320([
        0xb9b2411bfd0eafef,
        0xc61a10fadd9fecbd,
        0x89f128e59811f3fb,
        0x980c0f780adadabb,
        0x9ba1f11320,
    ]));

    /// COEFF_B =
    const COEFF_B: Fq = Fq::new(BigInteger320([
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
    const COFACTOR_INV: Fr = Fr::new(BigInteger320([
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
pub const G1_GENERATOR_X: Fq = Fq::new(BigInteger320([
    0x1a663562f74e1d24,
    0xc1d1d583fccd1b79,
    0xda077538a9763df2,
    0x70c4a4ea36aa01d9,
    0x86537578a8,
]));

/// G1_GENERATOR_Y =
pub const G1_GENERATOR_Y: Fq = Fq::new(BigInteger320([
    0x7ad5bfd16dcfffb2,
    0x88dd739252215070,
    0x43f137a8b517b339,
    0x9a7fac709a8c463c,
    0x3140fbc3593,
]));

#[derive(Eq, PartialEq, Copy, Clone, Debug)]
pub struct G1Prepared {
    pub x:       Fq,
    pub y:       Fq,
    pub x_twist: Fq3,
    pub y_twist: Fq3,
}

impl ToBytes for G1Prepared {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.x.write(&mut writer)?;
        self.y.write(&mut writer)?;
        self.x_twist.write(&mut writer)?;
        self.y_twist.write(&mut writer)
    }
}

impl G1Prepared {
    pub fn from_affine(point: &G1Affine) -> Self {
        MNT6::ate_precompute_g1(&point.into_projective())
    }
}

impl Default for G1Prepared {
    fn default() -> Self {
        Self::from_affine(&G1Affine::prime_subgroup_generator())
    }
}
