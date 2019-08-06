use crate::field_new;
use crate::{
    biginteger::{BigInteger256, BigInteger384},
    curves::{
        bls12::{G1Affine as Bls12G1Affine, G1Prepared, G1Projective as Bls12G1Projective},
        bls12_381::{g2::G2Affine, Bls12_381, Bls12_381Parameters},
        models::{ModelParameters, SWModelParameters},
        PairingCurve, PairingEngine,
    },
    fields::{
        bls12_381::{Fq, Fq12, Fr},
        Field,
    },
};

pub type G1Affine = Bls12G1Affine<Bls12_381Parameters>;
pub type G1Projective = Bls12G1Projective<Bls12_381Parameters>;

impl PairingCurve for G1Affine {
    type Engine = Bls12_381;
    type Prepared = G1Prepared<Bls12_381Parameters>;
    type PairWith = G2Affine;
    type PairingResult = Fq12;

    fn prepare(&self) -> Self::Prepared {
        Self::Prepared::from_affine(*self)
    }

    fn pairing_with(&self, other: &Self::PairWith) -> Self::PairingResult {
        Bls12_381::pairing(*self, *other)
    }
}

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct Bls12_381G1Parameters;

impl ModelParameters for Bls12_381G1Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl SWModelParameters for Bls12_381G1Parameters {
    /// COEFF_A = 0
    const COEFF_A: Fq = field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0]));

    /// COEFF_B = 4
    const COEFF_B: Fq = field_new!(Fq, BigInteger384([
        0xaa270000000cfff3,
        0x53cc0032fc34000a,
        0x478fe97a6b0a807f,
        0xb1d37ebee6ba24d7,
        0x8ec9733bbf78ab2f,
        0x9d645513d83de7e,
    ]));

    /// COFACTOR = (x - 1)^2 / 3  = 76329603384216526031706109802092473003
    const COFACTOR: &'static [u64] = &[0x8c00aaab0000aaab, 0x396c8c005555e156];

    /// COFACTOR_INV = COFACTOR^{-1} mod r
    /// = 52435875175126190458656871551744051925719901746859129887267498875565241663483
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger256([
        288839107172787499,
        1152722415086798946,
        2612889808468387987,
        5124657601728438008,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G1_GENERATOR_X, G1_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }
}

/// G1_GENERATOR_X =
/// 3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507
pub const G1_GENERATOR_X: Fq = field_new!(Fq, BigInteger384([
    0x5cb38790fd530c16,
    0x7817fc679976fff5,
    0x154f95c7143ba1c1,
    0xf0ae6acdf3d0e747,
    0xedce6ecc21dbf440,
    0x120177419e0bfb75,
]));

/// G1_GENERATOR_Y =
/// 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569
pub const G1_GENERATOR_Y: Fq = field_new!(Fq, BigInteger384([
    0xbaac93d50ce72271,
    0x8c22631a7918fd8e,
    0xdd595f13570725ce,
    0x51ac582950405194,
    0xe1c8c3fad0059c0,
    0xbbc3efc5008a26a,
]));
