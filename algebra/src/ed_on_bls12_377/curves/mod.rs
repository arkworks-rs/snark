use crate::ed_on_bls12_377::{fq::Fq, fr::Fr};
use algebra_core::{
    biginteger::BigInteger256,
    curves::{
        models::{ModelParameters, MontgomeryModelParameters, TEModelParameters},
        twisted_edwards_extended::{GroupAffine, GroupProjective},
    },
    field_new,
};

#[cfg(test)]
mod tests;

pub type EdwardsAffine = GroupAffine<EdwardsParameters>;
pub type EdwardsProjective = GroupProjective<EdwardsParameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct EdwardsParameters;

impl ModelParameters for EdwardsParameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl TEModelParameters for EdwardsParameters {
    /// COEFF_A = -1
    #[rustfmt::skip]
    const COEFF_A: Fq = field_new!(Fq, BigInteger256([
        0x8cf500000000000e,
        0xe75281ef6000000e,
        0x49dc37a90b0ba012,
        0x55f8b2c6e710ab9,
    ]));

    /// COEFF_D = 3021
    #[rustfmt::skip]
    const COEFF_D: Fq = field_new!(Fq, BigInteger256([
        0xd047ffffffff5e30,
        0xf0a91026ffff57d2,
        0x9013f560d102582,
        0x9fd242ca7be5700,
    ]));

    /// COFACTOR = 4
    const COFACTOR: &'static [u64] = &[4];

    /// COFACTOR_INV =
    /// 527778859339273151515551558673846658209717731602102048798421311598680340096
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger256([
        10836190823041854989,
        14880086764632731920,
        5023208332782666747,
        239524813690824359,
    ]));

    /// Generated randomly
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) = (GENERATOR_X, GENERATOR_Y);

    type MontgomeryModelParameters = EdwardsParameters;

    /// Multiplication by `a` is just negation.
    /// Is `a` 1 or -1?
    #[inline(always)]
    fn mul_by_a(elem: &Self::BaseField) -> Self::BaseField {
        -*elem
    }
}

impl MontgomeryModelParameters for EdwardsParameters {
    /// COEFF_A = 0x8D26E3FADA9010A26949031ECE3971B93952AD84D4753DDEDB748DA37E8F552
    #[rustfmt::skip]
    const COEFF_A: Fq = field_new!(Fq, BigInteger256([
        13800168384327121454u64,
        6841573379969807446u64,
        12529593083398462246u64,
        853978956621483129u64,
    ]));
    /// COEFF_B = 0x9D8F71EEC83A44C3A1FBCEC6F5418E5C6154C2682B8AC231C5A3725C8170AAD
    #[rustfmt::skip]
    const COEFF_B: Fq = field_new!(Fq, BigInteger256([
        7239382437352637935u64,
        14509846070439283655u64,
        5083066350480839936u64,
        1265663645916442191u64,
    ]));

    type TEModelParameters = EdwardsParameters;
}

/// GENERATOR_X =
/// 7810607721416582242904415504650443951498042435501746664987470571546413371306
#[rustfmt::skip]
const GENERATOR_X: Fq = field_new!(Fq, BigInteger256([
    0x5bbc9878d817221d,
    0xd2b03489424e720,
    0x6b66f128c16bb3c9,
    0xdd3bff78733576d,
]));

/// GENERATOR_Y =
/// 1867362672570137759132108893390349941423731440336755218616442213142473202417
#[rustfmt::skip]
const GENERATOR_Y: Fq = field_new!(Fq, BigInteger256([
    0x471517ae5e5e979e,
    0xd9c97f6a73a7ff83,
    0x85a95b45a5494402,
    0xfad27c9b545b1f0,
]));
