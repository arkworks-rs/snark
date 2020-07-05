use crate::ed_on_bn254::{Fq, Fr};
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

#[rustfmt::skip]
const GENERATOR_X: Fq = field_new!(Fq, BigInteger256([
    0x3db6612c2863cc99,
    0x8a9e4521b36347dc,
    0x310a1a625c16a534,
    0x23ceae2710df4a14,
]));
#[rustfmt::skip]
const GENERATOR_Y: Fq = field_new!(Fq, BigInteger256([
    0xb83342d20d0201aa,
    0x2ffef2f7cdcfeac7,
    0xbfa79a9425a6e625,
    0xdfb859dc3a44b70,
]));

/// `Baby-JubJub` is a twisted Edwards curve. These curves have equations of the
/// form: ax² + y² = 1 + dx²y².
/// over some base finite field Fq.
///
/// Baby-JubJub's curve equation: x² + y² = 1 + (168696/168700)x²y²
///
/// q = 21888242871839275222246405745257275088548364400416034343698204186575808495617
///
#[derive(Clone, Default, PartialEq, Eq)]
pub struct EdwardsParameters;

impl ModelParameters for EdwardsParameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl TEModelParameters for EdwardsParameters {
    /// COEFF_A = 1
    #[rustfmt::skip]
    const COEFF_A: Fq = field_new!(Fq, BigInteger256([
        0xac96341c4ffffffb,
        0x36fc76959f60cd29,
        0x666ea36f7879462e,
        0xe0a77c19a07df2f,
    ]));

    #[inline(always)]
    fn mul_by_a(elem: &Self::BaseField) -> Self::BaseField {
        *elem
    }

    /// COEFF_D = 168696/168700 mod q
    ///         = 9706598848417545097372247223557719406784115219466060233080913168975159366771
    #[rustfmt::skip]
    const COEFF_D: Fq = field_new!(Fq, BigInteger256([
        0xe7a66d1d9fb08e74,
        0xd775bbd5e17629dc,
        0x70ccd097286ef1e7,
        0x45809398fdf98,
    ]));

    /// COFACTOR = 8
    const COFACTOR: &'static [u64] = &[8];

    /// COFACTOR^(-1) mod r =
    /// 2394026564107420727433200628387514462817212225638746351800188703329891451411
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger256([
        0xfac308b2e25a3d4b,
        0xa7c55b66e25b59cb,
        0xeccdd46def0f28c5,
        0x1c14ef83340fbe5,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (GENERATOR_X, GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) = (GENERATOR_X, GENERATOR_Y);

    type MontgomeryModelParameters = EdwardsParameters;
}

impl MontgomeryModelParameters for EdwardsParameters {
    /// COEFF_A = 168698
    #[rustfmt::skip]
    const COEFF_A: Fq = field_new!(Fq, BigInteger256([
        9251058552732279275u64,
        16047179255329565110u64,
        14708493084570629864u64,
        2559515811206512830u64,
    ]));
    /// COEFF_B = 168700
    #[rustfmt::skip]
    const COEFF_B: Fq = field_new!(Fq, BigInteger256([
        10785223227458347488u64,
        2627865112663806840u64,
        16189334210225400552u64,
        1096023023792938739u64,
    ]));

    type TEModelParameters = EdwardsParameters;
}
