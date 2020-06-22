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
    0x0a8fc7bc1a89fa86,
    0xa7d9d786e9e48627,
    0xee6158b465bea369,
    0x14a0ff6d2f874519,
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
/// Baby-JubJub's curve equation: 168700x² + y² = 1 + 168696x²y²
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
    /// COEFF_A = 168700
    #[rustfmt::skip]
    const COEFF_A: Fq = field_new!(Fq, BigInteger256([
        0x95accf61fff261e0,
        0x24780d659df7d378,
        0xe0ac11b07e906ae8,
        0xf35db2216d3def3,
    ]));

    /// COEFF_D = 168696
    #[rustfmt::skip]
    const COEFF_D: Fq = field_new!(Fq, BigInteger256([
        0x2735f484aff261f5,
        0x70ba1b579a2e0f63,
        0xff41c9a91e2caa8c,
        0x7704a8e8fe6025f,
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

    /*
    /// Multiplication by `a` is simply negation here.
    #[inline(always)]
    fn mul_by_a(elem: &Self::BaseField) -> Self::BaseField {
        *elem
    }
    */
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
    /// COEFF_B = 1
    #[rustfmt::skip]
    const COEFF_B: Fq = field_new!(Fq, BigInteger256([
        12436184717236109307u64,
        3962172157175319849u64,
        7381016538464732718u64,
        1011752739694698287u64,
    ]));

    type TEModelParameters = EdwardsParameters;
}
