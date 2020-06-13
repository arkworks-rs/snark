use crate::ed_on_bls12_381::{Fq, Fr};
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
    14080349899812819339,
    4104857150246327429,
    8293216003873356624,
    7400363483732984990,
]));
#[rustfmt::skip]
const GENERATOR_Y: Fq = field_new!(Fq, BigInteger256([
    13388310974700241893,
    7654361511478576605,
    8037907163910805792,
    5188938133920569885,
]));

/// `JubJub` is a twisted Edwards curve. These curves have equations of the
/// form: ax² + y² = 1 - dx²y².
/// over some base finite field Fq.
///
/// JubJub's curve equation: -x² + y² = 1 - (10240/10241)x²y²
///
/// q = 52435875175126190479447740508185965837690552500527637822603658699938581184513.
///
/// a = -1.
/// d = (10240/10241) mod q
///   = 19257038036680949359750312669786877991949435402254120286184196891950884077233.
///
/// Sage script to calculate these:
///
/// ```text
/// q = 52435875175126190479447740508185965837690552500527637822603658699938581184513
/// Fq = GF(q)
/// d = -(Fq(10240)/Fq(10241))
/// ```
/// These parameters and the sage script obtained from:
/// <https://github.com/zcash/zcash/issues/2230#issuecomment-317182190>
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
        18446744060824649731,
        18102478225614246908,
        11073656695919314959,
        6613806504683796440,
    ]));

    /// COEFF_D = (10240/10241) mod q
    #[rustfmt::skip]
    const COEFF_D: Fq = field_new!(Fq, BigInteger256([
        3049539848285517488,
        18189135023605205683,
        8793554888777148625,
        6339087681201251886,
    ]));

    /// COFACTOR = 8
    const COFACTOR: &'static [u64] = &[8];

    /// COFACTOR^(-1) mod r =
    /// 819310549611346726241370945440405716213240158234039660170669895299022906775
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger256([
        6832491983681988242,
        12911748493335322362,
        17523939349049608702,
        217463794347581613,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (GENERATOR_X, GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) = (GENERATOR_X, GENERATOR_Y);

    type MontgomeryModelParameters = EdwardsParameters;

    /// Multiplication by `a` is simply negation here.
    #[inline(always)]
    fn mul_by_a(elem: &Self::BaseField) -> Self::BaseField {
        -(*elem)
    }
}

impl MontgomeryModelParameters for EdwardsParameters {
    /// COEFF_A = 0xA002
    #[rustfmt::skip]
    const COEFF_A: Fq = field_new!(Fq, BigInteger256([
        388496971701930u64,
        6855257088226130262u64,
        553476580979119549u64,
        6516741293351590684u64,
    ]));
    /// COEFF_B = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFEFFFF5FFD
    #[rustfmt::skip]
    const COEFF_B: Fq = field_new!(Fq, BigInteger256([
        18446355550968045916u64,
        10902955289292811939u64,
        3147092737149958754u64,
        6710871716016002197u64,
    ]));

    type TEModelParameters = EdwardsParameters;
}
