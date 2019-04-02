use crate::{
    biginteger::BigInteger256,
    curves::{
        models::{ModelParameters, TEModelParameters},
        twisted_edwards_extended::{GroupAffine, GroupProjective},
    },
    fields::jubjub::{fq::Fq, fr::Fr},
};
use std::str::FromStr;

#[cfg(test)]
mod tests;

pub type JubJubAffine = GroupAffine<JubJubParameters>;
pub type JubJubProjective = GroupProjective<JubJubParameters>;

const GENERATOR_X: Fq = Fq::new(BigInteger256([
    14080349899812819339,
    4104857150246327429,
    8293216003873356624,
    7400363483732984990,
]));
const GENERATOR_Y: Fq = Fq::new(BigInteger256([
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
#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct JubJubParameters;

impl ModelParameters for JubJubParameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl TEModelParameters for JubJubParameters {
    /// COEFF_A = -1
    const COEFF_A: Fq = Fq::new(BigInteger256([
        18446744060824649731,
        18102478225614246908,
        11073656695919314959,
        6613806504683796440,
    ]));

    /// COEFF_D = (10240/10241) mod q
    const COEFF_D: Fq = Fq::new(BigInteger256([
        3049539848285517488,
        18189135023605205683,
        8793554888777148625,
        6339087681201251886,
    ]));

    /// COFACTOR = 8
    const COFACTOR: &'static [u64] = &[8];

    /// COFACTOR^(-1) mod r =
    /// 819310549611346726241370945440405716213240158234039660170669895299022906775
    const COFACTOR_INV: Fr = Fr::new(BigInteger256([
        6832491983681988242,
        12911748493335322362,
        17523939349049608702,
        217463794347581613,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (GENERATOR_X, GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) = (GENERATOR_X, GENERATOR_Y);

    /// Multiplication by `a` is simply negation here.
    #[inline(always)]
    fn mul_by_a(elem: &Self::BaseField) -> Self::BaseField {
        -(*elem)
    }
}

impl FromStr for JubJubAffine {
    type Err = ();

    fn from_str(mut s: &str) -> Result<Self, Self::Err> {
        s = s.trim();
        if s.is_empty() {
            return Err(());
        }
        if s.len() < 3 {
            return Err(());
        }
        if !(s.starts_with('(') && s.ends_with(')')) {
            return Err(());
        }
        let mut point = Vec::new();
        for substr in s.split(|c| c == '(' || c == ')' || c == ',' || c == ' ') {
            if !substr.is_empty() {
                point.push(Fq::from_str(substr)?);
            }
        }
        if point.len() != 2 {
            return Err(());
        }
        let point = JubJubAffine::new(point[0], point[1]);

        if !point.is_on_curve() {
            Err(())
        } else {
            Ok(point)
        }
    }
}
