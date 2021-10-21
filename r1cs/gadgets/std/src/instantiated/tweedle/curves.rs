use crate::{
    groups::curves::short_weierstrass::short_weierstrass_jacobian::AffineGadget,
    instantiated::tweedle::{FqGadget, FrGadget},
};
use algebra::{
    curves::tweedle::{dee::TweedledeeParameters, dum::TweedledumParameters},
    fields::tweedle::{Fq, Fr},
};

pub type TweedleDeeGadget = AffineGadget<TweedledeeParameters, Fq, FqGadget>;
pub type TweedleDumGadget = AffineGadget<TweedledumParameters, Fr, FrGadget>;

#[test]
fn test_dee() {
    crate::groups::test::group_test_with_unsafe_add::<_, _, TweedleDeeGadget>();
}

#[test]
fn test_dum() {
    crate::groups::test::group_test_with_unsafe_add::<_, _, TweedleDumGadget>();
}
