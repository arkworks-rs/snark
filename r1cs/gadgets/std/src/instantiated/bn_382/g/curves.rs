use crate::{
    bn_382::g::FqGadget,
    groups::curves::short_weierstrass::AffineGadget
};
use algebra::{
    fields::bn_382::Fr,
    curves::bn_382::g::Bn382GParameters,
};

pub type Bn382GGadget = AffineGadget<Bn382GParameters, Fr, FqGadget>;


#[test]
fn test() {
    crate::groups::test::group_test_with_unsafe_add::<
        _, _, Bn382GGadget
    >();
}
