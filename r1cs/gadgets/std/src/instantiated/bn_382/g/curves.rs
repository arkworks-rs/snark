use crate::{bn_382::g::FqGadget, groups::curves::short_weierstrass::AffineGadget};
use algebra::{curves::bn_382::g::Bn382GParameters, fields::bn_382::Fr};

pub type Bn382GGadget = AffineGadget<Bn382GParameters, Fr, FqGadget>;

#[test]
fn test() {
    crate::groups::test::group_test_with_unsafe_add::<_, _, Bn382GGadget>();
}
