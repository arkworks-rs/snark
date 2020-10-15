use crate::groups::bn;
use algebra::curves::bn_382::g::Bn382GParameters;

pub type G1Gadget = bn::G1Gadget<Bn382GParameters>;
pub type G2Gadget = bn::G2Gadget<Bn382GParameters>;

pub type G1PreparedGadget = bn::G1PreparedGadget<Bn382GParameters>;
pub type G2PreparedGadget = bn::G2PreparedGadget<Bn382GParameters>;

#[test]
fn test() {
    crate::groups::test::group_test_with_unsafe_add::<
        _, _, G1Gadget
    >();
    crate::groups::test::group_test_with_unsafe_add::<
        _, _, G2Gadget
    >();
}
