use crate::groups::bls12::{
    G1Gadget as Bls12G1Gadget, G1PreparedGadget as Bls12G1PreparedGadget,
    G2Gadget as Bls12G2Gadget, G2PreparedGadget as Bls12G2PreparedGadget,
};
use algebra::curves::bls12_377::Bls12_377Parameters;

pub type G1Gadget = Bls12G1Gadget<Bls12_377Parameters>;
pub type G2Gadget = Bls12G2Gadget<Bls12_377Parameters>;

pub type G1PreparedGadget = Bls12G1PreparedGadget<Bls12_377Parameters>;
pub type G2PreparedGadget = Bls12G2PreparedGadget<Bls12_377Parameters>;

#[test]
fn test() {
    crate::groups::test::group_test_with_unsafe_add::<_, _, G1Gadget>();
    crate::groups::test::group_test_with_unsafe_add::<_, _, G2Gadget>();
}
