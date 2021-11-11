use crate::groups::bn;
use algebra::curves::bn_382::Bn382Parameters;

pub type G1Gadget = bn::G1Gadget<Bn382Parameters>;
pub type G2Gadget = bn::G2Gadget<Bn382Parameters>;

pub type G1PreparedGadget = bn::G1PreparedGadget<Bn382Parameters>;
pub type G2PreparedGadget = bn::G2PreparedGadget<Bn382Parameters>;

#[test]
fn test() {
    crate::groups::test::group_test_with_unsafe_add::<_, _, G1Gadget>();
    crate::groups::test::group_test_with_unsafe_add::<_, _, G2Gadget>();
}
