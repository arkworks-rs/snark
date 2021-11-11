use crate::groups::mnt::mnt4;
use algebra::curves::mnt4753::MNT4_753Parameters as Parameters;

pub type G1Gadget = mnt4::G1Gadget<Parameters>;
pub type G2Gadget = mnt4::G2Gadget<Parameters>;

pub type G1PreparedGadget = mnt4::G1PreparedGadget<Parameters>;
pub type G2PreparedGadget = mnt4::G2PreparedGadget<Parameters>;

#[test]
fn test() {
    crate::groups::test::group_test_with_unsafe_add::<_, _, G1Gadget>();
    crate::groups::test::group_test_with_unsafe_add::<_, _, G2Gadget>();
}
