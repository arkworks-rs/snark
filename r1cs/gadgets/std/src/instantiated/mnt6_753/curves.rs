use crate::groups::mnt::mnt6;
use algebra::curves::mnt6753::MNT6_753Parameters as Parameters;

pub type G1Gadget = mnt6::G1Gadget<Parameters>;
pub type G2Gadget = mnt6::G2Gadget<Parameters>;

pub type G1PreparedGadget = mnt6::G1PreparedGadget<Parameters>;
pub type G2PreparedGadget = mnt6::G2PreparedGadget<Parameters>;

#[test]
fn test() {
    crate::groups::test::group_test_with_unsafe_add::<_, _, G2Gadget>();
    crate::groups::test::group_test_with_unsafe_add::<_, _, G2Gadget>();
}
