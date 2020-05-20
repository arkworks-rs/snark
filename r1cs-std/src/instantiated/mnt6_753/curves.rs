use crate::groups::mnt6;
use algebra::mnt6_753::Parameters;

pub type G1Gadget = mnt6::G1Gadget<Parameters>;
pub type G2Gadget = mnt6::G2Gadget<Parameters>;

pub type G1PreparedGadget = mnt6::G1PreparedGadget<Parameters>;
pub type G2PreparedGadget = mnt6::G2PreparedGadget<Parameters>;

#[test]
fn test() {
    use algebra::curves::models::mnt6::MNT6Parameters;
    crate::groups::curves::short_weierstrass::test::<
        _,
        <Parameters as MNT6Parameters>::G1Parameters,
        G1Gadget,
    >();
    crate::groups::curves::short_weierstrass::test::<
        _,
        <Parameters as MNT6Parameters>::G2Parameters,
        G2Gadget,
    >();
}
