use crate::groups::mnt4;
use algebra::mnt4_753::Parameters;

pub type G1Gadget = mnt4::G1Gadget<Parameters>;
pub type G2Gadget = mnt4::G2Gadget<Parameters>;

pub type G1PreparedGadget = mnt4::G1PreparedGadget<Parameters>;
pub type G2PreparedGadget = mnt4::G2PreparedGadget<Parameters>;

#[test]
fn test() {
    use algebra::curves::models::mnt4::MNT4Parameters;
    crate::groups::curves::short_weierstrass::test::<
        _,
        <Parameters as MNT4Parameters>::G1Parameters,
        G1Gadget,
    >();
    crate::groups::curves::short_weierstrass::test::<
        _,
        <Parameters as MNT4Parameters>::G2Parameters,
        G2Gadget,
    >();
}
