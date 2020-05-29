use crate::groups::bls12;
use algebra::bls12_377::Parameters;

pub type G1Gadget = bls12::G1Gadget<Parameters>;
pub type G2Gadget = bls12::G2Gadget<Parameters>;

pub type G1PreparedGadget = bls12::G1PreparedGadget<Parameters>;
pub type G2PreparedGadget = bls12::G2PreparedGadget<Parameters>;

#[test]
fn test() {
    use algebra::curves::models::bls12::Bls12Parameters;
    crate::groups::curves::short_weierstrass::test::<
        _,
        <Parameters as Bls12Parameters>::G1Parameters,
        G1Gadget,
    >();
    crate::groups::curves::short_weierstrass::test::<
        _,
        <Parameters as Bls12Parameters>::G2Parameters,
        G2Gadget,
    >();
}
