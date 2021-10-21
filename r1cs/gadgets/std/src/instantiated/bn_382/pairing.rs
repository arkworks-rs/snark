use algebra::curves::bn_382::Bn382Parameters;

pub type PairingGadget = crate::pairing::bn::PairingGadget<Bn382Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::curves::bn_382::Bn382, _, PairingGadget>()
}
