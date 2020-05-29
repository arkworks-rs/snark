use algebra::bls12_377::Parameters;

pub type PairingGadget = crate::pairing::bls12::PairingGadget<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::Bls12_377, _, PairingGadget>()
}
