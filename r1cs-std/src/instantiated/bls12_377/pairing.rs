use algebra::bls12_377::Parameters;

pub type PairingVar = crate::pairing::bls12::PairingVar<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::Bls12_377, PairingVar>().unwrap()
}
