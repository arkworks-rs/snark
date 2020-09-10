use algebra::bls12_377::Parameters;

/// Specifies the constraints for computing a pairing in the BLS12-377 bilinear group.
pub type PairingVar = crate::pairing::bls12::PairingVar<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::Bls12_377, PairingVar>().unwrap()
}
