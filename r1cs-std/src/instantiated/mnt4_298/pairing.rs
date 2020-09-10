use algebra::mnt4_298::Parameters;

/// Specifies the constraints for computing a pairing in the MNT4-298 bilinear group.
pub type PairingVar = crate::pairing::mnt4::PairingVar<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::MNT4_298, PairingVar>().unwrap()
}
