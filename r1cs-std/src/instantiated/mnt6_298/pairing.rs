use algebra::mnt6_298::Parameters;

/// Specifies the constraints for computing a pairing in the MNT6-298 bilinear group.
pub type PairingVar = crate::pairing::mnt6::PairingVar<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::MNT6_298, PairingVar>().unwrap()
}
