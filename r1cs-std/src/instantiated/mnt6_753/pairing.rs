use algebra::mnt6_753::Parameters;

/// Specifies the constraints for computing a pairing in the MNT6-753 bilinear group.
pub type PairingVar = crate::pairing::mnt6::PairingVar<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::MNT6_753, PairingVar>().unwrap()
}
