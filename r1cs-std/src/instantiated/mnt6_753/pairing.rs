use algebra::mnt6_753::Parameters;

pub type PairingVar = crate::pairing::mnt6::PairingVar<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::MNT6_753, PairingVar>().unwrap()
}
