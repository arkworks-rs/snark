use algebra::mnt4_753::Parameters;

pub type PairingVar = crate::pairing::mnt4::PairingVar<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::MNT4_753, PairingVar>().unwrap()
}
