use algebra::mnt6_753::Parameters;

pub type PairingGadget = crate::pairing::mnt6::PairingGadget<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::MNT6_753, _, PairingGadget>()
}
