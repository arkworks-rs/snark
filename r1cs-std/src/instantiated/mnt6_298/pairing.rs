use algebra::mnt6_298::Parameters;

pub type PairingGadget = crate::pairing::mnt6::PairingGadget<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::MNT6_298, _, PairingGadget>()
}
