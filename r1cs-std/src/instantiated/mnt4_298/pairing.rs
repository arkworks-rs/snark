use algebra::mnt4_298::Parameters;

pub type PairingGadget = crate::pairing::mnt4::PairingGadget<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::MNT4_298, _, PairingGadget>()
}
