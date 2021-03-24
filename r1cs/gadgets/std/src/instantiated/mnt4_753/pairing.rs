use algebra::curves::mnt4753::MNT4_753Parameters as Parameters;

pub type PairingGadget = crate::pairing::mnt4::MNT4PairingGadget<Parameters>;

#[test]
fn test() {
    crate::pairing::tests::bilinearity_test::<algebra::curves::mnt4753::MNT4, _, PairingGadget>()
}
