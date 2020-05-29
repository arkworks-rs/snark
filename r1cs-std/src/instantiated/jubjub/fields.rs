use crate::fields::fp::FpGadget;

pub type FqGadget = FpGadget<algebra::jubjub::Fq>;

#[test]
fn test() {
    crate::fields::tests::field_test::<_, algebra::jubjub::Fq, FqGadget>();
}
