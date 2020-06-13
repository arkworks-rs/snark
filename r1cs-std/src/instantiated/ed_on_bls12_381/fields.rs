use crate::fields::fp::FpGadget;

pub type FqGadget = FpGadget<algebra::ed_on_bls12_381::Fq>;

#[test]
fn test() {
    crate::fields::tests::field_test::<_, algebra::ed_on_bls12_381::Fq, FqGadget>();
}
