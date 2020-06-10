use crate::fields::fp::FpGadget;

pub type FqGadget = FpGadget<algebra::edwards_on_bls12_381::Fq>;

#[test]
fn test() {
    crate::fields::tests::field_test::<_, algebra::edwards_on_bls12_381::Fq, FqGadget>();
}
