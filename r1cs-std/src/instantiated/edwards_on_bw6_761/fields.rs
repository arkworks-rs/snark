use crate::fields::fp::FpGadget;
use algebra::edwards_on_cp6_782::fq::Fq;

pub type FqGadget = FpGadget<Fq>;

#[test]
fn test() {
    crate::fields::tests::field_test::<_, Fq, FqGadget>();
}
