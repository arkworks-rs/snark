use crate::fields::fp::FpGadget;
use algebra::ed_on_cp6_782::fq::Fq;

pub type FqGadget = FpGadget<Fq>;

#[test]
fn test() {
    crate::fields::tests::field_test::<_, _, Fq, FqGadget>();
}
