use crate::fields::fp::FpGadget;
use algebra::ed_on_mnt4_298::fq::Fq;

pub type FqGadget = FpGadget<Fq>;

#[test]
fn test() {
    crate::fields::tests::field_test::<_, Fq, FqGadget>();
}
