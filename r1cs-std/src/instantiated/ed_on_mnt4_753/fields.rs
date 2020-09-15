use crate::fields::fp::FpVar;
use algebra::ed_on_mnt4_753::fq::Fq;

pub type FqVar = FpVar<Fq>;

#[test]
fn test() {
    crate::fields::tests::field_test::<_, _, FqVar>().unwrap();
}
