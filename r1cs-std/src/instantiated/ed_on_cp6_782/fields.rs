use crate::fields::fp::FpVar;
use algebra::ed_on_cp6_782::fq::Fq;

/// A variable that is the R1CS equivalent of `algebra::ed_on_cp6_782::Fq`.
pub type FqVar = FpVar<Fq>;

#[test]
fn test() {
    crate::fields::tests::field_test::<_, _, FqVar>().unwrap();
}
