use crate::fields::fp::FpVar;

pub type FqVar = FpVar<algebra::ed_on_bn254::Fq>;

#[test]
fn test() {
    crate::fields::tests::field_test::<_, _, FqVar>().unwrap();
}
