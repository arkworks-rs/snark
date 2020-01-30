use crate::{
    fields::tests::{field_test, primefield_test},
    test_rng,
};
use rand::Rng;

#[test]
fn test_edwards_bls12_fr() {
    use crate::fields::edwards_bls12::fr::Fr;

    let mut rng = test_rng();
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();
    field_test(a, b);
    primefield_test::<Fr>();
}

#[test]
fn test_edwards_bls12_fq() {
    use crate::fields::edwards_bls12::fq::Fq;

    let mut rng = test_rng();
    let a: Fq = rng.gen();
    let b: Fq = rng.gen();
    field_test(a, b);
    primefield_test::<Fq>();
}
