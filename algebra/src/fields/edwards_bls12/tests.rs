use crate::fields::tests::{field_test, primefield_test};

#[test]
fn test_edwards_bls12_fr() {
    use crate::fields::edwards_bls12::fr::Fr;

    let a: Fr = rand::random();
    let b: Fr = rand::random();
    field_test(a, b);
    primefield_test::<Fr>();
}

#[test]
fn test_edwards_bls12_fq() {
    use crate::fields::edwards_bls12::fq::Fq;

    let a: Fq = rand::random();
    let b: Fq = rand::random();
    field_test(a, b);
    primefield_test::<Fq>();
}
