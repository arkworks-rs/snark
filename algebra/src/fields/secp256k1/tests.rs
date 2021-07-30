use crate::fields::tests::{field_test, primefield_test, sqrt_field_test};

#[test]
fn test_secp256k1_fr() {
    use crate::fields::secp256k1::Fr;

    let a: Fr = rand::random();
    let b: Fr = rand::random();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fr>();
}

#[test]
fn test_secp256k1_fq() {
    use crate::fields::secp256k1::Fq;

    let a: Fq = rand::random();
    let b: Fq = rand::random();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fq>();
}