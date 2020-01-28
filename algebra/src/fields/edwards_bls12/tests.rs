use crate::fields::tests::{field_test, primefield_test};
use rand::{rngs::OsRng, Rng};

#[test]
fn test_edwards_bls12_fr() {
    use crate::fields::edwards_bls12::fr::Fr;

    let a: Fr = OsRng.gen();
    let b: Fr = OsRng.gen();
    field_test(a, b);
    primefield_test::<Fr>();
}

#[test]
fn test_edwards_bls12_fq() {
    use crate::fields::edwards_bls12::fq::Fq;

    let a: Fq = OsRng.gen();
    let b: Fq = OsRng.gen();
    field_test(a, b);
    primefield_test::<Fq>();
}
