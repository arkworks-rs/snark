use algebra_core::test_rng;
use rand::Rng;
use crate::tests::fields::{field_test, primefield_test};

use crate::edwards_bls12::{Fr, Fq};

#[test]
fn test_fr() {
    let mut rng = test_rng();
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();
    field_test(a, b);
    primefield_test::<Fr>();
}

#[test]
fn test_fq() {
    let mut rng = test_rng();
    let a: Fq = rng.gen();
    let b: Fq = rng.gen();
    field_test(a, b);
    primefield_test::<Fq>();
}
