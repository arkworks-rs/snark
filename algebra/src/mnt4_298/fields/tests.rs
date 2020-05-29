use algebra_core::{test_rng, Field};
use rand::Rng;

use crate::mnt4_298::*;

use crate::tests::fields::{field_test, frobenius_test, primefield_test, sqrt_field_test};

#[test]
fn test_fr() {
    let mut rng = test_rng();
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fr>();
}

#[test]
fn test_fq() {
    let mut rng = test_rng();
    let a: Fq = rng.gen();
    let b: Fq = rng.gen();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fq>();
}

#[test]
fn test_fq2() {
    let mut rng = test_rng();
    let a: Fq2 = rng.gen();
    let b: Fq2 = rng.gen();
    field_test(a, b);
    sqrt_field_test(a);
    frobenius_test::<Fq2, _>(Fq::characteristic(), 13);
}

#[test]
fn test_fq4() {
    let mut rng = test_rng();
    let a: Fq4 = rng.gen();
    let b: Fq4 = rng.gen();
    field_test(a, b);
    frobenius_test::<Fq4, _>(Fq::characteristic(), 13);
}
