use algebra_core::fields::models::fp6_2over3::*;
use algebra_core::fields::quadratic_extension::QuadExtParameters;
use algebra_core::{test_rng, Field};
use rand::Rng;

use crate::mnt6_298::*;

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
fn test_fq3() {
    let mut rng = test_rng();
    let a: Fq3 = rng.gen();
    let b: Fq3 = rng.gen();
    field_test(a, b);
    sqrt_field_test(a);
    frobenius_test::<Fq3, _>(Fq::characteristic(), 13);
    assert_eq!(
        a * Fq6Parameters::NONRESIDUE,
        <Fp6ParamsWrapper<Fq6Parameters>>::mul_base_field_by_nonresidue(&a)
    );
}

#[test]
fn test_fq6() {
    let mut rng = test_rng();
    let a: Fq6 = rng.gen();
    let b: Fq6 = rng.gen();
    field_test(a, b);
    frobenius_test::<Fq6, _>(Fq::characteristic(), 13);
}
