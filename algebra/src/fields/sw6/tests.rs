use crate::{
    fields::tests::{field_test, frobenius_test, primefield_test, sqrt_field_test},
    test_rng, Field,
};
use rand::Rng;

#[test]
fn test_sw6_fr() {
    use crate::fields::sw6::Fr;

    let mut rng = test_rng();
    let a: Fr = rng.gen();
    let b: Fr = rng.gen();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fr>();
}

#[test]
fn test_sw6_fq() {
    use crate::fields::sw6::Fq;

    let mut rng = test_rng();
    let a: Fq = rng.gen();
    let b: Fq = rng.gen();
    field_test(a, b);
    primefield_test::<Fq>();
    sqrt_field_test(a);
}

#[test]
fn test_sw6_fq3() {
    use crate::fields::sw6::{Fq, Fq3};

    let mut rng = test_rng();
    let a: Fq3 = rng.gen();
    let b: Fq3 = rng.gen();
    field_test(a, b);
    sqrt_field_test(a);
    frobenius_test::<Fq3, _>(Fq::characteristic(), 13);
}

#[test]
fn test_sw6_fq6() {
    use crate::fields::sw6::{Fq, Fq6};

    let mut rng = test_rng();
    let a: Fq6 = rng.gen();
    let b: Fq6 = rng.gen();
    field_test(a, b);
    frobenius_test::<Fq6, _>(Fq::characteristic(), 13);
}
