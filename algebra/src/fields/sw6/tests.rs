use crate::{
    fields::tests::{field_test, frobenius_test, primefield_test, sqrt_field_test},
    Field,
};

#[test]
fn test_sw6_fr() {
    use crate::fields::sw6::Fr;

    let a: Fr = rand::random();
    let b: Fr = rand::random();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fr>();
}

#[test]
fn test_sw6_fq() {
    use crate::fields::sw6::Fq;

    let a: Fq = rand::random();
    let b: Fq = rand::random();
    field_test(a, b);
    primefield_test::<Fq>();
    sqrt_field_test(a);
}

#[test]
fn test_sw6_fq3() {
    use crate::fields::sw6::{Fq, Fq3};

    let a: Fq3 = rand::random();
    let b: Fq3 = rand::random();
    field_test(a, b);
    sqrt_field_test(a);
    frobenius_test::<Fq3, _>(Fq::characteristic(), 13);
}

#[test]
fn test_sw6_fq6() {
    use crate::fields::sw6::{Fq, Fq6};

    let a: Fq6 = rand::random();
    let b: Fq6 = rand::random();
    field_test(a, b);
    frobenius_test::<Fq6, _>(Fq::characteristic(), 13);
}
