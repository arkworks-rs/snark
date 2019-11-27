use crate::{fields::tests::{field_test, frobenius_test, primefield_test, sqrt_field_test}, Field};

#[test]
fn test_mnt4753_fr() {
    use crate::fields::mnt4753::Fr;

    let a: Fr = rand::random();
    let b: Fr = rand::random();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fr>();
}

#[test]
fn test_mnt4753_fq() {
    use crate::fields::mnt4753::Fq;

    let a: Fq = rand::random();
    let b: Fq = rand::random();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fq>();
}

#[test]
fn test_mnt4753_fq2() {
    use crate::fields::mnt4753::{Fq2, Fq};

    let a: Fq2 = rand::random();
    let b: Fq2 = rand::random();
    field_test(a, b);
    sqrt_field_test(a);
    frobenius_test::<Fq2, _>(Fq::characteristic(), 13);
}

#[test]
fn test_mnt4753_fq4() {
    use crate::fields::mnt4753::Fq4;
    //use crate::fields::mnt4753::Fq;

    let a: Fq4 = rand::random();
    let b: Fq4 = rand::random();
    field_test(a, b);
    //frobenius_test::<Fq4, _>(Fq::characteristic(), 13);
}
