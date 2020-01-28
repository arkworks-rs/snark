use crate::{
    fields::tests::{field_test, frobenius_test, primefield_test, sqrt_field_test},
    Field,
};
use rand::{rngs::OsRng, Rng};

#[test]
fn test_mnt6_fr() {
    use crate::fields::mnt6::Fr;

    let a: Fr = OsRng.gen();
    let b: Fr = OsRng.gen();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fr>();
}

#[test]
fn test_mnt6_fq() {
    use crate::fields::mnt6::Fq;

    let a: Fq = OsRng.gen();
    let b: Fq = OsRng.gen();
    field_test(a, b);
    sqrt_field_test(a);
    primefield_test::<Fq>();
}

#[test]
fn test_mnt6_fq3() {
    use crate::fields::mnt6::{Fq, Fq3};

    let a: Fq3 = OsRng.gen();
    let b: Fq3 = OsRng.gen();
    field_test(a, b);
    sqrt_field_test(a);
    frobenius_test::<Fq3, _>(Fq::characteristic(), 13);
}

#[test]
fn test_mnt6_fq6() {
    use crate::fields::mnt6::{Fq, Fq6};

    let a: Fq6 = OsRng.gen();
    let b: Fq6 = OsRng.gen();
    field_test(a, b);
    frobenius_test::<Fq6, _>(Fq::characteristic(), 13);
}
