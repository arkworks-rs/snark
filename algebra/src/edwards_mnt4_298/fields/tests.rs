use crate::tests::fields::{field_test, primefield_test};
use algebra_core::{test_rng, FpParameters, fields::{PrimeField, SquareRootField, Field}, One};
use rand::Rng;

use crate::edwards_mnt4_298::{Fq, Fr, FrParameters};

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

#[test]
fn test_fr_root_of_unity() {
    assert_eq!(FrParameters::TWO_ADICITY, 1);
    assert_eq!(
        Fr::multiplicative_generator().pow([
            7767783825863817195u64,
            16719789556019334556u64,
            15662913863871949398u64,
            8380289145304910481u64,
            513768135310u64
        ]),
        Fr::root_of_unity()
    );
    assert_eq!(
        Fr::root_of_unity().pow([1 << FrParameters::TWO_ADICITY]),
        Fr::one()
    );
    assert!(Fr::multiplicative_generator().sqrt().is_none());
}
