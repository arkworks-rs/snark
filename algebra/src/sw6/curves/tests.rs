use algebra_core::{
    test_rng, AffineCurve, CanonicalSerialize, Field, One, PairingEngine, PrimeField,
    ProjectiveCurve,
};
use rand::Rng;

use crate::sw6::*;

use crate::tests::{curves::*, groups::*};

#[test]
fn test_g1_projective_curve() {
    curve_tests::<G1Projective>();

    let byte_size = <G1Affine as CanonicalSerialize>::buffer_size();
    sw_curve_serialization_test::<g1::Parameters>(byte_size);
}

#[test]
fn test_g1_projective_group() {
    let mut rng = test_rng();
    let a: G1Projective = rng.gen();
    let b: G1Projective = rng.gen();
    group_test(a, b);
}

#[test]
fn test_g1_generator() {
    let generator = G1Affine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_g2_projective_curve() {
    curve_tests::<G2Projective>();

    let byte_size = <G2Affine as CanonicalSerialize>::buffer_size();
    sw_curve_serialization_test::<g2::Parameters>(byte_size);
}

#[test]
fn test_g2_projective_group() {
    let mut rng = test_rng();
    let a: G2Projective = rng.gen();
    let b: G2Projective = rng.gen();
    group_test(a, b);
}

#[test]
fn test_g2_generator() {
    let generator = G2Affine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_bilinearity() {
    let mut rng = test_rng();
    let a: G1Projective = rng.gen();
    let b: G2Projective = rng.gen();
    let s: Fr = rng.gen();

    let sa = a.mul(s);
    let sb = b.mul(s);

    let ans1 = SW6::pairing(sa, b);
    let ans2 = SW6::pairing(a, sb);
    let ans3 = SW6::pairing(a, b).pow(s.into_repr());

    assert_eq!(ans1, ans2);
    assert_eq!(ans2, ans3);

    assert_ne!(ans1, Fq6::one());
    assert_ne!(ans2, Fq6::one());
    assert_ne!(ans3, Fq6::one());

    assert_eq!(ans1.pow(Fr::characteristic()), Fq6::one());
    assert_eq!(ans2.pow(Fr::characteristic()), Fq6::one());
    assert_eq!(ans3.pow(Fr::characteristic()), Fq6::one());
}
