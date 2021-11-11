use crate::curves::tests::sw_projective_tests;
use crate::{
    curves::{
        mnt6::{
            g1::MNT6G1Parameters, g2::MNT6G2Parameters, G1Affine, G1Projective, G2Affine,
            G2Projective, MNT6,
        },
        tests::curve_tests,
        AffineCurve, PairingEngine,
    },
    fields::mnt6::fr::Fr,
    groups::tests::group_test,
    SemanticallyValid,
};
use rand;

#[test]
fn test_g1_projective_curve() {
    curve_tests::<G1Projective>();
    sw_projective_tests::<MNT6G1Parameters>()
}

#[test]
fn test_g1_projective_group() {
    let a: G1Projective = rand::random();
    let b: G1Projective = rand::random();
    group_test(a, b);
}

#[test]
fn test_g1_generator() {
    let generator = G1Affine::prime_subgroup_generator();
    assert!(generator.is_valid());
}

#[test]
fn test_g2_projective_curve() {
    curve_tests::<G2Projective>();
    sw_projective_tests::<MNT6G2Parameters>()
}

#[test]
fn test_g2_projective_group() {
    let a: G2Projective = rand::random();
    let b: G2Projective = rand::random();
    group_test(a, b);
}

#[test]
fn test_g2_generator() {
    let generator = G2Affine::prime_subgroup_generator();
    assert!(generator.is_valid());
}

#[test]
fn test_bilinearity() {
    use crate::fields::{mnt6::fq6::Fq6, Field, PrimeField};

    let a: G1Projective = rand::random();
    let b: G2Projective = rand::random();
    let s: Fr = rand::random();

    let sa = a * &s;
    let sb = b * &s;

    let ans1 = MNT6::pairing(sa, b);
    let ans2 = MNT6::pairing(a, sb);
    let ans3 = MNT6::pairing(a, b).pow(s.into_repr());

    assert_eq!(ans1, ans2);
    assert_eq!(ans2, ans3);

    assert_ne!(ans1, Fq6::one());
    assert_ne!(ans2, Fq6::one());
    assert_ne!(ans3, Fq6::one());

    assert_eq!(ans1.pow(Fr::characteristic()), Fq6::one());
    assert_eq!(ans2.pow(Fr::characteristic()), Fq6::one());
    assert_eq!(ans3.pow(Fr::characteristic()), Fq6::one());
}

#[test]
fn test_product_of_pairings() {
    use crate::{ProjectiveCurve, UniformRand};
    let rng = &mut rand::thread_rng();

    let a = G1Projective::rand(rng).into_affine();
    let b = G2Projective::rand(rng).into_affine();
    let c = G1Projective::rand(rng).into_affine();
    let d = G2Projective::rand(rng).into_affine();
    let ans1 = MNT6::pairing(a, b) * &MNT6::pairing(c, d);
    let ans2 = MNT6::product_of_pairings(&[(a.into(), b.into()), (c.into(), d.into())]);
    assert_eq!(ans1, ans2);
}
