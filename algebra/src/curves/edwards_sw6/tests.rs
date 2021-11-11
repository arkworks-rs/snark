use crate::curves::tests::edwards_tests;
use crate::{
    curves::{
        edwards_sw6::*, models::twisted_edwards_extended::tests::montgomery_conversion_test,
        tests::curve_tests, AffineCurve, ProjectiveCurve,
    },
    groups::tests::group_test,
    SemanticallyValid,
};
use rand;

#[test]
fn test_projective_curve() {
    curve_tests::<EdwardsProjective>();
    edwards_tests::<EdwardsParameters>()
}

#[test]
fn test_projective_group() {
    let a = rand::random();
    let b = rand::random();
    for _i in 0..100 {
        group_test::<EdwardsProjective>(a, b);
    }
}

#[test]
fn test_affine_group() {
    let a: EdwardsAffine = rand::random();
    let b: EdwardsAffine = rand::random();
    for _i in 0..100 {
        group_test::<EdwardsAffine>(a, b);
    }
}

#[test]
fn test_generator() {
    let generator = EdwardsAffine::prime_subgroup_generator();
    assert!(generator.is_valid());
}

#[test]
fn test_conversion() {
    let a: EdwardsAffine = rand::random();
    let b: EdwardsAffine = rand::random();
    let a_b = {
        use crate::groups::Group;
        (a + &b).double().double()
    };
    let a_b2 = (a.into_projective() + &b.into_projective())
        .double()
        .double();
    assert_eq!(a_b, a_b2.into_affine());
    assert_eq!(a_b.into_projective(), a_b2);
}

#[test]
fn test_montgomery_conversion() {
    montgomery_conversion_test::<EdwardsParameters>();
}
