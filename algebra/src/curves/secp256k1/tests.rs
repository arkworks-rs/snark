use crate::curves::secp256k1::{Projective, Affine, Secp256k1Parameters};
use crate::curves::tests::{curve_tests, sw_jacobian_curve_serialization_test};
use crate::curves::AffineCurve;
use rand_xorshift::XorShiftRng;
use rand::{SeedableRng, Rng};
use crate::groups::tests::group_test;

#[test]
fn test_secp256k1_projective_curve() {
    curve_tests::<Projective>();
    sw_jacobian_curve_serialization_test::<Secp256k1Parameters>();
}

#[test]
fn test_secp256k1_projective_group() {
    let mut rng = XorShiftRng::seed_from_u64(1234567890u64);
    let a: Projective = rng.gen();
    let b: Projective = rng.gen();
    group_test(a, b);
}

#[test]
fn test_secp256k1_generator() {
    let generator = Affine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}