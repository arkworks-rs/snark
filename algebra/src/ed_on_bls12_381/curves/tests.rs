use algebra_core::{
    bytes::FromBytes,
    curves::{AffineCurve, ProjectiveCurve},
    test_rng, Zero,
};
use core::str::FromStr;
use rand::Rng;

use crate::ed_on_bls12_381::*;

use crate::tests::{curves::*, groups::*};

#[test]
fn test_projective_curve() {
    curve_tests::<EdwardsProjective>();

    edwards_tests::<EdwardsParameters>();
}

#[test]
fn test_projective_group() {
    let mut rng = test_rng();
    let a = rng.gen();
    let b = rng.gen();
    for _i in 0..100 {
        group_test::<EdwardsProjective>(a, b);
    }
}

#[test]
fn test_affine_group() {
    let mut rng = test_rng();
    let a: EdwardsAffine = rng.gen();
    let b: EdwardsAffine = rng.gen();
    for _i in 0..100 {
        group_test::<EdwardsAffine>(a, b);
    }
}

#[test]
fn test_generator() {
    let generator = EdwardsAffine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_conversion() {
    let mut rng = test_rng();
    let a: EdwardsAffine = rng.gen();
    let b: EdwardsAffine = rng.gen();
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
fn test_scalar_multiplication() {
    println!("Started getting field elements");
    let f1 = Fr::from_str(
        "4691331900926794624732159288782398864809513177368446695323460897088210774597",
    )
    .unwrap();
    let f2 = Fr::from_str(
        "1305028103380024953477151132159456965337646722479526711736847301646466538045",
    )
    .unwrap();

    println!("Finished getting field elements");
    let g = EdwardsAffine::from_str(
        "(1158870117176967269192899343636553522971009777237254192973081388797299308391, \
         36933624999642413792569726058244472742169727126562409632889593958355839948294)",
    )
    .unwrap();
    let f1f2g = EdwardsAffine::from_str(
        "(12638652891150111215300246576936483137884466359309882317048163368620501191944, \
         38385045634663742820428406709832518145724237919360177362175527604556651918148)",
    )
    .unwrap();

    println!("Finished getting group elements");

    assert!(!g.is_zero());
    assert!(!f1f2g.is_zero());

    let f1g = g.mul(f1).into_affine();
    println!("f1: {:?}", f1);
    println!("f2: {:?}", f2);
    println!("g: {:?}", g);
    println!("f1f2g: {:?}", f1f2g);
    assert_eq!(g.mul(f1 * &f2).into_affine(), f1f2g);
    assert_eq!(f1g.mul(f2).into_affine(), f1f2g);
}

#[test]
fn test_bytes() {
    let g_from_repr = EdwardsAffine::from_str(
        "(1158870117176967269192899343636553522971009777237254192973081388797299308391, \
         36933624999642413792569726058244472742169727126562409632889593958355839948294)",
    )
    .unwrap();

    let g_bytes = algebra_core::to_bytes![g_from_repr].unwrap();
    let g = EdwardsAffine::read(g_bytes.as_slice()).unwrap();
    assert_eq!(g_from_repr, g);
}

#[test]
fn test_montgomery_conversion() {
    montgomery_conversion_test::<EdwardsParameters>();
}
