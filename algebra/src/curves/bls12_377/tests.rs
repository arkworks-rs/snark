#![allow(unused_imports)]
use crate::{
    curves::{
        bls12_377::{
            g1::Bls12_377G1Parameters, Bls12_377, G1Affine, G1Projective, G2Affine, G2Projective,
        },
        models::SWModelParameters,
        tests::curve_tests,
        AffineCurve, PairingEngine, ProjectiveCurve,
    },
    fields::{
        bls12_377::{Fq, Fq12, Fq2, Fr},
        Field, FpParameters, PrimeField, SquareRootField,
    },
    groups::tests::group_test,
};
use std::ops::{AddAssign, MulAssign};

#[test]
fn test_g1_projective_curve() {
    curve_tests::<G1Projective>();
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
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_g2_projective_curve() {
    curve_tests::<G2Projective>();
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
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

//    #[test]
//    fn test_bilinearity() {
//        let a: G1Projective = rand::random();
//        let b: G2Projective = rand::random();
//        let s: Fr = rand::random();
//
//        let sa = a * &s;
//        let sb = b * &s;
//
//        let ans1 = Bls12_377::pairing(sa, b);
//        let ans2 = Bls12_377::pairing(a, sb);
//        let ans3 = Bls12_377::pairing(a, b).pow(s.into_repr());
//
//        assert_eq!(ans1, ans2);
//        assert_eq!(ans2, ans3);
//
//        assert_ne!(ans1, Fq12::one());
//        assert_ne!(ans2, Fq12::one());
//        assert_ne!(ans3, Fq12::one());
//
//        assert_eq!(ans1.pow(Fr::characteristic()), Fq12::one());
//        assert_eq!(ans2.pow(Fr::characteristic()), Fq12::one());
//        assert_eq!(ans3.pow(Fr::characteristic()), Fq12::one());
//    }

#[test]
fn test_bilinearity() {
    let a: G1Projective = G1Projective::prime_subgroup_generator();
    let b: G2Projective = G2Projective::prime_subgroup_generator();
    let s: Fr = Fr::one() + &Fr::one();

    let sa = a * &s;
    let sb = b * &s;

    println!("a\n{:?}\n", a.into_affine());
    println!("b\n{:?}\n", b.into_affine());
    println!("s\n{:?}\n", s);
    println!("sa\n{:?}\n", sa.into_affine());
    println!("sb\n{:?}\n", sb.into_affine());

    let ans1 = Bls12_377::pairing(sa, b);
    let ans2 = Bls12_377::pairing(a, sb);

    assert_eq!(ans1, ans2);

    assert_ne!(ans1, Fq12::one());
    assert_ne!(ans2, Fq12::one());
    assert_eq!(ans1.pow(Fr::characteristic()), Fq12::one());
    assert_eq!(ans2.pow(Fr::characteristic()), Fq12::one());
}

#[test]
fn test_g1_generator_raw() {
    let mut x = Fq::zero();
    let mut i = 0;
    loop {
        // y^2 = x^3 + b
        let mut rhs = x;
        rhs.square_in_place();
        rhs.mul_assign(&x);
        rhs.add_assign(&Bls12_377G1Parameters::COEFF_B);

        if let Some(y) = rhs.sqrt() {
            let p = G1Affine::new(x, if y < -y { y } else { -y }, false);
            assert!(!p.is_in_correct_subgroup_assuming_on_curve());

            let g1 = p.scale_by_cofactor();
            if !g1.is_zero() {
                assert_eq!(i, 1);
                let g1 = G1Affine::from(g1);

                assert!(g1.is_in_correct_subgroup_assuming_on_curve());

                assert_eq!(g1, G1Affine::prime_subgroup_generator());
                break;
            }
        }

        i += 1;
        x.add_assign(&Fq::one());
    }
}

#[test]
fn bls12_377_unique() {
    use crate::fields::bls12_377::fq::Fq;

    use std::str::FromStr;
    println!("{}", Fq::from_str("155198655607781456406391640216936120121836107652948796323930557600032281009004493664981332883744016074664192874906").unwrap());
}
