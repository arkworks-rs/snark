use crate::{
    biginteger::BigInteger384,
    curves::{
        bn_382::*,
        models::SWModelParameters,
        tests::{curve_tests, sw_jacobian_tests},
        AffineCurve, PairingEngine, ProjectiveCurve,
    },
    field_new,
    fields::{bn_382::*, Field, SquareRootField},
    groups::tests::group_test,
};
use std::ops::{AddAssign, MulAssign};

use rand::{Rng, SeedableRng};
use rand_xorshift::XorShiftRng;

#[test]
fn test_g1_projective_curve() {
    curve_tests::<G1Projective>();
    sw_jacobian_tests::<g1::Bn382G1Parameters>()
}

#[test]
fn test_g1_projective_group() {
    let mut rng = XorShiftRng::seed_from_u64(1234567890u64);
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
    sw_jacobian_tests::<g2::Bn382G2Parameters>()
}

#[test]
fn test_g2_projective_group() {
    let mut rng = XorShiftRng::seed_from_u64(1234567890u64);
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
    let a: G1Projective = G1Projective::prime_subgroup_generator();
    let b: G2Projective = G2Projective::prime_subgroup_generator();
    let s: Fr = Fr::one() + &Fr::one();

    let sa = a * &s;
    let sb = b * &s;

    let ans1 = Bn382::pairing(sa, b).unwrap();
    let ans2 = Bn382::pairing(a, sb).unwrap();

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
        rhs.add_assign(&g1::Bn382G1Parameters::COEFF_B);

        if let Some(y) = rhs.sqrt() {
            let p = G1Affine::new(x, if y < -y { y } else { -y }, false);
            assert!(p.is_in_correct_subgroup_assuming_on_curve());

            let g1 = p.scale_by_cofactor();
            assert_eq!(g1.into_affine(), p);

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
fn test_g1_addition_correctness() {
    let mut p = G1Projective::new(
        field_new!(
            Fq,
            BigInteger384([
                0xb93d80f690db69d0,
                0x4f067992e8332718,
                0x34158a73ed82a5b5,
                0x76579a71fa073da3,
                0xd18af844d8e19090,
                0x132f6baaf304be2d
            ])
        ),
        field_new!(
            Fq,
            BigInteger384([
                0x94fa65a8bcfb4667,
                0x40c1a252a0ebcd03,
                0xb015731e87b1e56d,
                0x48ddb455fce53c8b,
                0x927dd51c6d8710d9,
                0x180205df8eb25b6f
            ])
        ),
        Fq::one(),
    );

    p.add_assign(&G1Projective::new(
        field_new!(
            Fq,
            BigInteger384([
                0x83e3f6c6d507009e,
                0x117f1d6b6a8e8015,
                0x8b35968a3723188c,
                0xb0fb67dcdd7acad3,
                0x7b7c13f80f311614,
                0xa98bd8574ae9bc0
            ])
        ),
        field_new!(
            Fq,
            BigInteger384([
                0x16da0a94ca80a05d,
                0x21b711d1e0aab153,
                0xbc6aad515f9d2f13,
                0x61338d9dd736bbbf,
                0x46993eca6ea17f51,
                0x1905aa15a2782887
            ])
        ),
        Fq::one(),
    ));

    let p = G1Affine::from(p);

    assert_eq!(
        p,
        G1Affine::new(
            field_new!(
                Fq,
                BigInteger384([
                    0x66938afc964e6bb1,
                    0x404f672ef202544c,
                    0xd134de56fd2929da,
                    0xb36a988806affd22,
                    0x778395d12bfa7bb3,
                    0x1c7c8ce70cc8e6ba
                ])
            ),
            field_new!(
                Fq,
                BigInteger384([
                    0x2ac6ee03d57352e1,
                    0xa1797851137cedaf,
                    0x77c6655d51aa0e5f,
                    0x7b110de6353a0c7a,
                    0xf2838e6295889402,
                    0x1322f26a63efe9e4
                ])
            ),
            false,
        )
    );
}
