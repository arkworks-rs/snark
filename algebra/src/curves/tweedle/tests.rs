use crate::{
    curves::{
        models::SWModelParameters, tests::curve_tests, tweedle::*, AffineCurve, ProjectiveCurve,
    },
    fields::{tweedle::*, Field, SquareRootField},
    groups::tests::group_test,
};
use std::ops::{AddAssign, MulAssign};
use std::str::FromStr;

use crate::curves::tests::sw_jacobian_tests;
use crate::curves::tweedle::dee::TweedledeeParameters;
use crate::curves::tweedle::dum::TweedledumParameters;
use rand::{Rng, SeedableRng};
use rand_xorshift::XorShiftRng;

#[test]
fn test_dee_projective_curve() {
    curve_tests::<dee::Projective>();
    sw_jacobian_tests::<TweedledeeParameters>()
}

#[test]
fn test_dee_projective_group() {
    let mut rng = XorShiftRng::seed_from_u64(1234567890u64);
    let a: dee::Projective = rng.gen();
    let b: dee::Projective = rng.gen();
    group_test(a, b);
}

#[test]
fn test_dee_generator() {
    let generator = dee::Affine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_dum_projective_curve() {
    curve_tests::<dum::Projective>();
    sw_jacobian_tests::<TweedledumParameters>()
}

#[test]
fn test_dum_projective_group() {
    let mut rng = XorShiftRng::seed_from_u64(1234567890u64);
    let a: dum::Projective = rng.gen();
    let b: dum::Projective = rng.gen();
    group_test(a, b);
}

#[test]
fn test_dum_generator() {
    let generator = dum::Affine::prime_subgroup_generator();
    assert!(generator.is_on_curve());
    assert!(generator.is_in_correct_subgroup_assuming_on_curve());
}

#[test]
fn test_dee_generator_raw() {
    let mut x = Fq::zero();
    let mut i = 0;
    loop {
        // y^2 = x^3 + b
        let mut rhs = x;
        rhs.square_in_place();
        rhs.mul_assign(&x);
        rhs.add_assign(&dee::TweedledeeParameters::COEFF_B);

        if let Some(y) = rhs.sqrt() {
            let p = dee::Affine::new(x, if y < -y { y } else { -y }, false);
            assert!(p.is_in_correct_subgroup_assuming_on_curve());

            let dee = p.scale_by_cofactor();
            assert_eq!(dee.into_affine(), p);

            if !dee.is_zero() {
                assert_eq!(i, 1);
                let dee = dee::Affine::from(dee);

                assert!(dee.is_in_correct_subgroup_assuming_on_curve());

                assert_eq!(dee, dee::Affine::prime_subgroup_generator());
                break;
            }
        }

        i += 1;
        x.add_assign(&Fq::one());
    }
}

#[test]
fn test_dum_generator_raw() {
    let mut x = Fr::zero();
    let mut i = 0;
    loop {
        // y^2 = x^3 + b
        let mut rhs = x;
        rhs.square_in_place();
        rhs.mul_assign(&x);
        rhs.add_assign(&dum::TweedledumParameters::COEFF_B);

        if let Some(y) = rhs.sqrt() {
            let p = dum::Affine::new(x, if y < -y { y } else { -y }, false);
            assert!(p.is_in_correct_subgroup_assuming_on_curve());

            let dum = p.scale_by_cofactor();
            assert_eq!(dum.into_affine(), p);

            if !dum.is_zero() {
                assert_eq!(i, 1);
                let dum = dum::Affine::from(dum);

                assert!(dum.is_in_correct_subgroup_assuming_on_curve());

                assert_eq!(dum, dum::Affine::prime_subgroup_generator());
                break;
            }
        }

        i += 1;
        x.add_assign(&Fr::one());
    }
}

#[test]
fn test_dee_addition_correctness() {
    let mut p = dee::Projective::new(
        Fq::from_str(
            "17071515411234329267051251142008744532074161438140426170549136904789606209155",
        )
        .unwrap(),
        Fq::from_str(
            "9067370984564524093871625068725679070040168060994636121507153477916099620826",
        )
        .unwrap(),
        Fq::one(),
    );

    p.add_assign(&dee::Projective::new(
        Fq::from_str(
            "5902988235118225415057554152593109689819081116067139376852243422243422684655",
        )
        .unwrap(),
        Fq::from_str(
            "9237374262095944048575165674046716194558759078123659312337709713005853948132",
        )
        .unwrap(),
        Fq::one(),
    ));

    let p = dee::Affine::from(p);

    assert_eq!(
        p,
        dee::Affine::new(
            Fq::from_str(
                "17272972729543522859996365140537720509583378385403153153034405894416507370075"
            )
            .unwrap(),
            Fq::from_str(
                "10919319153241406943315020022865635527830995765162202572118118072098170575117"
            )
            .unwrap(),
            false,
        )
    );
}

#[test]
fn test_dum_addition_correctness() {
    let mut p = dum::Projective::new(
        Fr::from_str(
            "21118483776076764996122757821606091900059043860162004907989579660882026321197",
        )
        .unwrap(),
        Fr::from_str(
            "9025588652913915603174720117986570170395425582417356177673155554443430464689",
        )
        .unwrap(),
        Fr::one(),
    );

    p.add_assign(&dum::Projective::new(
        Fr::from_str(
            "20385173229981432379197513268506886433340219379830521001646291041798263137109",
        )
        .unwrap(),
        Fr::from_str(
            "16494790468966191765270742698088328193228887152919586743292725150337386016283",
        )
        .unwrap(),
        Fr::one(),
    ));

    let p = dum::Affine::from(p);

    assert_eq!(
        p,
        dum::Affine::new(
            Fr::from_str(
                "3707088439511374954709258634608802460084680838305626554041952787711711292620"
            )
            .unwrap(),
            Fr::from_str(
                "21427612888550306000000889405343941940930914059283626531936541886782117113518"
            )
            .unwrap(),
            false,
        )
    );
}
