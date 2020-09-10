#![allow(unused_imports)]
use algebra_core::{curves::models::SWModelParameters, fields::SquareRootField, Zero};
use core::ops::{AddAssign, MulAssign};

use crate::bls12_381::*;
std_curve_tests!(Bls12_381, Fq12);

#[test]
#[cfg(feature = "all_tests")]
fn test_g1_generator_raw() {
    let mut x = Fq::zero();
    let mut i = 0;
    loop {
        // y^2 = x^3 + b
        let mut rhs = x;
        rhs.square_in_place();
        rhs.mul_assign(&x);
        rhs.add_assign(&g1::Parameters::COEFF_B);

        if let Some(y) = rhs.sqrt() {
            let p = G1Affine::new(x, if y < -y { y } else { -y }, false);
            assert!(!p.is_in_correct_subgroup_assuming_on_curve());

            let g1 = p.scale_by_cofactor();
            if !g1.is_zero() {
                assert_eq!(i, 4);
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
