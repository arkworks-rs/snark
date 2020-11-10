#[allow(unused_macros)]
macro_rules! std_curve_tests {
    ($CURVE_IDENT: ident, $GTField: ident) => {
        use algebra_core::{
            test_rng, AffineCurve, Field, One, PairingEngine, PrimeField, ProjectiveCurve,
            UniformRand,
        };
        use rand::Rng;

        use crate::tests::{cuda::*, curves::*, groups::*, msm::*};

        #[test]
        #[cfg(feature = "curve")]
        fn test_g1_curve() {
            curve_tests::<G1Projective>();
        }

        #[test]
        #[cfg(any(
            feature = "serialisation",
            feature = "verify",
            feature = "random_bytes"
        ))]
        fn test_sw_g1() {
            sw_tests::<g1::Parameters>();
        }

        #[test]
        #[cfg(feature = "curve")]
        fn test_g2_curve() {
            curve_tests::<G2Projective>();
        }

        #[test]
        #[cfg(any(
            feature = "serialisation",
            feature = "verify",
            feature = "random_bytes"
        ))]
        fn test_sw_g2() {
            sw_tests::<g2::Parameters>();
        }

        #[test]
        #[cfg(feature = "batch_affine")]
        fn test_batch_affine_g1() {
            batch_affine_test::<G1Projective>();
        }

        #[test]
        #[cfg(feature = "batch_affine")]
        fn test_batch_affine_g2() {
            batch_affine_test::<G2Projective>();
        }

        #[test]
        #[cfg(feature = "curve")]
        fn test_g1_group() {
            let mut rng = test_rng();
            let a: G1Projective = rng.gen();
            let b: G1Projective = rng.gen();
            group_test(a, b);
        }

        #[test]
        #[cfg(feature = "curve")]
        fn test_g1_generator() {
            let generator = G1Affine::prime_subgroup_generator();
            assert!(generator.is_on_curve());
            assert!(generator.is_in_correct_subgroup_assuming_on_curve());
        }

        #[test]
        #[cfg(feature = "curve")]
        fn test_g2_group() {
            let mut rng = test_rng();
            let a: G2Projective = rng.gen();
            let b: G2Projective = rng.gen();
            group_test(a, b);
        }

        #[test]
        #[cfg(feature = "curve")]
        fn test_g2_generator() {
            let generator = G2Affine::prime_subgroup_generator();
            assert!(generator.is_on_curve());
            assert!(generator.is_in_correct_subgroup_assuming_on_curve());
        }

        #[test]
        #[cfg(feature = "msm")]
        fn test_g1_msm() {
            test_msm::<G1Affine>();
        }

        #[test]
        #[cfg(feature = "msm")]
        fn test_g2_msm() {
            test_msm::<G2Affine>();
        }

        #[test]
        #[cfg(any(feature = "curve", feature = "cuda_test"))]
        fn test_g1_cuda_scalar_mul() {
            test_cuda_scalar_mul::<G1Affine>();
        }

        #[test]
        #[cfg(any(feature = "curve", feature = "cuda_test"))]
        fn test_g2_cuda_scalar_mul() {
            test_cuda_scalar_mul::<G2Affine>();
        }

        #[test]
        #[cfg(feature = "pairing")]
        fn test_bilinearity() {
            let mut rng = test_rng();
            let a: G1Projective = rng.gen();
            let b: G2Projective = rng.gen();
            let s: Fr = rng.gen();

            let sa = a.mul(s);
            let sb = b.mul(s);

            let ans1 = $CURVE_IDENT::pairing(sa, b);
            let ans2 = $CURVE_IDENT::pairing(a, sb);
            let ans3 = $CURVE_IDENT::pairing(a, b).pow(s.into_repr());

            assert_eq!(ans1, ans2);
            assert_eq!(ans2, ans3);

            assert_ne!(ans1, $GTField::one());
            assert_ne!(ans2, $GTField::one());
            assert_ne!(ans3, $GTField::one());

            assert_eq!(ans1.pow(Fr::characteristic()), $GTField::one());
            assert_eq!(ans2.pow(Fr::characteristic()), $GTField::one());
            assert_eq!(ans3.pow(Fr::characteristic()), $GTField::one());
        }

        #[test]
        #[cfg(feature = "pairing")]
        fn test_product_of_pairings() {
            let rng = &mut test_rng();

            let a = G1Projective::rand(rng).into_affine();
            let b = G2Projective::rand(rng).into_affine();
            let c = G1Projective::rand(rng).into_affine();
            let d = G2Projective::rand(rng).into_affine();
            let ans1 = $CURVE_IDENT::pairing(a, b) * &$CURVE_IDENT::pairing(c, d);
            let ans2 =
                $CURVE_IDENT::product_of_pairings(&[(a.into(), b.into()), (c.into(), d.into())]);
            assert_eq!(ans1, ans2);
        }
    };
}

#[allow(unused_macros)]
macro_rules! edwards_curve_tests {
    () => {
        use algebra_core::{
            curves::{AffineCurve, ProjectiveCurve},
            test_rng,
        };
        use rand::Rng;

        use crate::tests::{cuda::*, curves::*, groups::*, msm::*};

        #[test]
        #[cfg(feature = "curve")]
        fn test_curve() {
            curve_tests::<EdwardsProjective>();
        }

        #[test]
        #[cfg(any(
            feature = "serialisation",
            feature = "verify",
            feature = "random_bytes"
        ))]
        fn test_edwards() {
            edwards_tests::<EdwardsParameters>();
        }

        #[test]
        #[cfg(feature = "curve")]
        fn test_group() {
            let mut rng = test_rng();
            let a = rng.gen();
            let b = rng.gen();

            for _i in 0..100 {
                group_test::<EdwardsProjective>(a, b);
            }
        }

        #[test]
        #[cfg(feature = "batch_affine")]
        fn test_batch_affine() {
            batch_affine_test::<EdwardsProjective>();
        }

        #[test]
        #[cfg(feature = "curve")]
        fn test_affine_group() {
            let mut rng = test_rng();
            let a: EdwardsAffine = rng.gen();
            let b: EdwardsAffine = rng.gen();
            for _i in 0..100 {
                group_test::<EdwardsAffine>(a, b);
            }
        }

        #[test]
        #[cfg(feature = "msm")]
        fn test_affine_msm() {
            test_msm::<EdwardsAffine>();
        }

        #[test]
        #[cfg(any(feature = "curve", feature = "cuda_test"))]
        fn test_edwards_cuda_scalar_mul() {
            test_cuda_scalar_mul::<EdwardsAffine>();
        }

        #[test]
        #[cfg(any(feature = "curve", feature = "cuda_test"))]
        fn test_generator() {
            let generator = EdwardsAffine::prime_subgroup_generator();
            assert!(generator.is_on_curve());
            assert!(generator.is_in_correct_subgroup_assuming_on_curve());
        }

        #[test]
        #[cfg(feature = "conversion")]
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
        #[cfg(feature = "conversion")]
        fn test_montgomery_conversion() {
            montgomery_conversion_test::<EdwardsParameters>();
        }
    };
}
