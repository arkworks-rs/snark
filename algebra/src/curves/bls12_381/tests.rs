use crate::{
    biginteger::BigInteger384,
    curves::{
        bls12_381::{
            g1::{Bls12_381G1Parameters, G1Affine, G1Projective},
            g2::{Bls12_381G2Parameters, G2Affine, G2Projective},
            Bls12_381,
        },
        models::SWModelParameters,
        tests::curve_tests,
        AffineCurve, PairingEngine, ProjectiveCurve,
    },
    fields::{
        bls12_381::{Fq, Fq12, Fq2, Fr},
        Field, PrimeField, SquareRootField,
    },
    groups::tests::group_test,
};
use rand;
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

#[test]
fn test_bilinearity() {
    let a: G1Projective = rand::random();
    let b: G2Projective = rand::random();
    let s: Fr = rand::random();

    let sa = a * &s;
    let sb = b * &s;

    let ans1 = Bls12_381::pairing(sa, b);
    let ans2 = Bls12_381::pairing(a, sb);
    let ans3 = Bls12_381::pairing(a, b).pow(s.into_repr());

    assert_eq!(ans1, ans2);
    assert_eq!(ans2, ans3);

    assert_ne!(ans1, Fq12::one());
    assert_ne!(ans2, Fq12::one());
    assert_ne!(ans3, Fq12::one());

    assert_eq!(ans1.pow(Fr::characteristic()), Fq12::one());
    assert_eq!(ans2.pow(Fr::characteristic()), Fq12::one());
    assert_eq!(ans3.pow(Fr::characteristic()), Fq12::one());
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
        rhs.add_assign(&Bls12_381G1Parameters::COEFF_B);

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

#[test]
fn test_g1_is_valid() {
    // Reject point on isomorphic twist (b = 24)
    {
        let p = G1Affine::new(
            Fq::from_repr(BigInteger384([
                0xc58d887b66c035dc,
                0x10cbfd301d553822,
                0xaf23e064f1131ee5,
                0x9fe83b1b4a5d648d,
                0xf583cc5a508f6a40,
                0xc3ad2aefde0bb13,
            ])),
            Fq::from_repr(BigInteger384([
                0x60aa6f9552f03aae,
                0xecd01d5181300d35,
                0x8af1cdb8aa8ce167,
                0xe760f57922998c9d,
                0x953703f5795a39e5,
                0xfe3ae0922df702c,
            ])),
            false,
        );
        assert!(!p.is_on_curve());
        assert!(p.is_in_correct_subgroup_assuming_on_curve());
    }

    // Reject point on a twist (b = 3)
    {
        let p = G1Affine::new(
            Fq::from_repr(BigInteger384([
                0xee6adf83511e15f5,
                0x92ddd328f27a4ba6,
                0xe305bd1ac65adba7,
                0xea034ee2928b30a8,
                0xbd8833dc7c79a7f7,
                0xe45c9f0c0438675,
            ])),
            Fq::from_repr(BigInteger384([
                0x3b450eb1ab7b5dad,
                0xa65cb81e975e8675,
                0xaa548682b21726e5,
                0x753ddf21a2601d20,
                0x532d0b640bd3ff8b,
                0x118d2c543f031102,
            ])),
            false,
        );
        assert!(!p.is_on_curve());
        assert!(!p.is_in_correct_subgroup_assuming_on_curve());
    }

    // Reject point in an invalid subgroup
    // There is only one r-order subgroup, as r does not divide the cofactor.
    {
        let p = G1Affine::new(
            Fq::from_repr(BigInteger384([
                0x76e1c971c6db8fe8,
                0xe37e1a610eff2f79,
                0x88ae9c499f46f0c0,
                0xf35de9ce0d6b4e84,
                0x265bddd23d1dec54,
                0x12a8778088458308,
            ])),
            Fq::from_repr(BigInteger384([
                0x8a22defa0d526256,
                0xc57ca55456fcb9ae,
                0x1ba194e89bab2610,
                0x921beef89d4f29df,
                0x5b6fda44ad85fa78,
                0xed74ab9f302cbe0,
            ])),
            false,
        );
        assert!(p.is_on_curve());
        assert!(!p.is_in_correct_subgroup_assuming_on_curve());
    }
}

#[test]
fn test_g1_addition_correctness() {
    let mut p = G1Projective::new(
        Fq::from_repr(BigInteger384([
            0x47fd1f891d6e8bbf,
            0x79a3b0448f31a2aa,
            0x81f3339e5f9968f,
            0x485e77d50a5df10d,
            0x4c6fcac4b55fd479,
            0x86ed4d9906fb064,
        ])),
        Fq::from_repr(BigInteger384([
            0xd25ee6461538c65,
            0x9f3bbb2ecd3719b9,
            0xa06fd3f1e540910d,
            0xcefca68333c35288,
            0x570c8005f8573fa6,
            0x152ca696fe034442,
        ])),
        Fq::one(),
    );

    p.add_assign(&G1Projective::new(
        Fq::from_repr(BigInteger384([
            0xeec78f3096213cbf,
            0xa12beb1fea1056e6,
            0xc286c0211c40dd54,
            0x5f44314ec5e3fb03,
            0x24e8538737c6e675,
            0x8abd623a594fba8,
        ])),
        Fq::from_repr(BigInteger384([
            0x6b0528f088bb7044,
            0x2fdeb5c82917ff9e,
            0x9a5181f2fac226ad,
            0xd65104c6f95a872a,
            0x1f2998a5a9c61253,
            0xe74846154a9e44,
        ])),
        Fq::one(),
    ));

    let p = G1Affine::from(p);

    assert_eq!(
        p,
        G1Affine::new(
            Fq::from_repr(BigInteger384([
                0x6dd3098f22235df,
                0xe865d221c8090260,
                0xeb96bb99fa50779f,
                0xc4f9a52a428e23bb,
                0xd178b28dd4f407ef,
                0x17fb8905e9183c69,
            ])),
            Fq::from_repr(BigInteger384([
                0xd0de9d65292b7710,
                0xf6a05f2bcf1d9ca7,
                0x1040e27012f20b64,
                0xeec8d1a5b7466c58,
                0x4bc362649dce6376,
                0x430cbdc5455b00a,
            ])),
            false,
        )
    );
}

#[test]
fn test_g1_doubling_correctness() {
    let mut p = G1Projective::new(
        Fq::from_repr(BigInteger384([
            0x47fd1f891d6e8bbf,
            0x79a3b0448f31a2aa,
            0x81f3339e5f9968f,
            0x485e77d50a5df10d,
            0x4c6fcac4b55fd479,
            0x86ed4d9906fb064,
        ])),
        Fq::from_repr(BigInteger384([
            0xd25ee6461538c65,
            0x9f3bbb2ecd3719b9,
            0xa06fd3f1e540910d,
            0xcefca68333c35288,
            0x570c8005f8573fa6,
            0x152ca696fe034442,
        ])),
        Fq::one(),
    );

    p.double_in_place();

    let p = G1Affine::from(p);

    assert_eq!(
        p,
        G1Affine::new(
            Fq::from_repr(BigInteger384([
                0xf939ddfe0ead7018,
                0x3b03942e732aecb,
                0xce0e9c38fdb11851,
                0x4b914c16687dcde0,
                0x66c8baf177d20533,
                0xaf960cff3d83833,
            ])),
            Fq::from_repr(BigInteger384([
                0x3f0675695f5177a8,
                0x2b6d82ae178a1ba0,
                0x9096380dd8e51b11,
                0x1771a65b60572f4e,
                0x8b547c1313b27555,
                0x135075589a687b1e,
            ])),
            false,
        )
    );
}

#[test]
fn test_g1_same_y() {
    // Test the addition of two points with different x coordinates
    // but the same y coordinate.

    // x1 = 128100205326445210408953809171070606737678357140298133325128175840781723996595026100005714405541449960643523234125
    // x2 = 3821408151224848222394078037104966877485040835569514006839342061575586899845797797516352881516922679872117658572470
    // y = 2291134451313223670499022936083127939567618746216464377735567679979105510603740918204953301371880765657042046687078

    let a = G1Affine::new(
        Fq::from_repr(BigInteger384([
            0xea431f2cc38fc94d,
            0x3ad2354a07f5472b,
            0xfe669f133f16c26a,
            0x71ffa8021531705,
            0x7418d484386d267,
            0xd5108d8ff1fbd6,
        ])),
        Fq::from_repr(BigInteger384([
            0xa776ccbfe9981766,
            0x255632964ff40f4a,
            0xc09744e650b00499,
            0x520f74773e74c8c3,
            0x484c8fc982008f0,
            0xee2c3d922008cc6,
        ])),
        false,
    );

    let b = G1Affine::new(
        Fq::from_repr(BigInteger384([
            0xe06cdb156b6356b6,
            0xd9040b2d75448ad9,
            0xe702f14bb0e2aca5,
            0xc6e05201e5f83991,
            0xf7c75910816f207c,
            0x18d4043e78103106,
        ])),
        Fq::from_repr(BigInteger384([
            0xa776ccbfe9981766,
            0x255632964ff40f4a,
            0xc09744e650b00499,
            0x520f74773e74c8c3,
            0x484c8fc982008f0,
            0xee2c3d922008cc6,
        ])),
        false,
    );

    // Expected
    // x = 52901198670373960614757979459866672334163627229195745167587898707663026648445040826329033206551534205133090753192
    // y = 1711275103908443722918766889652776216989264073722543507596490456144926139887096946237734327757134898380852225872709
    let c = G1Affine::new(
        Fq::from_repr(BigInteger384([
            0xef4f05bdd10c8aa8,
            0xad5bf87341a2df9,
            0x81c7424206b78714,
            0x9676ff02ec39c227,
            0x4c12c15d7e55b9f3,
            0x57fd1e317db9bd,
        ])),
        Fq::from_repr(BigInteger384([
            0x1288334016679345,
            0xf955cd68615ff0b5,
            0xa6998dbaa600f18a,
            0x1267d70db51049fb,
            0x4696deb9ab2ba3e7,
            0xb1e4e11177f59d4,
        ])),
        false,
    );

    assert!(a.is_on_curve() && a.is_in_correct_subgroup_assuming_on_curve());
    assert!(b.is_on_curve() && b.is_in_correct_subgroup_assuming_on_curve());
    assert!(c.is_on_curve() && c.is_in_correct_subgroup_assuming_on_curve());

    let mut tmp1 = a.into_projective();
    tmp1.add_assign(&b.into_projective());
    assert_eq!(tmp1.into_affine(), c);
    assert_eq!(tmp1, c.into_projective());

    let mut tmp2 = a.into_projective();
    tmp2.add_assign_mixed(&b);
    assert_eq!(tmp2.into_affine(), c);
    assert_eq!(tmp2, c.into_projective());
}

#[test]
fn test_g2_generator_raw() {
    let mut x = Fq2::zero();
    let mut i = 0;
    loop {
        // y^2 = x^3 + b
        let rhs = (x.square() * &x) + &Bls12_381G2Parameters::COEFF_B;
        if let Some(y) = rhs.sqrt() {
            let p = G2Affine::new(x, if y < -y { y } else { -y }, false);

            assert!(!p.is_in_correct_subgroup_assuming_on_curve());

            let g2 = p.scale_by_cofactor();
            if !g2.is_zero() {
                assert_eq!(i, 2);
                let g2 = G2Affine::from(g2);

                assert!(g2.is_in_correct_subgroup_assuming_on_curve());
                assert_eq!(g2, G2Affine::prime_subgroup_generator());
                break;
            }
        }

        i += 1;
        x += &Fq2::one();
    }
}

#[test]
fn test_g2_is_valid() {
    // Reject point on isomorphic twist (b = 3 * (u + 1))
    {
        let p = G2Affine::new(
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0xa757072d9fa35ba9,
                    0xae3fb2fb418f6e8a,
                    0xc1598ec46faa0c7c,
                    0x7a17a004747e3dbe,
                    0xcc65406a7c2e5a73,
                    0x10b8c03d64db4d0c,
                ])),
                Fq::from_repr(BigInteger384([
                    0xd30e70fe2f029778,
                    0xda30772df0f5212e,
                    0x5b47a9ff9a233a50,
                    0xfb777e5b9b568608,
                    0x789bac1fec71a2b9,
                    0x1342f02e2da54405,
                ])),
            ),
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0xfe0812043de54dca,
                    0xe455171a3d47a646,
                    0xa493f36bc20be98a,
                    0x663015d9410eb608,
                    0x78e82a79d829a544,
                    0x40a00545bb3c1e,
                ])),
                Fq::from_repr(BigInteger384([
                    0x4709802348e79377,
                    0xb5ac4dc9204bcfbd,
                    0xda361c97d02f42b2,
                    0x15008b1dc399e8df,
                    0x68128fd0548a3829,
                    0x16a613db5c873aaa,
                ])),
            ),
            false,
        );
        assert!(!p.is_on_curve());
        assert!(p.is_in_correct_subgroup_assuming_on_curve());
    }

    // Reject point on a twist (b = 2 * (u + 1))
    {
        let p = G2Affine::new(
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0xf4fdfe95a705f917,
                    0xc2914df688233238,
                    0x37c6b12cca35a34b,
                    0x41abba710d6c692c,
                    0xffcc4b2b62ce8484,
                    0x6993ec01b8934ed,
                ])),
                Fq::from_repr(BigInteger384([
                    0xb94e92d5f874e26,
                    0x44516408bc115d95,
                    0xe93946b290caa591,
                    0xa5a0c2b7131f3555,
                    0x83800965822367e7,
                    0x10cf1d3ad8d90bfa,
                ])),
            ),
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0xbf00334c79701d97,
                    0x4fe714f9ff204f9a,
                    0xab70b28002f3d825,
                    0x5a9171720e73eb51,
                    0x38eb4fd8d658adb7,
                    0xb649051bbc1164d,
                ])),
                Fq::from_repr(BigInteger384([
                    0x9225814253d7df75,
                    0xc196c2513477f887,
                    0xe05e2fbd15a804e0,
                    0x55f2b8efad953e04,
                    0x7379345eda55265e,
                    0x377f2e6208fd4cb,
                ])),
            ),
            false,
        );
        assert!(!p.is_on_curve());
        assert!(!p.is_in_correct_subgroup_assuming_on_curve());
    }

    // Reject point in an invalid subgroup
    // There is only one r-order subgroup, as r does not divide the cofactor.
    {
        let p = G2Affine::new(
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0x262cea73ea1906c,
                    0x2f08540770fabd6,
                    0x4ceb92d0a76057be,
                    0x2199bc19c48c393d,
                    0x4a151b732a6075bf,
                    0x17762a3b9108c4a7,
                ])),
                Fq::from_repr(BigInteger384([
                    0x26f461e944bbd3d1,
                    0x298f3189a9cf6ed6,
                    0x74328ad8bc2aa150,
                    0x7e147f3f9e6e241,
                    0x72a9b63583963fff,
                    0x158b0083c000462,
                ])),
            ),
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0x91fb0b225ecf103b,
                    0x55d42edc1dc46ba0,
                    0x43939b11997b1943,
                    0x68cad19430706b4d,
                    0x3ccfb97b924dcea8,
                    0x1660f93434588f8d,
                ])),
                Fq::from_repr(BigInteger384([
                    0xaaed3985b6dcb9c7,
                    0xc1e985d6d898d9f4,
                    0x618bd2ac3271ac42,
                    0x3940a2dbb914b529,
                    0xbeb88137cf34f3e7,
                    0x1699ee577c61b694,
                ])),
            ),
            false,
        );
        assert!(p.is_on_curve());
        assert!(!p.is_in_correct_subgroup_assuming_on_curve());
    }
}

#[test]
fn test_g2_addition_correctness() {
    let mut p = G2Projective::new(
        Fq2::new(
            Fq::from_repr(BigInteger384([
                0x6c994cc1e303094e,
                0xf034642d2c9e85bd,
                0x275094f1352123a9,
                0x72556c999f3707ac,
                0x4617f2e6774e9711,
                0x100b2fe5bffe030b,
            ])),
            Fq::from_repr(BigInteger384([
                0x7a33555977ec608,
                0xe23039d1fe9c0881,
                0x19ce4678aed4fcb5,
                0x4637c4f417667e2e,
                0x93ebe7c3e41f6acc,
                0xde884f89a9a371b,
            ])),
        ),
        Fq2::new(
            Fq::from_repr(BigInteger384([
                0xe073119472e1eb62,
                0x44fb3391fe3c9c30,
                0xaa9b066d74694006,
                0x25fd427b4122f231,
                0xd83112aace35cae,
                0x191b2432407cbb7f,
            ])),
            Fq::from_repr(BigInteger384([
                0xf68ae82fe97662f5,
                0xe986057068b50b7d,
                0x96c30f0411590b48,
                0x9eaa6d19de569196,
                0xf6a03d31e2ec2183,
                0x3bdafaf7ca9b39b,
            ])),
        ),
        Fq2::one(),
    );

    p.add_assign(&G2Projective::new(
        Fq2::new(
            Fq::from_repr(BigInteger384([
                0xa8c763d25910bdd3,
                0x408777b30ca3add4,
                0x6115fcc12e2769e,
                0x8e73a96b329ad190,
                0x27c546f75ee1f3ab,
                0xa33d27add5e7e82,
            ])),
            Fq::from_repr(BigInteger384([
                0x93b1ebcd54870dfe,
                0xf1578300e1342e11,
                0x8270dca3a912407b,
                0x2089faf462438296,
                0x828e5848cd48ea66,
                0x141ecbac1deb038b,
            ])),
        ),
        Fq2::new(
            Fq::from_repr(BigInteger384([
                0xf5d2c28857229c3f,
                0x8c1574228757ca23,
                0xe8d8102175f5dc19,
                0x2767032fc37cc31d,
                0xd5ee2aba84fd10fe,
                0x16576ccd3dd0a4e8,
            ])),
            Fq::from_repr(BigInteger384([
                0x4da9b6f6a96d1dd2,
                0x9657f7da77f1650e,
                0xbc150712f9ffe6da,
                0x31898db63f87363a,
                0xabab040ddbd097cc,
                0x11ad236b9ba02990,
            ])),
        ),
        Fq2::one(),
    ));

    let p = G2Affine::from(p);

    assert_eq!(
        p,
        G2Affine::new(
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0xcde7ee8a3f2ac8af,
                    0xfc642eb35975b069,
                    0xa7de72b7dd0e64b7,
                    0xf1273e6406eef9cc,
                    0xababd760ff05cb92,
                    0xd7c20456617e89,
                ])),
                Fq::from_repr(BigInteger384([
                    0xd1a50b8572cbd2b8,
                    0x238f0ac6119d07df,
                    0x4dbe924fe5fd6ac2,
                    0x8b203284c51edf6b,
                    0xc8a0b730bbb21f5e,
                    0x1a3b59d29a31274,
                ])),
            ),
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0x9e709e78a8eaa4c9,
                    0xd30921c93ec342f4,
                    0x6d1ef332486f5e34,
                    0x64528ab3863633dc,
                    0x159384333d7cba97,
                    0x4cb84741f3cafe8,
                ])),
                Fq::from_repr(BigInteger384([
                    0x242af0dc3640e1a4,
                    0xe90a73ad65c66919,
                    0x2bd7ca7f4346f9ec,
                    0x38528f92b689644d,
                    0xb6884deec59fb21f,
                    0x3c075d3ec52ba90,
                ])),
            ),
            false,
        )
    );
}

#[test]
fn test_g2_doubling_correctness() {
    let mut p = G2Projective::new(
        Fq2::new(
            Fq::from_repr(BigInteger384([
                0x6c994cc1e303094e,
                0xf034642d2c9e85bd,
                0x275094f1352123a9,
                0x72556c999f3707ac,
                0x4617f2e6774e9711,
                0x100b2fe5bffe030b,
            ])),
            Fq::from_repr(BigInteger384([
                0x7a33555977ec608,
                0xe23039d1fe9c0881,
                0x19ce4678aed4fcb5,
                0x4637c4f417667e2e,
                0x93ebe7c3e41f6acc,
                0xde884f89a9a371b,
            ])),
        ),
        Fq2::new(
            Fq::from_repr(BigInteger384([
                0xe073119472e1eb62,
                0x44fb3391fe3c9c30,
                0xaa9b066d74694006,
                0x25fd427b4122f231,
                0xd83112aace35cae,
                0x191b2432407cbb7f,
            ])),
            Fq::from_repr(BigInteger384([
                0xf68ae82fe97662f5,
                0xe986057068b50b7d,
                0x96c30f0411590b48,
                0x9eaa6d19de569196,
                0xf6a03d31e2ec2183,
                0x3bdafaf7ca9b39b,
            ])),
        ),
        Fq2::one(),
    );

    p.double_in_place();

    let p = G2Affine::from(p);

    assert_eq!(
        p,
        G2Affine::new(
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0x91ccb1292727c404,
                    0x91a6cb182438fad7,
                    0x116aee59434de902,
                    0xbcedcfce1e52d986,
                    0x9755d4a3926e9862,
                    0x18bab73760fd8024,
                ])),
                Fq::from_repr(BigInteger384([
                    0x4e7c5e0a2ae5b99e,
                    0x96e582a27f028961,
                    0xc74d1cf4ef2d5926,
                    0xeb0cf5e610ef4fe7,
                    0x7b4c2bae8db6e70b,
                    0xf136e43909fca0,
                ])),
            ),
            Fq2::new(
                Fq::from_repr(BigInteger384([
                    0x954d4466ab13e58,
                    0x3ee42eec614cf890,
                    0x853bb1d28877577e,
                    0xa5a2a51f7fde787b,
                    0x8b92866bc6384188,
                    0x81a53fe531d64ef,
                ])),
                Fq::from_repr(BigInteger384([
                    0x4c5d607666239b34,
                    0xeddb5f48304d14b3,
                    0x337167ee6e8e3cb6,
                    0xb271f52f12ead742,
                    0x244e6c2015c83348,
                    0x19e2deae6eb9b441,
                ])),
            ),
            false,
        )
    );
}
