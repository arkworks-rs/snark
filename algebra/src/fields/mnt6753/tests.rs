use crate::{
    bytes::ToBytes,
    fields::mnt6753::{Fq, Fq3, Fq3Parameters, Fq6, Fq6Parameters, FqParameters},
    fields::models::{Fp3Parameters, Fp6Parameters},
    fields::tests::{field_test, frobenius_test, primefield_test, sqrt_field_test},
    fields::FpParameters,
    to_bytes, BigInteger, BigInteger768, Field, PrimeField, SemanticallyValid, SquareRootField,
    ToBits, UniformRand,
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::{
    cmp::Ordering,
    ops::{AddAssign, MulAssign, SubAssign},
};

pub(crate) const ITERATIONS: usize = 5;

#[test]
fn test_mnt6753_fr() {
    use crate::fields::mnt6753::Fr;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let a: Fr = UniformRand::rand(&mut rng);
        let b: Fr = UniformRand::rand(&mut rng);
        field_test(a, b);
        primefield_test::<Fr>();
        sqrt_field_test(b);
    }
}

#[test]
fn test_mnt6753_fq() {
    use crate::fields::mnt6753::Fq;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let a: Fq = UniformRand::rand(&mut rng);
        let b: Fq = UniformRand::rand(&mut rng);
        field_test(a, b);
        primefield_test::<Fq>();
        sqrt_field_test(a);
    }
}

#[test]
fn test_mnt6753_fq3() {
    use crate::fields::mnt6753::{Fq, Fq3};

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let a: Fq3 = UniformRand::rand(&mut rng);
        let b: Fq3 = UniformRand::rand(&mut rng);
        field_test(a, b);
        sqrt_field_test(a);
    }
    frobenius_test::<Fq3, _>(Fq::characteristic(), 13);
}

#[test]
fn test_mnt6753_fq6() {
    use crate::fields::mnt6753::{Fq, Fq6};

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let g: Fq6 = UniformRand::rand(&mut rng);
        let h: Fq6 = UniformRand::rand(&mut rng);
        field_test(g, h);
    }
    frobenius_test::<Fq6, _>(Fq::characteristic(), 13);
}

#[test]
fn test_frob_coeffs() {
    //Fq3 coefficients tests
    let nqr = Fq::new(BigInteger768([
        0x4768931cfff9c7d4,
        0xc45e46d6ada96ca0,
        0x479b0bdb0b3c0107,
        0x362a089610f8d41b,
        0xdbafcec2c8a91aaf,
        0x78428b0ff9d96a06,
        0xf2e4472a9080c353,
        0xc9006ed33f0e971c,
        0x0794d9d10bdb7288,
        0x3c1e44cab5419e2c,
        0x49b5fc6c81f4560c,
        0x1c287777c30ba,
    ]));

    assert_eq!(Fq3Parameters::FROBENIUS_COEFF_FP3_C1[0], Fq::one());
    assert_eq!(
        Fq3Parameters::FROBENIUS_COEFF_FP3_C1[1],
        nqr.pow([
            0x9dad27a0c0000000,
            0x1a35885d0535be1a,
            0xf2412b3ebfaac7dd,
            0x3df5532766ddbe36,
            0x14f94386b8606212,
            0x90cec962ed733889,
            0x3345b6f335c8fd34,
            0xad54930ca2e04f2f,
            0x74e7f85324329d7b,
            0x3d5332701e853aa4,
            0x5ab6300ba4f6448f,
            0x96ecb9db96b0,
        ])
    );
    assert_eq!(
        Fq3Parameters::FROBENIUS_COEFF_FP3_C1[2],
        nqr.pow([
            0xeb5a4f4180000000,
            0x11ef14dc4d1d1a86,
            0xbcc65c98b6cd8a4d,
            0xaf02ccd5b8bbe4b6,
            0xe1c0ead251ed2f70,
            0xbd7f08fc1ab959c3,
            0x5ebf8d9546b52026,
            0xb931b7f5b7d71339,
            0x8215b63c8e7b5b88,
            0xf070a30875e8a2f3,
            0x57c26ad4f20a441e,
            0xb0f25e4ced223e74,
            0x5a202de8e5c64e7d,
            0x5ef5c9680032ac34,
            0xc8781f524d75840a,
            0x2ed25d673f721af8,
            0xe8a21855d90b0c39,
            0x6f239ff4c19bf775,
            0xa18f4849963ccc2d,
            0xe2aff6ba0cd7738a,
            0x1b4c5bfda65ab8de,
            0x651843808a3349ec,
            0x992ad2d0693fab2e,
            0x10aeece1d,
        ])
    );

    assert_eq!(Fq3Parameters::FROBENIUS_COEFF_FP3_C2[0], Fq::one());
    assert_eq!(
        Fq3Parameters::FROBENIUS_COEFF_FP3_C2[1],
        nqr.pow([
            0x3b5a4f4180000000,
            0x346b10ba0a6b7c35,
            0xe482567d7f558fba,
            0x7beaa64ecdbb7c6d,
            0x29f2870d70c0c424,
            0x219d92c5dae67112,
            0x668b6de66b91fa69,
            0x5aa9261945c09e5e,
            0xe9cff0a648653af7,
            0x7aa664e03d0a7548,
            0xb56c601749ec891e,
            0x12dd973b72d60,
        ])
    );

    assert_eq!(Fq3Parameters::FROBENIUS_COEFF_FP3_C2[0], Fq::one());
    assert_eq!(
        Fq3Parameters::FROBENIUS_COEFF_FP3_C2[2],
        nqr.pow([
            0xd6b49e8300000000,
            0x23de29b89a3a350d,
            0x798cb9316d9b149a,
            0x5e0599ab7177c96d,
            0xc381d5a4a3da5ee1,
            0x7afe11f83572b387,
            0xbd7f1b2a8d6a404d,
            0x72636feb6fae2672,
            0x42b6c791cf6b711,
            0xe0e14610ebd145e7,
            0xaf84d5a9e414883d,
            0x61e4bc99da447ce8,
            0xb4405bd1cb8c9cfb,
            0xbdeb92d000655868,
            0x90f03ea49aeb0814,
            0x5da4bace7ee435f1,
            0xd14430abb2161872,
            0xde473fe98337eeeb,
            0x431e90932c79985a,
            0xc55fed7419aee715,
            0x3698b7fb4cb571bd,
            0xca308701146693d8,
            0x3255a5a0d27f565c,
            0x215dd9c3b,
        ])
    );

    //Fq6 coefficients tests
    assert_eq!(Fq6Parameters::FROBENIUS_COEFF_FP6_C1[0], Fq::one());
    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C1[1],
        nqr.pow([
            0x4ed693d060000000,
            0x8d1ac42e829adf0d,
            0x7920959f5fd563ee,
            0x1efaa993b36edf1b,
            0x8a7ca1c35c303109,
            0x486764b176b99c44,
            0x99a2db799ae47e9a,
            0xd6aa498651702797,
            0x3a73fc2992194ebd,
            0x9ea999380f429d52,
            0x2d5b1805d27b2247,
            0x4b765cedcb58,
        ])
    );

    assert_eq!(
        Fq6Parameters::FROBENIUS_COEFF_FP6_C1[2],
        nqr.pow([
            0x75ad27a0c0000000,
            0x88f78a6e268e8d43,
            0x5e632e4c5b66c526,
            0x5781666adc5df25b,
            0xf0e0756928f697b8,
            0x5ebf847e0d5cace1,
            0xaf5fc6caa35a9013,
            0x5c98dbfadbeb899c,
            0xc10adb1e473dadc4,
            0x783851843af45179,
            0x2be1356a7905220f,
            0xd8792f2676911f3a,
            0x2d1016f472e3273e,
            0x2f7ae4b40019561a,
            0x643c0fa926bac205,
            0x97692eb39fb90d7c,
            0xf4510c2aec85861c,
            0xb791cffa60cdfbba,
            0x50c7a424cb1e6616,
            0x7157fb5d066bb9c5,
            0xda62dfed32d5c6f,
            0x328c21c04519a4f6,
            0xcc956968349fd597,
            0x8577670e,
        ])
    );

    let t: Vec<u64> = vec![
        0x7483bb7120000000,
        0x8d062427d1db0aa2,
        0x3b826fe19509943c,
        0x803af98e873186a1,
        0x44a537450bec8454,
        0xd5ada5dddb27253a,
        0x1606faafc48f607f,
        0xafb414db6cf23ea3,
        0x7ca96fe5870cf82f,
        0xf47533e58eeda80e,
        0x6ab085cda25e548,
        0xc70588a2f6aca3ec,
        0x4e7a82ea1fc07b79,
        0x239c19a96aa053df,
        0xc76424d2ce530d37,
        0x178682bb64cbad89,
        0x91eaf86f631b3025,
        0x61b790cb9c3b230f,
        0x8eab5a07d1bfc656,
        0xdbf7a1ce1d7b2d5a,
        0x9e635069d087e71d,
        0x1a2d4ef1d06912,
        0x2f95383190a1c6b1,
        0x4a5e0cd6bfab1545,
        0xce007f07298e8887,
        0x8b4c6084ce3db558,
        0xd8ff45c58e0350a0,
        0xfdb419b82e99aefe,
        0x2424e690e56a24f7,
        0x5d8ea3e32205d680,
        0xdfeb0fc2b1b9dc66,
        0xc2b28cd0856c8f63,
        0x1a6d183285a50e5a,
        0x4f2bd2076ce511bf,
        0x24125818ba511c50,
        0xec0e,
    ];

    assert_eq!(Fq6Parameters::FROBENIUS_COEFF_FP6_C1[3], nqr.pow(t));

    let t: Vec<u64> = vec![
        0x4b5a4f4180000000,
        0x783662c46a80572a,
        0x31a642537380b8be,
        0x31fdf7d2783357f9,
        0xa64bf80f29d4a380,
        0x68b8cf53ddd69d96,
        0x8e3c7b3ea0c0be37,
        0x36cb38751b0f77f1,
        0xa8a417be3370d640,
        0xf0c066c248813483,
        0x8c685a2ac4b83397,
        0xf00dbe827f89c59a,
        0x91168858a5f2a3,
        0x55c6534a8616a8a3,
        0xff9251ff410648a2,
        0xafea5a69638e6516,
        0xa5f4000f268eacad,
        0xafe8790437aabf46,
        0x53dce7024e46ccc,
        0x1493122a461ad407,
        0xfb431035eb5a2f3b,
        0xd0a6ea8220a35fca,
        0xec572e8a0b9fcadf,
        0x52e31949e769e68f,
        0x295a639156032388,
        0x80889685b9148b71,
        0x109af4f71e5eda89,
        0xd39e04c2971decbe,
        0x80c9a3088efc5d9a,
        0x32b7b5d36c4dde52,
        0x86ab1afba43c6854,
        0x8786fc69c6c89706,
        0x5c5e06e5c600b84d,
        0x6f8e49eff8d6266a,
        0x227824926819d14e,
        0xd53010c350ad1f3d,
        0x4a8383be03e79dbe,
        0x72510e49651532fa,
        0x53de452724de6bd2,
        0x5256bd99476fec9f,
        0xd4477261ad9e79de,
        0x36b00ab62f7320b7,
        0x4aa540e25039151c,
        0x14b46b50e1e196a,
        0xccf6320467271e3c,
        0x194918cffa199bad,
        0xa17fb4a61aee2779,
        0x1,
    ];

    assert_eq!(Fq6Parameters::FROBENIUS_COEFF_FP6_C1[4], nqr.pow(t));

    let t: Vec<u64> = vec![
        0xfa30e311e0000000,
        0x6ef817acd67e72da,
        0x45338206ce5c9133,
        0x8438e5505522ac6d,
        0x3ea99e081b654883,
        0xc58e2dba60f065e9,
        0xf9c736b998219ec2,
        0x8a7cada98eb6ad4b,
        0xfe8e8c18e441c138,
        0x56aa21cd75a7b0fc,
        0x29c788f106af584,
        0x3cc6d8eb83299dc6,
        0x46e1e69cc5500317,
        0x8406d92cc587504,
        0x453486afd6485d60,
        0x997ac2bc0efe8466,
        0xe0fc78cd0390a835,
        0x558cdc1d93d153e0,
        0x9459405791421871,
        0xce4b742bf235a9d0,
        0x113589f871521a06,
        0x936d759a69e7f120,
        0xab6259a1b5d8c1a9,
        0x7f6202ebec4e84cc,
        0x8efa2ad0205b0e50,
        0xf835c625d07a159e,
        0xdf446a2449285c34,
        0x864bc1a4e2f98422,
        0x79e1d8a346184e21,
        0x134ff19c36e84fe1,
        0x9ea41e9ba6f28d41,
        0xdcb1937741a9521f,
        0xcf1c4297ce1939d7,
        0xb1b593ead6925afa,
        0xd345ad63dc77da8,
        0x1e4351c475baf3d5,
        0x5c7dd5e318f30a66,
        0xc6ce35c17727f901,
        0x69033f0d8fd3d637,
        0xbcbcc848d6ff231b,
        0xbb7d7610675d36bd,
        0x2d3b28aeeb42d9d5,
        0x4124fa3a31712e5a,
        0x56c4bb63c9ebdc47,
        0x6a264b5bbbbabba6,
        0x565a57cb6f363b87,
        0x6037930dbb8d7260,
        0x5f6a5cc2377d7ad5,
        0xcc5b55c23062df09,
        0xcd8679e41a58f4b8,
        0x5b611a47f5cdf206,
        0x245fce8964046f4d,
        0x69671932dd937705,
        0xe1c838f724cee452,
        0x1605d803ca616eac,
        0x834a9c5bf260a369,
        0xd5463104b340a76e,
        0x25761608b8cc8150,
        0x2e268ae0dcd5d,
    ];

    assert_eq!(Fq6Parameters::FROBENIUS_COEFF_FP6_C1[5], nqr.pow(t));
}

#[test]
fn test_neg_one() {
    let neg_one = Fq::new(BigInteger768([
        0x1f70f6cdc00090bf,
        0xffef30ff5a176ba7,
        0x34d7aee33286761d,
        0xaa6d9cc76f4f79ca,
        0x93df7bad553a4b62,
        0x12afb31fea4cde38,
        0x67c4e9228e277305,
        0xae72762315b0e32b,
        0x1e431f2d6f0b3251,
        0xa85518752329c761,
        0x7add306fcee92c18,
        0x1497e8ec9e1ce,
    ]));
    assert_eq!(neg_one, -Fq::one());
}

#[test]
fn test_fq_is_valid() {
    let mut a = Fq::new(FqParameters::MODULUS);
    assert!(!a.is_valid());
    a.0.sub_noborrow(&BigInteger768::from(1));
    assert!(a.is_valid());
    assert!(Fq::new(BigInteger768::from(0)).is_valid());
    assert!(Fq::new(BigInteger768([
        4334945402112658761,
        5754338769440963294,
        213982681521013032,
        12086861433024916758,
        3907127509866713521,
        13945672019815712008,
        15986918099604157897,
        1610539633786093561,
        6468823346244563772,
        2229132487154770553,
        4774597994323744797,
        467097584308943,
    ]))
    .is_valid());

    assert!(!Fq::new(BigInteger768([
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
        0xffffffffffffffff,
    ]))
    .is_valid());

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        let a = Fq::rand(&mut rng);
        assert!(a.is_valid());
    }
}

#[test]
fn test_fq_add_assign() {
    {
        // Random number
        let mut tmp = Fq::new(BigInteger768([
            0xd71ecd9366caebf4,
            0xf79dc7c215b5603d,
            0xc48a95ce7d901364,
            0x37be56f69d9fd17b,
            0x309e25a2d9c8a0d5,
            0xb7fedecd78bc2e5b,
            0x56617000c0209a8e,
            0x398e53fdc65be214,
            0x2d36e4a81de6695e,
            0x6999415ac5be48f4,
            0x8c8376b88bb61a43,
            0x571a8191ad3f,
        ]));
        assert!(tmp.is_valid());
        // Test that adding zero has no effect.
        tmp.add_assign(&Fq::new(BigInteger768::from(0)));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0xd71ecd9366caebf4,
                0xf79dc7c215b5603d,
                0xc48a95ce7d901364,
                0x37be56f69d9fd17b,
                0x309e25a2d9c8a0d5,
                0xb7fedecd78bc2e5b,
                0x56617000c0209a8e,
                0x398e53fdc65be214,
                0x2d36e4a81de6695e,
                0x6999415ac5be48f4,
                0x8c8376b88bb61a43,
                0x571a8191ad3f,
            ]))
        );
        // Add one and test for the result.
        tmp.add_assign(&Fq::new(BigInteger768::from(1)));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0xd71ecd9366caebf5,
                0xf79dc7c215b5603d,
                0xc48a95ce7d901364,
                0x37be56f69d9fd17b,
                0x309e25a2d9c8a0d5,
                0xb7fedecd78bc2e5b,
                0x56617000c0209a8e,
                0x398e53fdc65be214,
                0x2d36e4a81de6695e,
                0x6999415ac5be48f4,
                0x8c8376b88bb61a43,
                0x571a8191ad3f,
            ]))
        );
        // Add another random number that exercises the reduction.
        tmp.add_assign(&Fq::new(BigInteger768([
            0xa833ac17d3ed787d,
            0x51b38b65f01cf4b,
            0x6f9726cd804f9fba,
            0xcee05b70be5e63e8,
            0x96e8e72359ce53fa,
            0xf0424c5bce9ffe0e,
            0xb88205d99d8c566c,
            0x51f7b750d7c59b83,
            0xec186e8312be7259,
            0x9a575ff427f8e3c6,
            0xb63befd8b7eb976,
            0x8328f53786bc,
        ])));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0x7f5279ab3ab86472,
                0xfcb9007874b72f89,
                0x3421bc9bfddfb31e,
                0x69eb2675bfe3564,
                0xc7870cc63396f4d0,
                0xa8412b29475c2c69,
                0xee375da5dacf0fb,
                0x8b860b4e9e217d98,
                0x194f532b30a4dbb7,
                0x3f0a14eedb72cbb,
                0x97e735b61734d3ba,
                0xda4376c933fb,
            ]))
        );
        // Add one to (q - 1) and test for the result.
        tmp = Fq::new(BigInteger768([
            0xD90776E240000000,
            0x4EA099170FA13A4F,
            0xD6C381BC3F005797,
            0xB9DFF97634993AA4,
            0x3EEBCA9429212636,
            0xB26C5C28C859A99B,
            0x99D124D9A15AF79D,
            0x07FDB925E8A0ED8D,
            0x5EB7E8F96C97D873,
            0xB7F997505B8FAFED,
            0x10229022EEE2CDAD,
            0x01C4C62D92C411,
        ]));
        tmp.add_assign(&Fq::new(BigInteger768::from(1)));
        assert!(tmp.0.is_zero());
        // Add a random number to another one such that the result is q - 1
        tmp = Fq::new(BigInteger768([
            0x87abddbd74a24e95,
            0x502e6d2817124c21,
            0x73f2fa5c2607b48a,
            0xf1f4db12bee0c496,
            0xfbb608395fba6b77,
            0xe4aaf10d8f0e93c1,
            0x6b1e008da5bcb156,
            0xc342430f055fe8b1,
            0x964c5700528ae61e,
            0xf83cce9c272e4c22,
            0x5c90f71b9258fa9f,
            0x1a6ca658e16d7,
        ]));
        tmp.add_assign(&Fq::new(BigInteger768([
            0x515b9924cb5db16b,
            0xfe722beef88eee2e,
            0x62d0876018f8a30c,
            0xc7eb1e6375b8760e,
            0x4335c25ac966babe,
            0xcdc16b1b394b15d9,
            0x2eb3244bfb9e4646,
            0x44bb7616e34104dc,
            0xc86b91f91a0cf254,
            0xbfbcc8b4346163ca,
            0xb39199075c89d30d,
            0x1dfbc804ad39,
        ])));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0xD90776E240000000,
                0x4EA099170FA13A4F,
                0xD6C381BC3F005797,
                0xB9DFF97634993AA4,
                0x3EEBCA9429212636,
                0xB26C5C28C859A99B,
                0x99D124D9A15AF79D,
                0x07FDB925E8A0ED8D,
                0x5EB7E8F96C97D873,
                0xB7F997505B8FAFED,
                0x10229022EEE2CDAD,
                0x01C4C62D92C411,
            ]))
        );
        // Add one to the result and test for it.
        tmp.add_assign(&Fq::new(BigInteger768::from(1)));
        assert!(tmp.0.is_zero());
    }

    // Test associativity

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        // Generate a, b, c and ensure (a + b) + c == a + (b + c).
        let a = Fq::rand(&mut rng);
        let b = Fq::rand(&mut rng);
        let c = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.add_assign(&b);
        tmp1.add_assign(&c);

        let mut tmp2 = b;
        tmp2.add_assign(&c);
        tmp2.add_assign(&a);

        assert!(tmp1.is_valid());
        assert!(tmp2.is_valid());
        assert_eq!(tmp1, tmp2);
    }
}

#[test]
fn test_fq_sub_assign() {
    {
        // Test arbitrary subtraction that tests reduction.
        let mut tmp = Fq::new(BigInteger768([
            0x498e02ae4388631c,
            0x4e46e93ca8740be4,
            0x2c375ca1f4ce59da,
            0xb6976fbf66abf4d8,
            0xc44700335d60a831,
            0xb98f3093987eb2b6,
            0xa1ae604739e8506b,
            0x9ee466288cc528c1,
            0xe4f8ae389ba8678b,
            0x7f8f04daec7d211b,
            0xccb0bd2a806c61a0,
            0x7799ba9ec5f9,
        ]));
        tmp.sub_assign(&Fq::new(BigInteger768([
            0x8c27de877e6dd70,
            0xc4f112f5ff1602de,
            0x2e61440d01e26550,
            0xd7a5a7a27f1dd489,
            0x18c873d5908c3e1b,
            0xbff6fd4ee069d1c9,
            0x7d27a79b256e5544,
            0xdc6484f9c2f7d42d,
            0xfc45abe1c7dcc1b8,
            0x1ddc807c50f98529,
            0x3a5c2cbda6a1c247,
            0x24a8b7761c01,
        ])));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0x40cb84c5cba185ac,
                0x8955d646a95e0906,
                0xfdd61894f2ebf489,
                0xdef1c81ce78e204e,
                0xab7e8c5dccd46a15,
                0xf9983344b814e0ed,
                0x2486b8ac1479fb26,
                0xc27fe12ec9cd5494,
                0xe8b30256d3cba5d2,
                0x61b2845e9b839bf1,
                0x9254906cd9ca9f59,
                0x52f10328a9f8,
            ]))
        );

        // Test the opposite subtraction which doesn't test reduction.
        tmp = Fq::new(BigInteger768([
            0x8c27de877e6dd70,
            0xc4f112f5ff1602de,
            0x2e61440d01e26550,
            0xd7a5a7a27f1dd489,
            0x18c873d5908c3e1b,
            0xbff6fd4ee069d1c9,
            0x7d27a79b256e5544,
            0xdc6484f9c2f7d42d,
            0xfc45abe1c7dcc1b8,
            0x1ddc807c50f98529,
            0x3a5c2cbda6a1c247,
            0x24a8b7761c01,
        ]));
        tmp.sub_assign(&Fq::new(BigInteger768([
            0x498e02ae4388631c,
            0x4e46e93ca8740be4,
            0x2c375ca1f4ce59da,
            0xb6976fbf66abf4d8,
            0xc44700335d60a831,
            0xb98f3093987eb2b6,
            0xa1ae604739e8506b,
            0x9ee466288cc528c1,
            0xe4f8ae389ba8678b,
            0x7f8f04daec7d211b,
            0xccb0bd2a806c61a0,
            0x7799ba9ec5f9,
        ])));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0x983bf21c745e7a55,
                0xc54ac2d066433149,
                0xd8ed69274c14630d,
                0xdaee31594d0b1a55,
                0x936d3e365c4cbc20,
                0xb8d428e41044c8ad,
                0x754a6c2d8ce0fc76,
                0x457dd7f71ed398f9,
                0x7604e6a298cc32a0,
                0x564712f1c00c13fb,
                0x7dcdffb615182e54,
                0x171d52a6a1a18,
            ]))
        );

        // Test for sensible results with zero
        tmp = Fq::new(BigInteger768::from(0));
        tmp.sub_assign(&Fq::new(BigInteger768::from(0)));
        assert!(tmp.is_zero());

        tmp = Fq::new(BigInteger768([
            0xf9566119e604dffd,
            0x5128b45b17b1d1a,
            0xf9c68449e46773bc,
            0x83299ded47cf7c3c,
            0xf012a23810f90638,
            0x7764504372af51f6,
            0xa8c4872246f123d7,
            0xb19941f9e9b8ebe6,
            0x5ca77716ab16ce47,
            0x1afadccedb97f7c9,
            0x8e7b2a674d82952b,
            0x1c34b30dad512,
        ]));
        tmp.sub_assign(&Fq::new(BigInteger768::from(0)));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0xf9566119e604dffd,
                0x5128b45b17b1d1a,
                0xf9c68449e46773bc,
                0x83299ded47cf7c3c,
                0xf012a23810f90638,
                0x7764504372af51f6,
                0xa8c4872246f123d7,
                0xb19941f9e9b8ebe6,
                0x5ca77716ab16ce47,
                0x1afadccedb97f7c9,
                0x8e7b2a674d82952b,
                0x1c34b30dad512,
            ]))
        );
    }

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        // Ensure that (a - b) + (b - a) = 0.
        let a = Fq::rand(&mut rng);
        let b = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.sub_assign(&b);

        let mut tmp2 = b;
        tmp2.sub_assign(&a);

        tmp1.add_assign(&tmp2);
        assert!(tmp1.is_zero());
    }
}

#[test]
fn test_fq_mul_assign() {
    let mut tmp = Fq::new(BigInteger768([
        0xd13c2d67ff927355,
        0x48a3ed07a9255a92,
        0x982dc83dc987c355,
        0x95101640f58d0561,
        0x8c3b1401350d9d04,
        0x57cfe2bcc2bcea3e,
        0xbe3d2c34db7d9eca,
        0xe177d2eeaa96423b,
        0xce933a42b85f0ffb,
        0xab63f8cfb8d7325d,
        0xd7bbaa9d7be4f1e2,
        0x1b93ec2058d83,
    ]));
    tmp.mul_assign(&Fq::new(BigInteger768([
        0x4c91701e19e27ed,
        0x60f7ac0ade842b38,
        0xe21eaa804f7da8f5,
        0xd00d4726be1ec9d9,
        0x5cd2ce1f29861e99,
        0x2955e6fdee8d2147,
        0x57550b212d00d16e,
        0x14ad5270ee634ea2,
        0x11261b2c50ab3ce5,
        0x243adeb987ff1db3,
        0xed362b88a41ff6eb,
        0x13a7c7a59137,
    ])));
    assert_eq!(
        tmp,
        Fq::new(BigInteger768([
            0x42d9ca5556458579,
            0x27a6b6355cbe1d0a,
            0xc903f28df9887807,
            0xca107724944437d0,
            0xc89e7003b7975cb8,
            0x9a710305830486ef,
            0xcc5e827d897dd989,
            0x514efb8818c8f976,
            0x1dc8cff9cdef10ee,
            0x1fd9341478822f4d,
            0x528e7ff14e1b3445,
            0xfe7cea18aace,
        ]))
    );

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000000 {
        // Ensure that (a * b) * c = a * (b * c)
        let a = Fq::rand(&mut rng);
        let b = Fq::rand(&mut rng);
        let c = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.mul_assign(&b);
        tmp1.mul_assign(&c);

        let mut tmp2 = b;
        tmp2.mul_assign(&c);
        tmp2.mul_assign(&a);

        assert_eq!(tmp1, tmp2);
    }

    for _ in 0..1000000 {
        // Ensure that r * (a + b + c) = r*a + r*b + r*c

        let r = Fq::rand(&mut rng);
        let mut a = Fq::rand(&mut rng);
        let mut b = Fq::rand(&mut rng);
        let mut c = Fq::rand(&mut rng);

        let mut tmp1 = a;
        tmp1.add_assign(&b);
        tmp1.add_assign(&c);
        tmp1.mul_assign(&r);

        a.mul_assign(&r);
        b.mul_assign(&r);
        c.mul_assign(&r);

        a.add_assign(&b);
        a.add_assign(&c);

        assert_eq!(tmp1, a);
    }
}

#[test]
fn test_fq_squaring() {
    let mut a = Fq::new(BigInteger768([
        0x9584b092cee7ad6e,
        0x612640abc76c7a89,
        0xa149278cf5a0dbf4,
        0x6cca1f556c978932,
        0xaa2147d8f6d82c95,
        0xec1ced81b01d59ae,
        0x3d26125f76e2e9c,
        0xd4e4ca0a8a03f769,
        0x1627eb466dac06f0,
        0x7e41fd42d3d1e7fe,
        0xb5cc058c94c33a12,
        0x304222fb7e55,
    ]));
    assert!(a.is_valid());
    a.square_in_place();
    assert_eq!(
        a,
        Fq::from_repr(BigInteger768([
            0x8508bbd8f2bc49c1,
            0x3df121fd16e9636,
            0xa91b04bd08bf3e3e,
            0xcfac2ec8a861466a,
            0x364a448d4574e219,
            0x163623609ca88080,
            0x4a645805dc7788be,
            0xf6da93019b94efe9,
            0x2456f275e0477c59,
            0x886cb05190f97e75,
            0xee362a46476f91de,
            0xe066089e656c,
        ]))
    );

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000000 {
        // Ensure that (a * a) = a^2
        let a = Fq::rand(&mut rng);

        let mut tmp = a;
        tmp.square_in_place();

        let mut tmp2 = a;
        tmp2.mul_assign(&a);

        assert_eq!(tmp, tmp2);
    }
}

#[test]
fn test_fq_inverse() {
    assert!(Fq::zero().inverse().is_none());

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let one = Fq::one();

    for _ in 0..1000 {
        // Ensure that a * a^-1 = 1
        let mut a = Fq::rand(&mut rng);
        let ainv = a.inverse().unwrap();
        a.mul_assign(&ainv);
        assert_eq!(a, one);
    }
}

#[test]
fn test_fq_double_in_place() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        // Ensure doubling a is equivalent to adding a to itself.
        let mut a = Fq::rand(&mut rng);
        let mut b = a;
        b.add_assign(&a);
        a.double_in_place();
        assert_eq!(a, b);
    }
}

#[test]
fn test_fq_negate() {
    {
        let a = -Fq::zero();

        assert!(a.is_zero());
    }

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        // Ensure (a - (-a)) = 0.
        let mut a = Fq::rand(&mut rng);
        let b = -a;
        a.add_assign(&b);

        assert!(a.is_zero());
    }
}

#[test]
fn test_fq_pow() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for i in 0..1000 {
        // Exponentiate by various small numbers and ensure it consists with repeated
        // multiplication.
        let a = Fq::rand(&mut rng);
        let target = a.pow(&[i]);
        let mut c = Fq::one();
        for _ in 0..i {
            c.mul_assign(&a);
        }
        assert_eq!(c, target);
    }

    for _ in 0..1000 {
        // Exponentiating by the modulus should have no effect in a prime field.
        let a = Fq::rand(&mut rng);

        assert_eq!(a, a.pow(Fq::characteristic()));
    }
}

#[test]
fn test_fq_sqrt() {
    let a_squared = Fq::new(BigInteger768([
        0x815f0a0b6846238c,
        0x5949c2aef4191aac,
        0x7dd3ce5c3e2aca9b,
        0x33626ad4f94ccca5,
        0xef5d495e9555b4ff,
        0x8414d5bd8de49ef4,
        0x5b16424c19676079,
        0x58a9a5ebe1c0a51f,
        0xe5fd980f6f8e1385,
        0x4dab45076384cd54,
        0x5fb2e86d6a180aa0,
        0x18cd3d7d0d32c,
    ]));
    let a = a_squared.sqrt().unwrap();
    assert_eq!(
        a,
        Fq::new(BigInteger768([
            0x9696cbd70f683946,
            0x2fbe984a7f99fb5e,
            0x4152df84cdbbce30,
            0x5ebe5c628d5a355e,
            0xd448a5598db83394,
            0x5c8ba2124ebab55b,
            0xe8ef67a51612340,
            0x3c7382380f2f323f,
            0x93e740b12eec5af3,
            0x89a625f9546373f2,
            0xb2c2be3f0e9c1b51,
            0xafc892d9432f,
        ]))
    );

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    assert_eq!(Fq::zero().sqrt().unwrap(), Fq::zero());

    for _ in 0..1000 {
        // Ensure sqrt(a^2) = a or -a
        let a = Fq::rand(&mut rng);
        let nega = -a;
        let mut b = a;
        b.square_in_place();

        let b = b.sqrt().unwrap();

        assert!(a == b || nega == b);
    }

    for _ in 0..1000 {
        // Ensure sqrt(a)^2 = a for random a
        let a = Fq::rand(&mut rng);

        if let Some(mut tmp) = a.sqrt() {
            tmp.square_in_place();

            assert_eq!(a, tmp);
        }
    }
}

#[test]
fn test_fq_num_bits() {
    assert_eq!(FqParameters::MODULUS_BITS, 753);
    assert_eq!(FqParameters::CAPACITY, 752);
}

#[test]
fn test_fq_root_of_unity() {
    assert_eq!(FqParameters::TWO_ADICITY, 30);
    assert_eq!(
        Fq::multiplicative_generator(),
        Fq::from_repr(BigInteger768::from(17))
    );
    assert_eq!(
        Fq::multiplicative_generator().pow([
            0x3E84E93F641DDB89,
            0xFC015E5D3A82645C,
            0xD264EA935B0E06F0,
            0xA48498DAE77FE5D8,
            0x2166A66CFBAF2A50,
            0x856BDE76C9B170A3,
            0xA283B63667449366,
            0xB25F61CC1FF6E497,
            0x6E3EBFB57ADFA3E5,
            0xBB8B36B6DFE65D41,
            0xB64B1044408A408B,
            0x71318,
        ]),
        Fq::root_of_unity()
    );
    assert_eq!(
        Fq::root_of_unity().pow([1 << FqParameters::TWO_ADICITY]),
        Fq::one()
    );
    assert!(Fq::multiplicative_generator().sqrt().is_none());
}

#[test]
fn test_fq_ordering() {
    // BigInteger768's ordering is well-tested, but we still need to make sure the
    // Fq elements aren't being compared in Montgomery form.
    for i in 0..100 {
        assert!(Fq::from_repr(BigInteger768::from(i + 1)) > Fq::from_repr(BigInteger768::from(i)));
    }
}

#[test]
fn test_fq_legendre() {
    use crate::fields::LegendreSymbol::*;

    assert_eq!(QuadraticResidue, Fq::one().legendre());
    assert_eq!(Zero, Fq::zero().legendre());

    assert_eq!(
        QuadraticNonResidue,
        Fq::new(BigInteger768::from(11)).legendre()
    );
    assert_eq!(
        QuadraticResidue,
        Fq::new(BigInteger768::from(121)).legendre()
    );

    let e = BigInteger768([
        0xcd3165e89d5a6aec,
        0xfe794ff2c42eb6e2,
        0xb0fca85ac55700a1,
        0x436d3571cec6563a,
        0x21f200cb4d7bcfeb,
        0xc9ef4b7f2c8e322f,
        0x47ee13c44f7f103c,
        0xa90deb8adf454d57,
        0xd46263471ead205d,
        0xb5c0ca287e04c8ef,
        0x6d39092796528ed0,
        0x4ded2b341ce9,
    ]);
    assert_eq!(QuadraticNonResidue, Fq::from_repr(e).legendre());
    let e = BigInteger768([
        0x3804b3df0abf369c,
        0x5e5c63ebc1f00fb5,
        0x7d48d4105cbc998b,
        0x78139f68d36ffea2,
        0x76bb92a44562f42d,
        0xd7e07480810fac6d,
        0xf8341b3074bba776,
        0xf8fd6de7df4decaf,
        0xa5387dd904bc8425,
        0x590b46d9358adc1b,
        0x7d72acdc191ac04f,
        0xbd6fdb57393b,
    ]);
    assert_eq!(QuadraticResidue, Fq::from_repr(e).legendre());
}

#[test]
fn test_fq_bytes() {
    let a = Fq::from_repr(BigInteger768([
        0xecb23e0a2f1bd245,
        0x31b80c631edc39ae,
        0x590554649da78d85,
        0xc12a5b2366697d8c,
        0x1c4d1238bc71b6b2,
        0x8d0b219cdb2c21d4,
        0xa4a93f04ec97df4a,
        0x67ca57a5e41af3f0,
        0xd8ac205142b50568,
        0xc2a674c3ef59f6de,
        0x40067f9930800bb5,
        0x129bc5bbcc439,
    ]));
    let a_b = to_bytes!(a).unwrap();
    let a_b_read = std::fs::read("src/fields/mnt6753/test_vec/mnt6753_tobyte").unwrap();
    assert_eq!(a_b, a_b_read);
}

#[test]
fn test_convert_fq_fr() {
    use crate::fields::{
        convert,
        mnt6753::{FqParameters, Fr},
    };

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        // Safely convert a random Fq into a Fr
        let q: Fq = UniformRand::rand(&mut rng);
        let q_bits = &q.write_bits()[1..]; //Skip 1 bit, in order to perform a safe conversion
        let conv = convert::<Fr>(q_bits.to_vec()).unwrap();
        assert_eq!(conv.pow(Fr::characteristic()), conv);

        // Safely convert a random Fr into a Fq
        let r: Fr = UniformRand::rand(&mut rng);
        let r_bits = &r.write_bits()[1..]; //Skip 1 bit, in order to perform a safe conversion
        let conv = convert::<Fq>(r_bits.to_vec()).unwrap();
        assert_eq!(conv.pow(Fq::characteristic()), conv);
    }

    //Attempting to convert a bit array that exceeds other field's modulus will result in an error
    let modulus_q = Fq::new(FqParameters::MODULUS);
    assert!(convert::<Fr>((modulus_q - &Fq::one()).write_bits()).is_err()); //Fq_Modulus - 1 is bigger than Fr modulus
}

#[test]
fn test_fq3_ordering() {
    let mut a = Fq3::new(Fq::zero(), Fq::zero(), Fq::zero());

    let mut b = a.clone();

    assert!(a.cmp(&b) == Ordering::Equal);
    b.c0.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Less);
    a.c0.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Equal);
    b.c1.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Less);
    a.c0.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Less);
    a.c1.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Greater);
    b.c0.add_assign(&Fq::one());
    assert!(a.cmp(&b) == Ordering::Equal);
}

#[test]
fn test_fq3_basics() {
    assert_eq!(Fq3::new(Fq::zero(), Fq::zero(), Fq::zero()), Fq3::zero());
    assert_eq!(Fq3::new(Fq::one(), Fq::zero(), Fq::zero()), Fq3::one());
    assert!(Fq3::zero().is_zero());
    assert!(!Fq3::one().is_zero());
    assert!(!Fq3::new(Fq::zero(), Fq::one(), Fq::zero()).is_zero());
}

#[test]
fn test_fq3_squaring() {
    // i = sqrt(11) in mnt6_753 fq3

    //(11 + 0i + 0j)^2 = 121 + 0i + 0j = 121
    let a = Fq3::new(
        Fq::from_repr(BigInteger768::from(11)),
        Fq::zero(),
        Fq::zero(),
    )
    .square();
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768::from(121)),
            Fq::zero(),
            Fq::zero()
        )
    );

    let mut a = Fq3::new(
        Fq::from_repr(BigInteger768([
            0x4021ea91c1d7ed85,
            0x78bd65a80e548b88,
            0xb15b43140075906,
            0x143f963245fdce5f,
            0x3b07cd86d0be3997,
            0xbc8ccddb928234fc,
            0x80c7ce910c156ccc,
            0x117349672833c222,
            0x949d3f355fb23ca3,
            0x39c56800c5cb7865,
            0x33bcb002a225b9cb,
            0xcfdcc84af403,
        ])),
        Fq::from_repr(BigInteger768([
            0x9f5f269a1851cdc,
            0xbbe507b0cd3c14f9,
            0xaf9b1ffd3d0bc566,
            0x75321f979bd32876,
            0x8da033b1cebcfefc,
            0x26e23be5fd1d73c,
            0x761d56f0d23101f6,
            0x96bd1d6551c00c65,
            0x5e17d1b99b7a938,
            0x149c5088da3e94d1,
            0xb4771eea63454794,
            0x1a2b05501eeac,
        ])),
        Fq::from_repr(BigInteger768([
            0xe3e78d907b757a8c,
            0x17c324937f5735f2,
            0x5a0e6bc2ced16f28,
            0xb37843eece4d27b3,
            0xb8b225093ee90d4a,
            0xf80398bf27c851cc,
            0xb59f5ad72ad7b97,
            0x4f1e23b8c8aba546,
            0xe1bb6e79103f6e5,
            0xc6030187e4a8b103,
            0x49b7488076cb6d63,
            0x3cd9059f596d,
        ])),
    );
    a.square_in_place();
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0xf88598c734a0ebea,
                0xc4d16ff251e0d4ea,
                0xfdc0993bea96656b,
                0x42d66c4b9fb4c54c,
                0x895ad48fed7ebc9b,
                0x7c356895984e4190,
                0x8e33757f7dbe3f2f,
                0x8bf91ab0fb58265,
                0x31398876d267714e,
                0x60b52f79cb1b6ef2,
                0x9dac6bdcc46ffcca,
                0x16971309611df,
            ])),
            Fq::from_repr(BigInteger768([
                0xf3421c6cf4aaefc2,
                0x53c6bf7bbef81e58,
                0x8e7cb542c36ebc1d,
                0xbdb4b657dbd84fe2,
                0xb8e7b475066e8b42,
                0x1c0c3e2ec27607f,
                0x96c04d86f1e68651,
                0xeb63b5ec8c258fd1,
                0x79adf6a1e68f0cfb,
                0x3ac7dd569185f706,
                0x4aed0fb556acdcfe,
                0x5b97eb0ec575,
            ])),
            Fq::from_repr(BigInteger768([
                0xe941d8679493cd4c,
                0x32f1a39375bfb2f8,
                0x146259a5646c77f1,
                0xc0a32a02fa2b005a,
                0xa8f3bbf89b5b7ca4,
                0xd40dec29cbed7259,
                0xfd621ece49ba14db,
                0xcafd071c105f9545,
                0xa3f7d1d45413c2a9,
                0x3365e865115b992,
                0x2703e89cb154b091,
                0x2750df057935,
            ])),
        )
    );
}

#[test]
fn test_fq3_mul() {
    let mut a = Fq3::new(
        Fq::from_repr(BigInteger768([
            0x3246ad9d5f2e162,
            0xb1a35a4ae6524cf4,
            0x43c67ef1fc963fe9,
            0xbdd7d3a9279e5133,
            0xf674346a7c629056,
            0x4210f13baa1290d3,
            0xcb3b4dca901c3526,
            0x1e1531abf7129c57,
            0x5e28ee42c19296e8,
            0x30ab88ba97c1c54b,
            0xeab7ccc94a470cac,
            0x19c5d753734e8,
        ])),
        Fq::from_repr(BigInteger768([
            0x47a591f3099b6e84,
            0x17aa7b9c0701cb9c,
            0x7afb9e4c0bae4299,
            0x6e6c706d1e24a31f,
            0xd71d5501ca9eb295,
            0x209e568e6d68be58,
            0xf414daf62310083e,
            0x6d7d60ae249d9902,
            0xb12a28d2d941f691,
            0x30c08fde4dc9e7cd,
            0x9c2cf5cb5a054f10,
            0xd8ae87a6c6f9,
        ])),
        Fq::from_repr(BigInteger768([
            0xace91a518173d79e,
            0x85b4b32d9e2c29f1,
            0xc23398b50ae36fb1,
            0xe4ad44193654931b,
            0x206960e627a3db5,
            0x3b5b892e388858c3,
            0x784e2a0e051839b6,
            0xeea4ed5c465934b3,
            0xe22149dfeadd89f3,
            0x31f8225a1ce56854,
            0x6cacc0691fc39968,
            0x4c451bf14739,
        ])),
    );
    a.mul_assign(&Fq3::new(
        Fq::from_repr(BigInteger768([
            0xea3d71dace0032be,
            0x9e32fec5ecf8847,
            0xdc741dcaaa904521,
            0x17cf97ee724d55f3,
            0xeb0c6f944ccdccf2,
            0xbe2d742884ee1ade,
            0xed72d008693e424d,
            0x29d565346140503b,
            0xae6d1807f4ddb996,
            0xc2c87df1290f12b3,
            0x25c9e97ea9a6168c,
            0x14bcbd0c05740,
        ])),
        Fq::from_repr(BigInteger768([
            0xa7c248bf356ea80f,
            0xc9915a5e4713b051,
            0x96bfd4bedee3d075,
            0x86774975113b713,
            0x5cb302606c76bce7,
            0x20eb52e0ec2b5a90,
            0xce74054a1d7cda7d,
            0xcd5114716dd466fc,
            0xd92b87655b255f7f,
            0xa612d7a0a2cdfa50,
            0xdcbfdf74619fd088,
            0xe83492a23080,
        ])),
        Fq::from_repr(BigInteger768([
            0x6878848d258bd4e5,
            0x21957e8c41feefb9,
            0x3fcf69a08f1e9410,
            0xe1e3056e0fa30e3c,
            0x5808b8487bbe7fc4,
            0xb5cfe60844f7d6a2,
            0x6f706a5d16a5aa41,
            0x7d438e19d0a302e9,
            0xacccbb506f192a14,
            0xd7683e432f46ad80,
            0xbd42d09e48f32be3,
            0xd000d94a201b,
        ])),
    ));
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0x546b6559354ad33a,
                0x7cb3fc50b5a89374,
                0xdbeeecf57666280a,
                0xc2a90279a880bb6f,
                0x93d97a513df51f25,
                0x5c64fa3416f7f5fd,
                0x5eb4cf7680c263b7,
                0x536d8cf5eb09798b,
                0x68329591d9896a9b,
                0x3e692f0eadcb46a1,
                0xb87a30e08cbcba96,
                0x12819a4a1b8b1,
            ])),
            Fq::from_repr(BigInteger768([
                0x5dd5aa31a378481,
                0x2c308bb900870e56,
                0x63f20733e3579087,
                0x47153c90a1180256,
                0x5b7d3967f14a60a8,
                0x4995effcc846e2ea,
                0x1e854b6391135556,
                0x14aa2a716b4d9dbc,
                0xc0ab4f9869284465,
                0x6b661a8addf630a0,
                0x48ebf9a3fb386998,
                0x15cb12c006579,
            ])),
            Fq::from_repr(BigInteger768([
                0x43052d439c9f6567,
                0x481cdd48b1ed431c,
                0x790e3dddabd89525,
                0x69764daa04dd82ea,
                0x68f8eddae35c7197,
                0x889d19480480d2ba,
                0xfec70cb30b86c860,
                0xb8c82f5b93784ec6,
                0x2badf2fecc0092e2,
                0x4cf2853fe24e12d5,
                0xa1027b754aa3b1b7,
                0x1bead97aeff44,
            ])),
        )
    );
}

#[test]
fn test_fq3_inverse() {
    assert!(Fq3::zero().inverse().is_none());

    let a = Fq3::new(
        Fq::from_repr(BigInteger768([
            0xea2f5a5b285605a9,
            0x9b9811086eb9a36c,
            0xeb4155b27eb0fd70,
            0xd78395c5e2265dfb,
            0x38c7271a32d8d316,
            0x16d69e418c964f89,
            0x391e37e6cb9f7517,
            0x14a5d0cff50ee3de,
            0xbdbf5f1b2c60e32a,
            0xe86653ca83ceacb9,
            0xcb0c8450b3a1d50c,
            0x9eea8f1bd345,
        ])),
        Fq::from_repr(BigInteger768([
            0xe96321d6846bd938,
            0xc44457d924ebf799,
            0x64a87c9840b9f92f,
            0x7ec7476a0f7ac1d0,
            0xa7358171b0e872e0,
            0x46454b2e1fda0a3c,
            0x9f9d75e878378870,
            0xf1c0de78d4ced8e5,
            0x14c8681949bdbf8e,
            0x890ec974d4890f7a,
            0x8b6154ccda8486d2,
            0x54817d11357b,
        ])),
        Fq::from_repr(BigInteger768([
            0xe0b69c3f609e009e,
            0x32cd8f0daadaee,
            0x4098e710b4abdeb,
            0xa70f14273d01899d,
            0x7a72bf882a1b3019,
            0x1bcae3482320a23c,
            0xf0830d313f84efe7,
            0x17e00a2bedc5d2c2,
            0x912329932f17d4fe,
            0x8ab61f747de901a2,
            0x8b5728c254c913a1,
            0x12478efb39124,
        ])),
    );
    let a = a.inverse().unwrap();
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0xf6b123e3017ea59c,
                0xaae1b9ecd171bf54,
                0xed292157532e6b9,
                0xa29971bef1c7fd64,
                0x2754fddf487b12e,
                0xc2071ba19a8078d1,
                0x1f70013a7b8b523f,
                0xa8eabba8d7e3e06,
                0x51c8b51fdd68c5c0,
                0x467942fdc4d667db,
                0x6a26b87a2e55958a,
                0x63ea13b8d65f,
            ])),
            Fq::from_repr(BigInteger768([
                0x4f3c45fa94981f88,
                0xab180727571e1d30,
                0x400c3c0c47b831be,
                0x7f1be9868780fef9,
                0xe644700a947c200a,
                0x4b18c27e6ae5358f,
                0x3347724766e20ca4,
                0xd4a6c139c5a8d60,
                0x112e45da5f8e4b2d,
                0x5139fc7a469c4cb5,
                0x6a1f4a3d81d00c0,
                0x9385286edaea,
            ])),
            Fq::from_repr(BigInteger768([
                0xcc69730b69a27ce8,
                0xb20d449bf5cc5206,
                0xfcc8e81422a637d6,
                0xb8aa02f107456a00,
                0x98a2ee22af2c87bc,
                0x673097c76eb3a89b,
                0x90f8a8331c814b14,
                0x20735cbe81a9f60d,
                0x4f7fa65e8eff3455,
                0x2ef59d3d073df835,
                0xd4ed3a61398b57d9,
                0x23924edf3ce6,
            ])),
        )
    );
}

#[test]
fn test_fq3_addition() {
    let mut a = Fq3::new(
        Fq::from_repr(BigInteger768([
            0x70a9f6f7d7cfc496,
            0x84e1e992071a39a3,
            0xcce6af0aa9d42cd4,
            0x16da3ccd8aa3bce7,
            0x3abf78b7a0ccdedb,
            0x64533c073dc56616,
            0x57d48514426d003a,
            0x76c57931a08a7b94,
            0x408797ae7e5bf40c,
            0xb5a8219fe218369a,
            0x4871ad4c341c20dc,
            0x1b2d12ecc2d26,
        ])),
        Fq::from_repr(BigInteger768([
            0xdde2779130322cd8,
            0x3a5c2c1d565026a,
            0xa1cfca0f4ac1f2e7,
            0xc713963da24dacca,
            0xd6f7d36691cdd4a6,
            0xc67ebdcc0748810d,
            0x5d116e4ce8bdec44,
            0x5a1be86ecedb6931,
            0xbbb7b3a08620dca9,
            0xab47ba8d382bbc47,
            0x8276191222a801ea,
            0x167415192e699,
        ])),
        Fq::from_repr(BigInteger768([
            0x43311da439d6bda8,
            0x4a4c249c61ad4da5,
            0xfa2add6d7da7e98b,
            0x22735466fa682a2,
            0x9a8497a2aac7c4d3,
            0xf470d9c28a935608,
            0x30201634f7b3fee6,
            0xae48ae0d58656f61,
            0x45fdf3d91ae8a9e0,
            0xf4e0a80511532a6d,
            0x5e4da4b94811844d,
            0xf5315ce200d,
        ])),
    );
    a.add_assign(&Fq3::new(
        Fq::from_repr(BigInteger768([
            0x8a5b3caa7f25f292,
            0x319ab2fa79515ecb,
            0x422065fe8ece17ab,
            0xc85e6d69508cd1d9,
            0x57ff297e60fe88fb,
            0x317a6603ff60fb38,
            0x887f460c565c190,
            0x846e7d862cac4884,
            0xc5ef1c4f1a978a07,
            0xc1ae30a0818ef210,
            0x5aac8601312cac55,
            0xc382c8f64875,
        ])),
        Fq::from_repr(BigInteger768([
            0x40a20c619e79826c,
            0x69a3132c23d60833,
            0x96c322803a1669fb,
            0x96a26754fb90c64d,
            0x25bafad53d14fbdb,
            0xd23fcfe959d5d73c,
            0x61f2474a938bd25f,
            0xac66b0233da24bf5,
            0x5b8b6bf66139b5ec,
            0x592cb115aa3292ad,
            0x98a061fbe7a35eda,
            0x19059b30bba93,
        ])),
        Fq::from_repr(BigInteger768([
            0xa0450c373073d290,
            0x83566a6c362b7610,
            0xc05579d6cbc72938,
            0x20f59bc313cf6f25,
            0x7feb5a333872c8c4,
            0x3dfd55628b2d3f9f,
            0x435dd011644a05e7,
            0xf7059c8a97467555,
            0xcd8fe71709865ebc,
            0x248f8286f84fc456,
            0x4dce63167743ac32,
            0xff16fdcc042,
        ])),
    ));
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0x21fdbcc016f5b727,
                0x67dc037570ca5e1f,
                0x3843934cf9a1ece8,
                0x2558b0c0a697541c,
                0x53d2d7a1d8aa41a0,
                0xe36145e274ccb7b3,
                0xc68b549b6677ca2c,
                0xf3363d91e495d68a,
                0xa7becb042c5ba5a0,
                0xbf5cbaf0081778bd,
                0x92fba32a7665ff84,
                0xb18dca2fb18a,
            ])),
            Fq::from_repr(BigInteger768([
                0x457d0d108eabaf43,
                0x1ea83cd6e999d04e,
                0x61cf6ad345d8054b,
                0xa3d6041c69453873,
                0xbdc703a7a5c1aa4b,
                0xe652318c98c4aeae,
                0x253290bddaeec706,
                0xfe84df6c23dcc799,
                0xb88b369d7ac2ba22,
                0x4c7ad45286ce9f07,
                0xaf3eaeb1b689317,
                0x132d4d70bdd1c,
            ])),
            Fq::from_repr(BigInteger768([
                0xe37629db6a4a9038,
                0xcda28f0897d8c3b5,
                0xba805744496f12c3,
                0x231cd1098375f1c8,
                0x1a6ff1d5e33a8d97,
                0x326e2f2515c095a8,
                0x737de6465bfe04ce,
                0xa54e4a97efabe4b6,
                0x138ddaf0246f089d,
                0x19702a8c09a2eec4,
                0xac1c07cfbf553080,
                0x1f4485aae04f,
            ])),
        )
    );
}

#[test]
fn test_fq3_subtraction() {
    let mut a = Fq3::new(
        Fq::from_repr(BigInteger768([
            0x466f53ee62857c35,
            0xc2db6a3ac5d73b17,
            0x4e4b3b2667890ce,
            0xc439a945528dfb6d,
            0x658b45180b928613,
            0x6ba00c3473316448,
            0xe69d313374ac8a1a,
            0x74909447e6beb138,
            0x35b02d0d917d55a2,
            0x7a262f207791f762,
            0xc3f96812fa360cc4,
            0xc589fccf404a,
        ])),
        Fq::from_repr(BigInteger768([
            0x62079cdb754f6283,
            0xac492f5b26fcd792,
            0xa90a15721c31e5c3,
            0x2aa1a00708752976,
            0x25ffd74d046f9865,
            0xacbea3e3977ed71a,
            0x4c9dd2c1740d23ff,
            0x5439242985eedace,
            0x20015c9662207818,
            0x7efc2900c8dcd296,
            0x79b705f651b64b37,
            0x241a6f488ef3,
        ])),
        Fq::from_repr(BigInteger768([
            0x7421f4deec25cdf,
            0xe333e3b3f7ad9635,
            0xf4e0e90a81ecdfa4,
            0x30b664c73744090d,
            0xa2c939517669ca3,
            0x26c843d7e091fd84,
            0x12e2b5d09078573c,
            0xd65c1fcba1ce8535,
            0x53ad9c33664134d5,
            0x510d8bb445fb7364,
            0x65d97212aad39715,
            0x13e02b8c098aa,
        ])),
    );
    a.sub_assign(&Fq3::new(
        Fq::from_repr(BigInteger768([
            0xae66135f4c0de452,
            0xec320d88395e5e0,
            0x8cba68af85c9ed0f,
            0x996d5b81871aae5a,
            0x8eefe077a5f8904b,
            0x7063aa08173520ec,
            0x8e724a8c7a2cdb52,
            0x2220f0d274d9280d,
            0xe7faf2c8b08a53d6,
            0x3d7adea7f6f2feba,
            0x7388689327cc879e,
            0x12204bc594c27,
        ])),
        Fq::from_repr(BigInteger768([
            0xc0f4c0c54cff46e4,
            0xbd5d05139c379672,
            0x42e975e7d48f0bf8,
            0x854b54e796fceead,
            0xf70666bd557b9d3a,
            0x4761fd7d60b15402,
            0x58b6c761b51c5c38,
            0x297198dfad1e9ad0,
            0xd98876ab040eaa20,
            0x49d6ff00de45a75c,
            0x783f56ea7eeab38a,
            0xa75818e565b,
        ])),
        Fq::from_repr(BigInteger768([
            0x357041f1e563353a,
            0x4d7884c47ebd50d6,
            0x571fa308fb790e92,
            0x132a83fe50aa0d34,
            0x8a5615b190de6c96,
            0xf540a5d8a1c27249,
            0x50b718961ba058fc,
            0xe3f754eb96ffff41,
            0xc936c2b501095e6f,
            0x2b81b529966a348c,
            0x583520d69bc45417,
            0xeda15eca386f,
        ])),
    ));
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0x7110b771567797e4,
                0x2b8e27951e28f86,
                0x4eedccbf1faefb57,
                0xe4ac473a000c87b7,
                0x15872f348ebb1bfe,
                0xada8be552455ecf7,
                0xf1fc0b809bdaa665,
                0x5a6d5c9b5a8676b8,
                0xac6d233e4d8ada3f,
                0xf4a4e7c8dc2ea894,
                0x60938fa2c14c52d3,
                0x1684b6e08b834,
            ])),
            Fq::from_repr(BigInteger768([
                0xa112dc1628501b9f,
                0xeeec2a478ac5411f,
                0x66209f8a47a2d9ca,
                0xa5564b1f71783ac9,
                0x2ef9708faef3fb2a,
                0x655ca66636cd8317,
                0xf3e70b5fbef0c7c7,
                0x2ac78b49d8d03ffd,
                0x4678e5eb5e11cdf8,
                0x352529ffea972b39,
                0x177af0bd2cb97ad,
                0x19a4edba3898,
            ])),
            Fq::from_repr(BigInteger768([
                0xd1d1dd5c095f27a5,
                0x95bb5eef78f0455e,
                0x9dc146018673d112,
                0x1d8be0c8e699fbd9,
                0x7fd67de38688300d,
                0x31879dff3ecf8b3a,
                0xc22b9d3a74d7fe3f,
                0xf264cae00ace85f3,
                0x8a76d97e6537d665,
                0x258bd68aaf913ed7,
                0xda4513c0f0f42fe,
                0x506159f6603b,
            ])),
        )
    );
}

#[test]
fn test_fq3_negation() {
    let mut a = Fq3::new(
        Fq::from_repr(BigInteger768([
            0x45a4987bd3e0d56d,
            0xa001d1f786af5915,
            0x3b4e2747265dd474,
            0x6e4f1694c5e2c3c5,
            0xe2a2e49750d269,
            0x7ad3a81da9108d86,
            0x26b38ef00c89a2fb,
            0x3f549ebf55695b14,
            0x30556837056daf5e,
            0x90d6790c0b3928d0,
            0x6d91bbc81e4db043,
            0xa6357fe3f93d,
        ])),
        Fq::from_repr(BigInteger768([
            0x431babec51a9126f,
            0x7010cbe98cb9317e,
            0x68441948e9dd2030,
            0x5062859a46134908,
            0xf780c38412a2eca7,
            0xeccf583f8b368e3c,
            0x7f41ac709c0f00fc,
            0x1001838301b18f15,
            0x7366d796a5c2ac58,
            0x96df478644ed2323,
            0xe70fb635feae00bd,
            0x10d776a6661d,
        ])),
        Fq::from_repr(BigInteger768([
            0x876ea7e4a25c62f1,
            0x630ce5add7dcc4ee,
            0x411fa998cf05f747,
            0x7f777e93c2bb82a9,
            0xd45a316491697a3c,
            0x73d4c32b5d32cd3d,
            0x4d8837c1b004e6ed,
            0x697249a5886df985,
            0x99256513a6c0a93f,
            0x78295b06f24b2acc,
            0x674b3104ee4258bf,
            0x8439c1521467,
        ])),
    );
    a = -a;
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0x9362de666c1f2a94,
                0xae9ec71f88f1e13a,
                0x9b755a7518a28322,
                0x4b90e2e16eb676df,
                0x3e0927af91d053cd,
                0x3798b40b1f491c15,
                0x731d95e994d154a2,
                0xc8a91a6693379279,
                0x2e6280c2672a2914,
                0x27231e445056871d,
                0xa290d45ad0951d6a,
                0x11e90adaecad3,
            ])),
            Fq::from_repr(BigInteger768([
                0x95ebcaf5ee56ed92,
                0xde8fcd2d82e808d1,
                0x6e7f687355233766,
                0x697d73dbee85f19c,
                0x476b0710167e398f,
                0xc59d03e93d231b5e,
                0x1a8f7869054bf6a0,
                0xf7fc35a2e6ef5e78,
                0xeb511162c6d52c1a,
                0x211a4fca16a28cc9,
                0x2912d9ecf034ccf0,
                0x1b3eeb6ec5df3,
            ])),
            Fq::from_repr(BigInteger768([
                0x5198cefd9da39d10,
                0xeb93b36937c47561,
                0x95a3d8236ffa604f,
                0x3a687ae271ddb7fb,
                0x6a91992f97b7abfa,
                0x3e9798fd6b26dc5d,
                0x4c48ed17f15610b0,
                0x9e8b6f806032f408,
                0xc59283e5c5d72f33,
                0x3fd03c4969448520,
                0xa8d75f1e00a074ee,
                0x1408c6c40afa9,
            ])),
        )
    );
}

#[test]
fn test_fq3_doubling() {
    let mut a = Fq3::new(
        Fq::from_repr(BigInteger768([
            0xd192ae0d65f172e5,
            0xa973e6d47e71a73e,
            0x794ef7d80e97a231,
            0x388d30fad5bce669,
            0x828c216522e230a0,
            0xb1b9903e98e589a3,
            0x3bab2d4dfc3b4d98,
            0x511bb8398504956,
            0xd995de739cd417c3,
            0x5da89657b5dfde2f,
            0x932e48738bf2ead7,
            0xd35bfbd24ffa,
        ])),
        Fq::from_repr(BigInteger768([
            0x5c83fc6f4c683523,
            0xb37a157d260571fd,
            0xa4e5ba9d2c827a0a,
            0x81fd741407dd59,
            0x40a916dedb48c79a,
            0x9094f99a14b5ddd5,
            0xa86f628c0e0f313f,
            0xf0159069ad16e6cc,
            0x1739f0de0a7be25,
            0x7921849db4a258fb,
            0x3188837a5d4dda25,
            0x58d43cab04bc,
        ])),
        Fq::from_repr(BigInteger768([
            0x6160e374a4eb163a,
            0xadf9a96ea1729143,
            0xbd0c2b5aa0006c03,
            0xc20fa3537d5b2672,
            0x1b230f48f7f4cbf4,
            0xc4b47d6737ef05f4,
            0xdf8bd07bb9d562e0,
            0x893163812adfac13,
            0x13f4bdfb094a057b,
            0xd7ad1d85b9953c25,
            0xdf524e2d7f802a5d,
            0x79893fa4c427,
        ])),
    );
    a.double_in_place();
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0xa3255c1acbe2e5ca,
                0x52e7cda8fce34e7d,
                0xf29defb01d2f4463,
                0x711a61f5ab79ccd2,
                0x51842ca45c46140,
                0x6373207d31cb1347,
                0x77565a9bf8769b31,
                0xa23770730a092ac,
                0xb32bbce739a82f86,
                0xbb512caf6bbfbc5f,
                0x265c90e717e5d5ae,
                0x1a6b7f7a49ff5,
            ])),
            Fq::from_repr(BigInteger768([
                0xb907f8de98d06a46,
                0x66f42afa4c0ae3fa,
                0x49cb753a5904f415,
                0x103fae8280fbab3,
                0x81522dbdb6918f34,
                0x2129f334296bbbaa,
                0x50dec5181c1e627f,
                0xe02b20d35a2dcd99,
                0x2e73e1bc14f7c4b,
                0xf243093b6944b1f6,
                0x631106f4ba9bb44a,
                0xb1a879560978,
            ])),
            Fq::from_repr(BigInteger768([
                0xc2c1c6e949d62c74,
                0x5bf352dd42e52286,
                0x7a1856b54000d807,
                0x841f46a6fab64ce5,
                0x36461e91efe997e9,
                0x8968face6fde0be8,
                0xbf17a0f773aac5c1,
                0x1262c70255bf5827,
                0x27e97bf612940af7,
                0xaf5a3b0b732a784a,
                0xbea49c5aff0054bb,
                0xf3127f49884f,
            ])),
        )
    );
}

#[test]
fn test_fq3_frobenius_map() {
    let mut a = Fq3::new(
        Fq::from_repr(BigInteger768([
            0x5ae08dbcf05a7e94,
            0x7a3ae54580a38126,
            0xa44b6b87c98acd54,
            0x12fbbf7aac1c4257,
            0x83574728a5dd4cd,
            0x94f4abd8fb7572aa,
            0xb0cc279b783eb119,
            0xf7e82969a91de53a,
            0xcdf097798212fc90,
            0xb3a1873b082f9c90,
            0x7fbd55a49371f883,
            0x1281fb018fed1,
        ])),
        Fq::from_repr(BigInteger768([
            0x6fc7aff59a9b9fca,
            0x38d88519ba8e4994,
            0x3b0ce22e5d10e284,
            0xec07648b1d926326,
            0x4ad75a653684893,
            0x8787f090de330059,
            0x4ba554a712f7017c,
            0x48f9347afc8eafc1,
            0x1ff13912b5782b6e,
            0x7878e03cb713f699,
            0xaf28937e1c342778,
            0x12d05bdbb114d,
        ])),
        Fq::from_repr(BigInteger768([
            0x26738434e6f99141,
            0x21733ea49be3f6cb,
            0xcefb33846962f4f0,
            0x9ab4e36cdb571399,
            0x2b74ea5a7e0217dc,
            0x3adb3a4f7f7084e1,
            0x89a0611519e1952f,
            0xb69a42ce200b8f7,
            0x580f9df838107915,
            0xda957fb3b319f8af,
            0x73ce2e50e6d0fa3,
            0x8a5add60133,
        ])),
    );
    a.frobenius_map(0);
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0x5ae08dbcf05a7e94,
                0x7a3ae54580a38126,
                0xa44b6b87c98acd54,
                0x12fbbf7aac1c4257,
                0x83574728a5dd4cd,
                0x94f4abd8fb7572aa,
                0xb0cc279b783eb119,
                0xf7e82969a91de53a,
                0xcdf097798212fc90,
                0xb3a1873b082f9c90,
                0x7fbd55a49371f883,
                0x1281fb018fed1,
            ])),
            Fq::from_repr(BigInteger768([
                0x6fc7aff59a9b9fca,
                0x38d88519ba8e4994,
                0x3b0ce22e5d10e284,
                0xec07648b1d926326,
                0x4ad75a653684893,
                0x8787f090de330059,
                0x4ba554a712f7017c,
                0x48f9347afc8eafc1,
                0x1ff13912b5782b6e,
                0x7878e03cb713f699,
                0xaf28937e1c342778,
                0x12d05bdbb114d,
            ])),
            Fq::from_repr(BigInteger768([
                0x26738434e6f99141,
                0x21733ea49be3f6cb,
                0xcefb33846962f4f0,
                0x9ab4e36cdb571399,
                0x2b74ea5a7e0217dc,
                0x3adb3a4f7f7084e1,
                0x89a0611519e1952f,
                0xb69a42ce200b8f7,
                0x580f9df838107915,
                0xda957fb3b319f8af,
                0x73ce2e50e6d0fa3,
                0x8a5add60133,
            ])),
        )
    );
    a.frobenius_map(1);
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0x5ae08dbcf05a7e94,
                0x7a3ae54580a38126,
                0xa44b6b87c98acd54,
                0x12fbbf7aac1c4257,
                0x83574728a5dd4cd,
                0x94f4abd8fb7572aa,
                0xb0cc279b783eb119,
                0xf7e82969a91de53a,
                0xcdf097798212fc90,
                0xb3a1873b082f9c90,
                0x7fbd55a49371f883,
                0x1281fb018fed1,
            ])),
            Fq::from_repr(BigInteger768([
                0x3110a46cb0ed77f1,
                0xe1e033e34bd5fa34,
                0xc5f2da0bf19faee6,
                0xaa4b967c01a0b0a5,
                0xa491f6d1abe3bea5,
                0x24c45b5a04c13a77,
                0x8ec36da8681eaa21,
                0x15e5cc898998c2,
                0x83e9fa5d6864f771,
                0xbd110fb01b05d02,
                0x1e3d1d4ce9b15d5c,
                0x4c25d114bc14,
            ])),
            Fq::from_repr(BigInteger768([
                0x57f33e76b7aa6e5e,
                0x3296ee1b127f5a4d,
                0x4223426a6bccfa89,
                0xd36f0b46c0347d4d,
                0xd746bc9a42d72ff6,
                0x12f2f98837410346,
                0x68600cbaadcdae61,
                0xf04297331b1b5c9f,
                0xf4b1a2b62b5218d1,
                0x90b1594c69063bc,
                0xd2d70a164f85328b,
                0x187c61d058bf8,
            ])),
        )
    );
    a.frobenius_map(1);
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0x5ae08dbcf05a7e94,
                0x7a3ae54580a38126,
                0xa44b6b87c98acd54,
                0x12fbbf7aac1c4257,
                0x83574728a5dd4cd,
                0x94f4abd8fb7572aa,
                0xb0cc279b783eb119,
                0xf7e82969a91de53a,
                0xcdf097798212fc90,
                0xb3a1873b082f9c90,
                0x7fbd55a49371f883,
                0x1281fb018fed1,
            ])),
            Fq::from_repr(BigInteger768([
                0x382f227ff476e846,
                0x33e7e01a093cf687,
                0xd5c3c581f04fc62c,
                0x238cfe6f156626d8,
                0x95ac5e1c29d51efd,
                0x620103de5656eca,
                0xbf68628a26454c00,
                0xbeee9ede6288a509,
                0xbadcb5894ebab593,
                0x33afa618a2cb5c51,
                0x42bcdf57e8fd48d9,
                0x4b9a9ec2f6af,
            ])),
            Fq::from_repr(BigInteger768([
                0x5aa0b436a15c0062,
                0xfa966c57613de937,
                0xc5a50bcd69d0681d,
                0x4bbc0ac2990da9bd,
                0x3c30239f6847de63,
                0x649e285111a82173,
                0xa7d0b709d9abb40d,
                0xc517dc5eb84d7f6,
                0x11f6a84b0935468c,
                0xd4590207e1e55381,
                0x360ea32790f08b7e,
                0x345a62b736e5,
            ])),
        )
    );
    a.frobenius_map(2);
    assert_eq!(
        a,
        Fq3::new(
            Fq::from_repr(BigInteger768([
                0x5ae08dbcf05a7e94,
                0x7a3ae54580a38126,
                0xa44b6b87c98acd54,
                0x12fbbf7aac1c4257,
                0x83574728a5dd4cd,
                0x94f4abd8fb7572aa,
                0xb0cc279b783eb119,
                0xf7e82969a91de53a,
                0xcdf097798212fc90,
                0xb3a1873b082f9c90,
                0x7fbd55a49371f883,
                0x1281fb018fed1,
            ])),
            Fq::from_repr(BigInteger768([
                0x3110a46cb0ed77f1,
                0xe1e033e34bd5fa34,
                0xc5f2da0bf19faee6,
                0xaa4b967c01a0b0a5,
                0xa491f6d1abe3bea5,
                0x24c45b5a04c13a77,
                0x8ec36da8681eaa21,
                0x15e5cc898998c2,
                0x83e9fa5d6864f771,
                0xbd110fb01b05d02,
                0x1e3d1d4ce9b15d5c,
                0x4c25d114bc14,
            ])),
            Fq::from_repr(BigInteger768([
                0x57f33e76b7aa6e5e,
                0x3296ee1b127f5a4d,
                0x4223426a6bccfa89,
                0xd36f0b46c0347d4d,
                0xd746bc9a42d72ff6,
                0x12f2f98837410346,
                0x68600cbaadcdae61,
                0xf04297331b1b5c9f,
                0xf4b1a2b62b5218d1,
                0x90b1594c69063bc,
                0xd2d70a164f85328b,
                0x187c61d058bf8,
            ])),
        )
    );
}

#[test]
fn test_fq3_legendre() {
    use crate::fields::LegendreSymbol::*;

    assert_eq!(Zero, Fq3::zero().legendre());
    assert_eq!(QuadraticNonResidue, Fq3Parameters::NONRESIDUE.legendre());
    // i^2 = -1
    let mut m1 = -Fq3::one();
    assert_eq!(QuadraticResidue, m1.legendre());
    m1 = Fq6Parameters::mul_fp3_by_nonresidue(&m1);
    assert_eq!(QuadraticNonResidue, m1.legendre());

    assert_eq!(QuadraticNonResidue, Fq6Parameters::NONRESIDUE.legendre());
}

#[test]
fn test_fq3_mul_nonresidue() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let nqr = Fq3::new(Fq::zero(), Fq::one(), Fq::zero());

    for _ in 0..1000 {
        let mut a = Fq3::rand(&mut rng);
        let mut b = a;
        a = Fq6Parameters::mul_fp3_by_nonresidue(&a);
        b.mul_assign(&nqr);

        assert_eq!(a, b);
    }
}

#[test]
fn test_fq6_mul_by_2345() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        let c2 = Fq::rand(&mut rng);
        let c3 = Fq::rand(&mut rng);
        let c4 = Fq::rand(&mut rng);
        let c5 = Fq::rand(&mut rng);
        let to_mul = Fq6::new(Fq3::new(Fq::zero(), Fq::zero(), c2), Fq3::new(c3, c4, c5));
        let a = Fq6::rand(&mut rng);
        let mut b = a;

        b.mul_assign(&to_mul);

        assert_eq!(a.mul_by_2345(&to_mul), b);
    }
}
