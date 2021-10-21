use crate::{
    biginteger::{BigInteger, BigInteger768},
    bytes::{FromBytes, ToBytes},
    fields::mnt4753::{Fq, Fq2, Fq2Parameters, Fq4, Fq4Parameters, FqParameters, Fr},
    fields::models::{fp2::Fp2Parameters, fp4::Fp4Parameters},
    fields::tests::{field_test, frobenius_test, primefield_test, sqrt_field_test},
    fields::FpParameters,
    to_bytes, Field, PrimeField, SemanticallyValid, SquareRootField, ToBits, UniformRand,
};

use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use std::{
    cmp::Ordering,
    ops::{AddAssign, MulAssign, SubAssign},
};

pub(crate) const ITERATIONS: usize = 5;

#[test]
fn test_mnt4753_fr() {
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
fn test_mnt4753_fq() {
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
fn test_mnt4753_fq2() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let a: Fq2 = UniformRand::rand(&mut rng);
        let b: Fq2 = UniformRand::rand(&mut rng);
        field_test(a, b);
        sqrt_field_test(a);
    }
    frobenius_test::<Fq2, _>(Fq::characteristic(), 13);
}

#[test]
fn test_mnt4753_fq4() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    for _ in 0..ITERATIONS {
        let g: Fq4 = UniformRand::rand(&mut rng);
        let h: Fq4 = UniformRand::rand(&mut rng);
        field_test(g, h);
    }
    frobenius_test::<Fq4, _>(Fq::characteristic(), 13);
}

#[test]
fn test_frob_coeffs() {
    //Fq2 coefficients test
    let nqr = Fq::new(BigInteger768([
        11881297496860141143,
        13588356353764843511,
        9969398190777826186,
        17325157081734070311,
        16341533986183788031,
        8322434028726676858,
        13631157743146294957,
        8365783422740577875,
        3010239015809771096,
        11776256826687733591,
        7214251687253691272,
        268626707558702,
    ]));

    assert_eq!(Fq2Parameters::FROBENIUS_COEFF_FP2_C1[0], Fq::one());
    assert_eq!(
        Fq2Parameters::FROBENIUS_COEFF_FP2_C1[1],
        nqr.pow([
            0xAF4831EF122F4000,
            0x71CEAA29166E88CF,
            0x31C40838CD6212F8,
            0x342D6674BB392A52,
            0xDC0786D2E5A9BF1C,
            0xD88BF3BB790C02CE,
            0xCCE8926CD0AD7BCE,
            0x83FEDC92F45076C6,
            0xAF5BF47CB64BEC39,
            0xDBFCCBA82DC7D7F6,
            0x88114811777166D6,
            0xE26316C96208,
        ])
    );

    //Fq4 coefficients test
    assert_eq!(Fq4Parameters::FROBENIUS_COEFF_FP4_C1[0], Fq::one());
    assert_eq!(
        Fq4Parameters::FROBENIUS_COEFF_FP4_C1[1],
        nqr.pow([
            0xd7a418f78917a000,
            0x38e755148b374467,
            0x18e2041c66b1097c,
            0x1a16b33a5d9c9529,
            0x6e03c36972d4df8e,
            0x6c45f9ddbc860167,
            0x667449366856bde7,
            0xc1ff6e497a283b63,
            0x57adfa3e5b25f61c,
            0x6dfe65d416e3ebfb,
            0x4408a408bbb8b36b,
            0x71318b64b104,
        ])
    );

    assert_eq!(
        Fq4Parameters::FROBENIUS_COEFF_FP4_C1[2],
        nqr.pow([
            0xbb4c5fa7a22f4000,
            0xaa577656adec411c,
            0x818ea97cade6ed97,
            0xe20199288067443e,
            0x8e9c8ed556dc1767,
            0xb137c8ec23877dd9,
            0xe894d15ac94b88c1,
            0x31c2bf97498dcd49,
            0x481449239e4ea77,
            0x56e583b5fd2b5720,
            0x7f82bdde428d7c41,
            0x4f26223cae3a7324,
            0x52975e00b52b86b5,
            0xfe534fe59436d1ab,
            0xbd730a4837ec1719,
            0x4e8b576b5d85308d,
            0x5841c9f6ec780a78,
            0x135ab6ca7e688bd3,
            0xf92b763730ad9922,
            0x2a03f90b89a196a7,
            0x147944fe3cc40aa7,
            0xcbd232a067a67771,
            0x32e01e1c4eefc062,
            0xc8331a96,
        ])
    );

    let t: Vec<u64> = vec![
        0x5fd9fc104b46e000,
        0x1e345e64bc3f81d4,
        0xb34098ec25558a87,
        0xb282f2aa887e70ba,
        0x1b87c5580fb5fa31,
        0xb9fbe059f4d517f5,
        0xbd7d5cda58619ae0,
        0x890a0f3fa910756e,
        0x6e333bb25edc0865,
        0xee2630cae3709e25,
        0xb97f44a556e8e77c,
        0xf77783b54e8a1b16,
        0xeb017c6e0e6a9414,
        0x1e168ed166486490,
        0x731eccb21637216c,
        0xe0db563bbcf86ea3,
        0xe7ca94bbe8ffd7bb,
        0xf3fce5d2b115dc57,
        0x80a5f12d11994856,
        0x9be3b3a9d59d0c2c,
        0x88e2ef22bc9aee8e,
        0x81ee17c36e6941dd,
        0x23b8ae1b66484554,
        0xd6d1fe64efb02ca5,
        0xb1ed04968b94ce52,
        0xad7a8ff0bffcf3dc,
        0x2bd711b36662ed31,
        0xb848a251c83b23f3,
        0x96f629d6d9d61eb2,
        0x8c55f5d4afea045b,
        0xcfe097a40a96ca99,
        0x240bd338c822d715,
        0xa7a3a44bc8779588,
        0x76c1bb0b23579a9e,
        0x361b84251779aa78,
        0x16215,
    ];
    assert_eq!(Fq4Parameters::FROBENIUS_COEFF_FP4_C1[3], nqr.pow(t));
}

#[test]
fn test_neg_one() {
    let neg_one = Fq::new(BigInteger768([
        0xc5e777324a8210bf,
        0x51d0228bd2d9cb18,
        0xcbc42bd0cdafcec2,
        0xef0234cfaee99ea2,
        0xcae87111aa4ae6c8,
        0x930899ec2314e834,
        0x67c4e9228e277244,
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
        0x20fc924b28d2f7d6,
        0xeee2288b24070b7f,
        0xbd14dcce936d92bf,
        0x7705edb97ebcd3f3,
        0x2c497d412bd2c3e8,
        0x9363f538ef90135d,
        0xb0109742cc4add3f,
        0x577389b8e8af372e,
        0xbb1fec3e1ab79a25,
        0xcc9c980eac0222e2,
        0xf738570ed0a42ffa,
        0x1c3b43f4ef84d,
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
        println!("{:?}", a);
        assert!(a.is_valid());
    }
}

#[test]
fn test_fq_add_assign() {
    {
        // Random number
        let mut tmp = Fq::new(BigInteger768([
            0x6b73b160e812b6a3,
            0x920f4bde3d1e6d70,
            0xad7c43767007402a,
            0xeec097f4edb6d94e,
            0xedd649af7ff8f8d9,
            0xd487b3c97fa1aff1,
            0x7ee5aa4cab2095c1,
            0x1682796bc0d18747,
            0xb13abeedc98acc1e,
            0x5df407ccca403f0c,
            0xef8fc932df51be4d,
            0x188f263e1b224,
        ]));
        assert!(tmp.is_valid());
        // Test that adding zero has no effect.
        tmp.add_assign(&Fq::new(BigInteger768::from(0)));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0x6b73b160e812b6a3,
                0x920f4bde3d1e6d70,
                0xad7c43767007402a,
                0xeec097f4edb6d94e,
                0xedd649af7ff8f8d9,
                0xd487b3c97fa1aff1,
                0x7ee5aa4cab2095c1,
                0x1682796bc0d18747,
                0xb13abeedc98acc1e,
                0x5df407ccca403f0c,
                0xef8fc932df51be4d,
                0x188f263e1b224,
            ]))
        );
        // Add one and test for the result.
        tmp.add_assign(&Fq::new(BigInteger768::from(1)));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0x6b73b160e812b6a4,
                0x920f4bde3d1e6d70,
                0xad7c43767007402a,
                0xeec097f4edb6d94e,
                0xedd649af7ff8f8d9,
                0xd487b3c97fa1aff1,
                0x7ee5aa4cab2095c1,
                0x1682796bc0d18747,
                0xb13abeedc98acc1e,
                0x5df407ccca403f0c,
                0xef8fc932df51be4d,
                0x188f263e1b224,
            ]))
        );
        // Add another random number that exercises the reduction.
        tmp.add_assign(&Fq::new(BigInteger768([
            0xe95e2d43caa35471,
            0x1bd18c806ebb4160,
            0xcde4889fb2596a9e,
            0x4a5e38f927c76670,
            0xb3ea6fcf75bd204,
            0x9a9206b28ec054e9,
            0xc8e1add798955218,
            0x70d33702049247f5,
            0xe9e48e1b44c79f63,
            0x7758f3a861973a17,
            0x97bc179c7b96dda7,
            0x3caf1bf9214d,
        ])));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0xf6417ac68e578b14,
                0xca43840c7efc9d31,
                0x17d8bba4879c84d7,
                0xd0c404049f0beb1b,
                0x4105e306ac014ca5,
                0xbe01d3051c49ff3d,
                0xadf6334aa25af03c,
                0x7f57f747dcc2e1af,
                0x3c67640fa1ba930e,
                0x1d536424d047c937,
                0x772950ac6c05ce47,
                0xdb52480f61,
            ]))
        );
        // Add one to (q - 1) and test for the result.
        tmp = Fq::new(BigInteger768([
            0x5E9063DE245E8000,
            0xE39D54522CDD119F,
            0x638810719AC425F0,
            0x685ACCE9767254A4,
            0xB80F0DA5CB537E38,
            0xB117E776F218059D,
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
            0x4655d097df1b68a1,
            0x4b03ad48db48a072,
            0x7922ccd99130dd8a,
            0x829a090815c0b81c,
            0xddffc641db15e8e1,
            0xdf8850a0111c699b,
            0xef1e3e7153698182,
            0xe89676fbf981d7c7,
            0x768dcc47e764ac0f,
            0x9666e3be5241e70,
            0x438440b524e356d4,
            0x11762e7161ca2,
        ]));
        tmp.add_assign(&Fq::new(BigInteger768([
            0x183a93464543175f,
            0x9899a7095194712d,
            0xea65439809934866,
            0xe5c0c3e160b19c87,
            0xda0f4763f03d9556,
            0xd18f96d6e0fb9c01,
            0xaab2e6684df1761a,
            0x1f674229ef1f15c5,
            0xe82a1cb185332c63,
            0xae932914766b917c,
            0xcc9e4f6dc9ff76d9,
            0xad63467ca76e,
        ])));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0x5E9063DE245E8000,
                0xE39D54522CDD119F,
                0x638810719AC425F0,
                0x685ACCE9767254A4,
                0xB80F0DA5CB537E38,
                0xB117E776F218059D,
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
            0xc7564a475de9a839,
            0x28ddfb67731c1d10,
            0x3969ef83e0a5be75,
            0x6c8f5a5217f1a34,
            0xa9ab64a80f4d0044,
            0x125885735ebeab75,
            0xbb3f84d5ab0025a8,
            0xeddd57a2d71ba40a,
            0x70c8164857c17c3c,
            0xe85fd4c0ef646fda,
            0xe5527d7bc0d7e068,
            0x77786c31b9fe,
        ]));
        tmp.sub_assign(&Fq::new(BigInteger768([
            0x7cdb8adc095c68e4,
            0xef7e55bfa2c7d6fe,
            0x612b4c9e017b6ab2,
            0x1bd77344474e7180,
            0x7717e7ff0fe957fb,
            0x68db8609c8a3733f,
            0x5f85ed476f99b66c,
            0x13a1c595a7d28470,
            0x9e659a162f4b77ab,
            0x38e5f1a907dd8c03,
            0x2be8ddc3c2b9f5e4,
            0x5b292b78ef11,
        ])));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0x4a7abf6b548d3f55,
                0x395fa5a7d0544612,
                0xd83ea2e5df2a53c2,
                0xeaf18260da30a8b3,
                0x32937ca8ff63a848,
                0xa97cff69961b3836,
                0x5bb9978e3b666f3b,
                0xda3b920d2f491f9a,
                0xd2627c3228760491,
                0xaf79e317e786e3d6,
                0xb9699fb7fe1dea84,
                0x1c4f40b8caed,
            ]))
        );

        // Test the opposite subtraction which doesn't test reduction.
        tmp = Fq::new(BigInteger768([
            0x7cdb8adc095c68e4,
            0xef7e55bfa2c7d6fe,
            0x612b4c9e017b6ab2,
            0x1bd77344474e7180,
            0x7717e7ff0fe957fb,
            0x68db8609c8a3733f,
            0x5f85ed476f99b66c,
            0x13a1c595a7d28470,
            0x9e659a162f4b77ab,
            0x38e5f1a907dd8c03,
            0x2be8ddc3c2b9f5e4,
            0x5b292b78ef11,
        ]));
        tmp.sub_assign(&Fq::new(BigInteger768([
            0xc7564a475de9a839,
            0x28ddfb67731c1d10,
            0x3969ef83e0a5be75,
            0x6c8f5a5217f1a34,
            0xa9ab64a80f4d0044,
            0x125885735ebeab75,
            0xbb3f84d5ab0025a8,
            0xeddd57a2d71ba40a,
            0x70c8164857c17c3c,
            0xe85fd4c0ef646fda,
            0xe5527d7bc0d7e068,
            0x77786c31b9fe,
        ])));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0x1415a472cfd140ac,
                0xaa3daeaa5c88cb8d,
                0x8b496d8bbb99d22e,
                0x7d694a889c41abf0,
                0x857b90fccbefd5ef,
                0x79ae80d5bfccd67,
                0x3e178d4b65f48862,
                0x2dc22718b957cdf3,
                0x8c556cc74421d3e1,
                0x87fb4387408cc16,
                0x56b8f06af0c4e329,
                0x1a876ecd9f923,
            ]))
        );

        // Test for sensible results with zero
        tmp = Fq::new(BigInteger768::from(0));
        tmp.sub_assign(&Fq::new(BigInteger768::from(0)));
        assert!(tmp.is_zero());

        tmp = Fq::new(BigInteger768([
            0xb2b74f8c784b6de,
            0x1d5ae51c4d7a2afc,
            0xf424485fa9ab9789,
            0x39948fcd5f1ba445,
            0xa8635673ec8e6c81,
            0x19b75b83a122ee72,
            0x31e7f890d12e9443,
            0x3dbca4d50c6a02e6,
            0xa4bbd149268eda35,
            0x648728f79da906ef,
            0xa59d876a69cd2af8,
            0x554c1f5e7873,
        ]));
        tmp.sub_assign(&Fq::new(BigInteger768::from(0)));
        assert_eq!(
            tmp,
            Fq::new(BigInteger768([
                0xb2b74f8c784b6de,
                0x1d5ae51c4d7a2afc,
                0xf424485fa9ab9789,
                0x39948fcd5f1ba445,
                0xa8635673ec8e6c81,
                0x19b75b83a122ee72,
                0x31e7f890d12e9443,
                0x3dbca4d50c6a02e6,
                0xa4bbd149268eda35,
                0x648728f79da906ef,
                0xa59d876a69cd2af8,
                0x554c1f5e7873,
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
        0xf90599b5974382d2,
        0xe7581c1924d6b303,
        0x55d7d7228dd30eb2,
        0xa47c9f5d998f2b51,
        0x411c00a50673af12,
        0x181d1518a9c7b25f,
        0x3c64f5fd46039bcf,
        0xf2e55ae09dbc1241,
        0xb7f12d1fb9d5c945,
        0x967c968236916e02,
        0x717e420231853795,
        0x636d6103f527,
    ]));
    tmp.mul_assign(&Fq::new(BigInteger768([
        0x43ea1602974be9df,
        0x9d89a1653778ac89,
        0x3d241871d23271cd,
        0x423a8c8dc1ec87fd,
        0xbd1b39df736ddb58,
        0x1d797c82d55bfa7e,
        0x37e34ae333d830a,
        0x86c6146b1b283b29,
        0x83834e34a16c2ac4,
        0x3ab52e777269366b,
        0xab8bf157064f27ed,
        0x1292685293a9a,
    ])));
    assert_eq!(
        tmp,
        Fq::new(BigInteger768([
            0x14ae353c059de348,
            0xeb6633d6a82837f5,
            0xf41cd69d4fbc892c,
            0xdb444089367e0b4a,
            0x961b7b8786357e80,
            0x1fbe8be6536a0371,
            0x465ac4f7242a2243,
            0xecd5870ff01c825,
            0x332f77d2b0c14161,
            0x6ad5a4ecbc496483,
            0x4cf4822137e9d7a2,
            0x14c07b7ffa456,
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
        0xd185e57f9afb8fb5,
        0x2b7fd95edee5f66f,
        0x9ea544d6b3e7c8e5,
        0x2301ffe498f19737,
        0xfc1034b27fe57524,
        0x9d51b22c45752994,
        0x6555c3b6097f83c7,
        0x9d504f54fb0a7b27,
        0x9b45f041912c19e9,
        0xc08807443d74051f,
        0xe49806cdd5773372,
        0x1c3d520d338be,
    ]));
    assert!(a.is_valid());
    a.square_in_place();
    assert_eq!(
        a,
        Fq::from_repr(BigInteger768([
            0x69661e7ce51d9de1,
            0x80e7a38a90970c77,
            0xd7b7136d8b7a3cb6,
            0x1f85d9a72700c5f1,
            0x23c598bde8f3bd79,
            0x72993e5df0b896c0,
            0x3d745a0701458c74,
            0x847527793ef4edcf,
            0x3e683b7b96c452f6,
            0xda49b7613adf8939,
            0xe426f68da1e89cc6,
            0x569dd56b8303,
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
        0xd9ddf9cba96cc287,
        0xd4a37d9a7f28d94c,
        0x2fbe515c4feaefd3,
        0x39cf73ff8cf508d7,
        0xbc4da230bc9f52c,
        0x3e1f0132e1f0851c,
        0x7e7c04b6a099574e,
        0xa7b18147273defee,
        0x4d41983dd4323832,
        0xf85193e73b78a121,
        0x9b111ff50e57db2d,
        0x186a7dede838e,
    ]));
    let a = a_squared.sqrt().unwrap();
    assert_eq!(
        a,
        Fq::new(BigInteger768([
            0x2246981b0859aa51,
            0x2c27b2a6d58c0be4,
            0xd12541e352f9bff1,
            0x70401d9ca2890cde,
            0xfe3a5678bfaeb0f7,
            0x7c1f5e5cfd935e01,
            0x4f7f6b949a430333,
            0x31a49135470aeee2,
            0xffff8d5b9eab0d02,
            0x989ffec98fb0ed77,
            0xccfebe585ad372c8,
            0x13dc68aa6edec,
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
fn test_fq_bytes() {
    let a = Fq::from_repr(BigInteger768([
        0xc2eb1a79dcc80fb,
        0x8b74b0048ccc0f85,
        0x1395ff13d91ce297,
        0x651b4e825addab48,
        0x8bb827faf35476a9,
        0x4139332c620ccbbc,
        0x81fac1a457bf1d18,
        0xc9e4a3f3eda6e267,
        0x7204a9a538c2a1e0,
        0xa08a8a70b10f136a,
        0x6fb145bda6920c42,
        0x2c646dfb3edc,
    ]));
    let a_b = to_bytes!(a).unwrap();
    let a_b_read = std::fs::read("src/fields/mnt4753/test_vec/mnt4753_tobyte").unwrap();
    assert_eq!(a_b, a_b_read);
    let a_read = Fq::read(a_b_read.as_slice()).unwrap();
    assert_eq!(a, a_read);
}

#[test]
fn test_convert_fq_fr() {
    use crate::fields::{
        convert,
        mnt4753::{Fr, FrParameters},
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
    let modulus_r = Fr::new(FrParameters::MODULUS);
    assert!(convert::<Fq>((modulus_r - &Fr::one()).write_bits()).is_err()); //Fr_Modulus - 1 is bigger than Fq modulus
}

#[test]
fn test_fq_root_of_unity() {
    assert_eq!(FqParameters::TWO_ADICITY, 15);
    assert_eq!(
        Fq::multiplicative_generator(),
        Fq::from_repr(BigInteger768::from(17))
    );
    assert_eq!(
        Fq::multiplicative_generator().pow([
            0x233EBD20C7BC48BD,
            0x4BE1C73AA8A459BA,
            0xA948C71020E33588,
            0xFC70D0B599D2ECE4,
            0x0B3B701E1B4B96A6,
            0xEF3B622FCEEDE430,
            0xDB1B33A249B342B5,
            0xB0E60FFB724BD141,
            0x5FDABD6FD1F2D92F,
            0x9B5B6FF32EA0B71F,
            0x882220452045DDC5,
            0x3898C5B25,
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
        Fq::from_repr(BigInteger768::from(13)).legendre()
    );
    assert_eq!(
        QuadraticResidue,
        Fq::from_repr(BigInteger768::from(169)).legendre()
    );

    let e = BigInteger768([
        0x489302efc996adf7,
        0x7b4bb81ad0f8d9ea,
        0x831b945e1cb94c65,
        0xde6cbbddcb71a21c,
        0xc4c288920781396c,
        0x1f510e8a5d0f9204,
        0x137c3afdd9394bc5,
        0x7a9b5336fea79b3b,
        0x5d045d7cf7e6e740,
        0x78ce9be361f75af2,
        0x72442a1e6ff0a47f,
        0xa813136f81ec,
    ]);
    assert_eq!(QuadraticNonResidue, Fq::from_repr(e).legendre());
    let e = BigInteger768([
        0xc467a286665a3a01,
        0x746bb22010770da0,
        0x199fd7b97ae3dde8,
        0x5f5803abb402a9a1,
        0xea7b59755662a360,
        0xb0389a63076a2e8d,
        0x20e406e2cbb7362f,
        0x50c0fbbcf08074db,
        0x66d856be449cdfbb,
        0x567eadc74aa00b15,
        0x4412e4c5b9ce9aae,
        0x17580790e3633,
    ]);
    assert_eq!(QuadraticResidue, Fq::from_repr(e).legendre());
}

#[test]
fn test_fq2_ordering() {
    let mut a = Fq2::new(Fq::zero(), Fq::zero());

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
fn test_fq2_basics() {
    assert_eq!(Fq2::new(Fq::zero(), Fq::zero(),), Fq2::zero());
    assert_eq!(Fq2::new(Fq::one(), Fq::zero(),), Fq2::one());
    assert!(Fq2::zero().is_zero());
    assert!(!Fq2::one().is_zero());
    assert!(!Fq2::new(Fq::zero(), Fq::one(),).is_zero());
}

#[test]
fn test_fq2_squaring() {
    // i = sqrt(13) in mnt4_753 fq2

    //(8+i)^2 = 77 + 16i
    let a = Fq2::new(Fq::from_repr(BigInteger768::from(8)), Fq::one()).square();
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768::from(77)),
            Fq::from_repr(BigInteger768::from(16))
        )
    );

    //i^2 = 13
    let a = Fq2::new(Fq::zero(), Fq::one()).square();
    assert_eq!(
        a,
        Fq2::new(Fq::from_repr(BigInteger768::from(13)), Fq::zero())
    );

    let mut a = Fq2::new(
        Fq::from_repr(BigInteger768([
            0xff128f8b944c48c5,
            0x73351fc9610b2fc7,
            0x2a7ec9853b6149c2,
            0x829544b8e70c8324,
            0x90ca8df680dbb3cb,
            0x97b890988d408de,
            0xad34ee124fc9e1e5,
            0x28112f0b0c052e9,
            0x7efb3d6e48b56c29,
            0x42ca7e6b53f59fe4,
            0x59902c4c0c7c2794,
            0x3a53cfef3e95,
        ])),
        Fq::from_repr(BigInteger768([
            0x97e9cfc0c257c70b,
            0x7df17f22f3a7cd5,
            0xa17378ac19414061,
            0x19b57a91e3078396,
            0xa89301cd176713b0,
            0x779fd05b05cc6f7e,
            0xa78d0b554b342fe7,
            0xb6c8fdc260726a59,
            0xd46bb7e9849a0674,
            0x32135600cc7c33ef,
            0x9651b872b7c8c88d,
            0x11cc6249d7977,
        ])),
    );
    a.square_in_place();
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0x7e776a6249e3692f,
                0xe94d32ccee6a0fa0,
                0x6f7be94f4856f904,
                0xa32707f3e6c8681e,
                0xb703912c9df1b826,
                0x8695d8a2b3f9a2d7,
                0x4fc12ff6b2adc9e3,
                0x60c24677cc0ff766,
                0xbfc9fe3df21beb2f,
                0x51fb8916485165d8,
                0xb5dd8fc1fc3b6b62,
                0x104064032c6cb,
            ])),
            Fq::from_repr(BigInteger768([
                0x55da822db96624ac,
                0x55b289dc1366a1de,
                0xd668c50654d6f919,
                0xbf3c0956bb907919,
                0xf1830ee028bfbfc0,
                0xa86fbe0e8aa0a18c,
                0xaefb71371ed4a03e,
                0x3d39ed069544eefb,
                0xfaf8197459cb76dc,
                0xf3e574a54c55b2d,
                0x67da7270a707f7f5,
                0xbd0be9a5a08a,
            ])),
        )
    );
}

#[test]
fn test_fq2_mul() {
    let mut a = Fq2::new(
        Fq::from_repr(BigInteger768([
            0x2afdb2e15999761d,
            0xbe74b32b3fc6c1a5,
            0x98b5af1a5ae56540,
            0xbba4a5a2d14a9496,
            0x399e0ec60182c40,
            0x1a9fdc1d6c67c284,
            0x65b4509da003e4cd,
            0x7f1b656f0f977c2f,
            0x3844e026770c0ca9,
            0xd7a204e04badafd0,
            0x19cf945757c30fb0,
            0x18ac8efa22b1d,
        ])),
        Fq::from_repr(BigInteger768([
            0x859a004149732671,
            0x375fb8c3d741b4ec,
            0x8abe05fa793ec11a,
            0xa2d18fddf5c30ba0,
            0x5376af17184bb5cc,
            0x98559a5bc9071c94,
            0x91c9e09032a76bb,
            0x12d08be52bab51be,
            0x4ca7b42ec1c58def,
            0x283e9cecfeaae3e,
            0x82b3528d7699d664,
            0x15f0238583afa,
        ])),
    );
    a.mul_assign(&Fq2::new(
        Fq::from_repr(BigInteger768([
            0x182ee151232e14cf,
            0x5dd5c035904f4c4b,
            0x72655ce8807a17e1,
            0x9fd5dd0dd3e1f830,
            0xf28ce0175b26c7d4,
            0x123ba7eeeecc6ef9,
            0xa4e8cc882555844a,
            0xb734c8a6ee7f5732,
            0x2589347b7c6b2a4a,
            0x28866daaeef933a3,
            0xa5ec12d10c6c461,
            0x16772bdde2b5f,
        ])),
        Fq::from_repr(BigInteger768([
            0x90e899987662bb27,
            0xf05be6e5e0216763,
            0xc814bc0792075341,
            0x112e6ef80d2378f9,
            0x6d14622921e4e62,
            0xfaeea70c1909193e,
            0x1d48b9e0ed0434d0,
            0x892966119014a798,
            0xc7dc187e4e1d5461,
            0x5eb7beebd479c963,
            0x332549b47f987f3a,
            0xe34afdcf1d1d,
        ])),
    ));
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0xecb889331ad1571e,
                0x82678d008679c1cb,
                0x11741dc59f073796,
                0x747bb020f4f07faa,
                0xa5fe03aa50ea7ac0,
                0x43f301fb870ec899,
                0xbe67409dcc2a5cc4,
                0xb28630629c388aec,
                0xb0ac2892f1d8533d,
                0x545103503e23eb1,
                0xa66dfe734b1c7559,
                0x18aa7cd804a0,
            ])),
            Fq::from_repr(BigInteger768([
                0xaad2909b79f04c49,
                0x1a0985bb1bb7df0b,
                0xd308b5e6a70a91,
                0x7c50582215541503,
                0x4ebd6cbcfff753d8,
                0x8b3951d333a1e62b,
                0xe36bd3445f646239,
                0xc38e569980cd3b0a,
                0x20971efaf96a7221,
                0xd9732b87dd87e2f2,
                0x89fd7bc4b1e0d206,
                0x8ad6da1cadf6,
            ])),
        )
    );
}

#[test]
fn test_fq2_inverse() {
    assert!(Fq2::zero().inverse().is_none());

    let a = Fq2::new(
        Fq::from_repr(BigInteger768([
            0x964ee7cea60eb57a,
            0x6a20aa04ade93fef,
            0xacad401a731d63b8,
            0x9b17c20453c8ae50,
            0x90a161e54914c977,
            0xb47208162dbabb01,
            0xc9fb2d426c6efeb0,
            0x6d8dd6c568d71255,
            0xa1d593606b258533,
            0x2ab40dc82fb3e5b2,
            0xb3fd8dc7652f42fd,
            0xdd826c52197d,
        ])),
        Fq::from_repr(BigInteger768([
            0x3f8dc7345415c891,
            0x6e349462399200a9,
            0x5c6297acf72ec991,
            0xbb86c426e88a4703,
            0xf79cb1fdb656a004,
            0x3add3e19fdef1ca1,
            0x93d52b31be1bb9ab,
            0x51b69e5a6351de11,
            0x2d4c21d857c7d54,
            0x6133841cfe93454a,
            0x25e4324f4bed09d,
            0xd4a9ceb9414e,
        ])),
    );
    let a = a.inverse().unwrap();
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0xa9170fd89427009f,
                0x3007b9970b056fe0,
                0xb6fdba1a0e158619,
                0x8a2b68fb7c527ea3,
                0x5e27ea1d2fa42e38,
                0x273c25a3af332812,
                0xfb4c5e4663edf31f,
                0xc622624050d88ce0,
                0xc4bc1e449e9330f0,
                0x299d388979feb63a,
                0x9e33166d1029311d,
                0x15c016a382d93,
            ])),
            Fq::from_repr(BigInteger768([
                0x54cb5db3c0e0dfd,
                0xcc69f7964588fe33,
                0x8d1746178ee824fe,
                0x2ad9cfb2c866b0f2,
                0x7ce0518badf35a70,
                0x4908b6bff1aa2b0f,
                0x7700ee3e86ebe18f,
                0xb824ee7d86a120e3,
                0xeb69a93c7c547c8,
                0x658da3ee4426f138,
                0x2ff3cf1bc2054fbb,
                0x10913f2d35c0b,
            ])),
        )
    );
}

#[test]
fn test_fq2_addition() {
    let mut a = Fq2::new(
        Fq::from_repr(BigInteger768([
            0x9c6a036d6c0e6b2,
            0xe9622ccb4be29f0,
            0xd1090e2133506233,
            0xe2be2eca38e9223b,
            0xb215917925a0b6bb,
            0x9b4543407e543059,
            0xf1ac1204181afa12,
            0x14d664ce1b1ec9fa,
            0xd3965b381be984e9,
            0xe858e69aca0e179,
            0x3d0573773a04373c,
            0xad05b8045478,
        ])),
        Fq::from_repr(BigInteger768([
            0x6e61832ccce96f9,
            0xe7e230cea04d41e9,
            0x36feb392cebdfe0,
            0xb80d4df2196fa340,
            0x5849ab55060d640e,
            0x92332e56d682f3b,
            0x73678f1c037bfc3,
            0xf90690d8831b1b1d,
            0xfc9aa8b7c28a18cb,
            0x89725c4fcb551e4f,
            0xb9a41c05d4daae61,
            0xa6cc7717bc8d,
        ])),
    );
    a.add_assign(&Fq2::new(
        Fq::from_repr(BigInteger768([
            0x989b61d995737933,
            0xd749df4a5ff4f0b,
            0x48fbe68eb6a94b7a,
            0xec7ac71f4db92487,
            0xa785e0a0ad64a59a,
            0x5c954cb8bfa22fe6,
            0xe3a2c18a1e689e2c,
            0x798e452560b52cbb,
            0xe829a8718cc5fd11,
            0x197b35482311da0e,
            0x9602e4eb487dc7d1,
            0x1703bbbf3925b,
        ])),
        Fq::from_repr(BigInteger768([
            0x33eb5194f9e480e8,
            0x3802e435600e705f,
            0xce1b3a556e69d362,
            0x17b3496b34b071b3,
            0x961047d03fabb74,
            0xe17e03b8314d400d,
            0x836ac49d3810035d,
            0x27069fab144e67e,
            0x170dddb329465830,
            0xaee5d09dcf6bc44,
            0xa6f2b42ef5824b55,
            0x8f00b6755561,
        ])),
    ));
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0x43d19e3247d5dfe4,
                0x386d6c6f2de0675c,
                0xb67ce43e4f3587bc,
                0x66de2900102ff21e,
                0xa18c647407b1de1e,
                0x46c2a8824bde5aa2,
                0x3b7daeb49528a0a1,
                0x8666f0cd93330929,
                0x5d081ab03c17a987,
                0x70072c6174230b9b,
                0xc2e5c83f939f315f,
                0x587b466522c2,
            ])),
            Fq::from_repr(BigInteger768([
                0x3ad169c7c6b317e1,
                0x1fe51504005bb248,
                0xd18b258e9b55b343,
                0xcfc0975d4e2014f3,
                0x61aaafd20a081f82,
                0xeaa1369d9eb56f48,
                0x8aa13d8ef847c320,
                0xfb76fad33460019b,
                0x13a8866aebd070fb,
                0x9460b959a84bda94,
                0x6096d034ca5cf9b6,
                0x135cd2d8d11ef,
            ])),
        )
    );
}

#[test]
fn test_fq2_subtraction() {
    let mut a = Fq2::new(
        Fq::from_repr(BigInteger768([
            0xd223db5aa01be301,
            0xb43df8745c408868,
            0xda58342db6a6e003,
            0x3551f61a382a5c2e,
            0xa2eb7284b58e8994,
            0xb26cb3d1474a47aa,
            0xc2fd220690df29bc,
            0x65ce308f11635aab,
            0xb996006a2d77ee54,
            0x6f6bb892d2ae31f6,
            0xc8af668ef80502b6,
            0xf3665d7d89e5,
        ])),
        Fq::from_repr(BigInteger768([
            0x95cb398cb87c0a77,
            0x8f3c32f14300d583,
            0xb6a93e2d244afc4c,
            0xc9180ec333d1582b,
            0x81858a29598551f2,
            0xcece72a28273e79e,
            0x24fd7d0687680c73,
            0xc8e2dab05c00b126,
            0x1b3626091256f084,
            0x31356f806cce35ef,
            0x62b546ca1c5c9625,
            0x11e6a92ae3475,
        ])),
    );
    a.sub_assign(&Fq2::new(
        Fq::from_repr(BigInteger768([
            0x58e8fabc65a52629,
            0x813ab8eb0453c060,
            0x4ed6e1b7df50afd8,
            0xc924a861da12942f,
            0x9ec615461d151df7,
            0x9aa300afbf8f0c9f,
            0x5aca78259b54103e,
            0xb5a526fa22b2ba5b,
            0x9d493bf936d88b81,
            0xf2b01cb9c763b92,
            0x99144cf7e765cd17,
            0x651ce9fa1806,
        ])),
        Fq::from_repr(BigInteger768([
            0x3583e181dd94c2b1,
            0xd60e335793fda76f,
            0xd640106688f4c053,
            0xe4dde19abd9bf22c,
            0x1b50bcd3b8e65405,
            0x483e6b331c38ddd4,
            0x2b5a5514814add9d,
            0x1460bb60835d99b6,
            0x3d7edd36ad73cba2,
            0xd2ba566386ceff49,
            0x2510954f48b81a0,
            0x1163a5829bd3a,
        ])),
    ));
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0x793ae09e3a76bcd8,
                0x33033f8957ecc808,
                0x8b815275d756302b,
                0x6c2d4db85e17c7ff,
                0x4255d3e98796b9c,
                0x17c9b32187bb3b0b,
                0x6832a9e0f58b197e,
                0xb0290994eeb0a050,
                0x1c4cc470f69f62d2,
                0x6040b6c73637f664,
                0x2f9b1997109f359f,
                0x8e49738371df,
            ])),
            Fq::from_repr(BigInteger768([
                0x6047580adae747c6,
                0xb92dff99af032e14,
                0xe0692dc69b563bf8,
                0xe43a2d28763565fe,
                0x6634cd55a09efdec,
                0x8690076f663b09ca,
                0xf9a327f2061d2ed6,
                0xb4821f4fd8a3176f,
                0xddb748d264e324e2,
                0x5e7b191ce5ff36a5,
                0x60643d7527d11484,
                0x8303a84773b,
            ])),
        )
    );
}

#[test]
fn test_fq2_negation() {
    let mut a = Fq2::new(
        Fq::from_repr(BigInteger768([
            0xacfb56ae6b5e56c5,
            0x25b0948724d89058,
            0xfb36d5a1676d8cf0,
            0xb323d8e527ee2e3a,
            0x381dcb3112b73661,
            0xe223c17603e0f4cd,
            0xd1892e0ca7b20e00,
            0x1fec661b9b11f52e,
            0x80eab1b3db98720d,
            0x35d5144f89e05606,
            0xee744410b3f8a3f3,
            0xf8ff68170059,
        ])),
        Fq::from_repr(BigInteger768([
            0x26c4d98b160b718d,
            0x8d394f31abe4784,
            0xfb98c2b6c2052e84,
            0x4382195183aa598,
            0x3186453d2058019,
            0x9889296f1d8715f0,
            0x86437bcffbe3c11c,
            0x9bc7d14c4a915d69,
            0xdc71f4fdc3a37922,
            0x8df9c47d66790331,
            0xbe517812bd9a9752,
            0x1b7cb0803f1c0,
        ])),
    );
    a = -a;
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0xb1950d2fb900293c,
                0xbdecbfcb08048146,
                0x68513ad033569900,
                0xb536f4044e842669,
                0x7ff14274b89c47d6,
                0xcef42600ee3710d0,
                0xc847f6ccf9a8e99c,
                0xe811530a4d8ef85e,
                0xddcd374590ff6665,
                0x82248300d1af59e6,
                0x21ae4c123aea29ba,
                0xcbc6c57bc3b7,
            ])),
            Fq::from_repr(BigInteger768([
                0x37cb8a530e530e74,
                0xdac9bf5f121eca1b,
                0x67ef4dbad8bef76c,
                0x6422ab545e37af0b,
                0xb4f6a951f94dfe1f,
                0x188ebe07d490efad,
                0x138da909a5773681,
                0x6c35e7d99e0f9024,
                0x8245f3fba8f45f50,
                0x29ffd2d2f516acbb,
                0x51d118103148365b,
                0xcfb258ed250,
            ])),
        )
    );
}

#[test]
fn test_fq2_doubling() {
    let mut a = Fq2::new(
        Fq::from_repr(BigInteger768([
            0x5de6a90e1dabcd6f,
            0x91232ff1f83871e1,
            0xa3b7cf9905943cda,
            0x2ae2448c879b402b,
            0xae8dbfbb420e533e,
            0x37e2505227b2488a,
            0x1b22f1a4d29c13e4,
            0xa6ffb961595f2f41,
            0xc0e32085f70a1843,
            0xab5734b40aeafeaa,
            0x4e35a1ed7b854060,
            0x1b9215305af80,
        ])),
        Fq::from_repr(BigInteger768([
            0x496546eafb30bf92,
            0x9c79448c3429e603,
            0x9d987088402f672f,
            0x21271066822a7fe3,
            0xed1ba13befaeb1e6,
            0x329f98ccf9f342c2,
            0xacf3614d1bb82214,
            0x5934a9705c3c91cc,
            0xfdfb138fb457b4fd,
            0x6eac21a419a97545,
            0xdf52d3010537ce54,
            0x10b44e2b445c4,
        ])),
    );
    a.double_in_place();
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0x5d3cee3e16f91add,
                0x3ea90b91c393d223,
                0xe3e78ec0706453c4,
                0xed69bc2f98c42bb2,
                0xa50c71d0b8c92843,
                0xbeacb92d5d4c8b77,
                0x9c74be7003dd302a,
                0x4601b99cca1d70f4,
                0x230e5812817c5814,
                0x9eb4d217ba464d68,
                0x8c48b3b80827b313,
                0x1ad7c78789aef,
            ])),
            Fq::from_repr(BigInteger768([
                0x343a29f7d202ff23,
                0x555534c63b76ba67,
                0xd7a8d09ee59aa86e,
                0xd9f353e38de2ab22,
                0x222834d21409e593,
                0xb4274a2301ce7fe8,
                0xc0159dc096154c8a,
                0xaa6b99bacfd8360b,
                0x9d3e3e25fc179187,
                0x255eabf7d7c33a9e,
                0xae8315df1b8ccefb,
                0x51c397d5c778,
            ])),
        )
    );
}

#[test]
fn test_fq2_frobenius_map() {
    let mut a = Fq2::new(
        Fq::from_repr(BigInteger768([
            0xe4cc1c461ee5b8af,
            0xbfdde425f805dba7,
            0x6583adec0dff5830,
            0x595619386f074c4c,
            0x70bf2940f4b13a47,
            0x9da1accdeaa0ffed,
            0xb4463c694091e28d,
            0x8099a1ea8431ee7a,
            0x68c85b07984826ba,
            0x29318418cb52a2fd,
            0xd79b651f5823264d,
            0x12cd3bc49a35b,
        ])),
        Fq::from_repr(BigInteger768([
            0x3fbc766f5be800c,
            0x9f560a6739cda688,
            0xec0e383d4c2b9af7,
            0x7fed2657aa03f5f1,
            0x3d4f91a837a35136,
            0xcfd09e93ffbee91a,
            0x7636125662e8020f,
            0xb4303ceb7a4cfb00,
            0x4204c128b8670d5b,
            0x63f63e46eb59aaa3,
            0x987831d5a07208cc,
            0x180bca61c183b,
        ])),
    );
    a.frobenius_map(0);
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0xe4cc1c461ee5b8af,
                0xbfdde425f805dba7,
                0x6583adec0dff5830,
                0x595619386f074c4c,
                0x70bf2940f4b13a47,
                0x9da1accdeaa0ffed,
                0xb4463c694091e28d,
                0x8099a1ea8431ee7a,
                0x68c85b07984826ba,
                0x29318418cb52a2fd,
                0xd79b651f5823264d,
                0x12cd3bc49a35b,
            ])),
            Fq::from_repr(BigInteger768([
                0x3fbc766f5be800c,
                0x9f560a6739cda688,
                0xec0e383d4c2b9af7,
                0x7fed2657aa03f5f1,
                0x3d4f91a837a35136,
                0xcfd09e93ffbee91a,
                0x7636125662e8020f,
                0xb4303ceb7a4cfb00,
                0x4204c128b8670d5b,
                0x63f63e46eb59aaa3,
                0x987831d5a07208cc,
                0x180bca61c183b,
            ])),
        )
    );
    a.frobenius_map(1);
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0xe4cc1c461ee5b8af,
                0xbfdde425f805dba7,
                0x6583adec0dff5830,
                0x595619386f074c4c,
                0x70bf2940f4b13a47,
                0x9da1accdeaa0ffed,
                0xb4463c694091e28d,
                0x8099a1ea8431ee7a,
                0x68c85b07984826ba,
                0x29318418cb52a2fd,
                0xd79b651f5823264d,
                0x12cd3bc49a35b,
            ])),
            Fq::from_repr(BigInteger768([
                0x5a949c772e9ffff5,
                0x444749eaf30f6b17,
                0x7779d8344e988af9,
                0xe86da691cc6e5eb2,
                0x7abf7bfd93b02d01,
                0xe14748e2f2591c83,
                0x239b12833e72f58d,
                0x53cd7c3a6e53f28d,
                0x1cb327d0b430cb17,
                0x540359097036054a,
                0x77aa5e4d4e70c4e1,
                0x44098776abd5,
            ])),
        )
    );
    a.frobenius_map(1);
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0xe4cc1c461ee5b8af,
                0xbfdde425f805dba7,
                0x6583adec0dff5830,
                0x595619386f074c4c,
                0x70bf2940f4b13a47,
                0x9da1accdeaa0ffed,
                0xb4463c694091e28d,
                0x8099a1ea8431ee7a,
                0x68c85b07984826ba,
                0x29318418cb52a2fd,
                0xd79b651f5823264d,
                0x12cd3bc49a35b,
            ])),
            Fq::from_repr(BigInteger768([
                0x3fbc766f5be800c,
                0x9f560a6739cda688,
                0xec0e383d4c2b9af7,
                0x7fed2657aa03f5f1,
                0x3d4f91a837a35136,
                0xcfd09e93ffbee91a,
                0x7636125662e8020f,
                0xb4303ceb7a4cfb00,
                0x4204c128b8670d5b,
                0x63f63e46eb59aaa3,
                0x987831d5a07208cc,
                0x180bca61c183b,
            ])),
        )
    );
    a.frobenius_map(2);
    assert_eq!(
        a,
        Fq2::new(
            Fq::from_repr(BigInteger768([
                0xe4cc1c461ee5b8af,
                0xbfdde425f805dba7,
                0x6583adec0dff5830,
                0x595619386f074c4c,
                0x70bf2940f4b13a47,
                0x9da1accdeaa0ffed,
                0xb4463c694091e28d,
                0x8099a1ea8431ee7a,
                0x68c85b07984826ba,
                0x29318418cb52a2fd,
                0xd79b651f5823264d,
                0x12cd3bc49a35b,
            ])),
            Fq::from_repr(BigInteger768([
                0x3fbc766f5be800c,
                0x9f560a6739cda688,
                0xec0e383d4c2b9af7,
                0x7fed2657aa03f5f1,
                0x3d4f91a837a35136,
                0xcfd09e93ffbee91a,
                0x7636125662e8020f,
                0xb4303ceb7a4cfb00,
                0x4204c128b8670d5b,
                0x63f63e46eb59aaa3,
                0x987831d5a07208cc,
                0x180bca61c183b,
            ])),
        )
    );
}

#[test]
fn test_fq2_legendre() {
    use crate::fields::LegendreSymbol::*;

    assert_eq!(Zero, Fq2::zero().legendre());
    assert_eq!(QuadraticNonResidue, Fq2Parameters::NONRESIDUE.legendre());
    assert_eq!(QuadraticNonResidue, Fq4Parameters::NONRESIDUE.legendre());

    // i^2 = -1
    let mut m1 = -Fq2::one();
    assert_eq!(QuadraticResidue, m1.legendre());
    m1 = Fq4Parameters::mul_fp2_by_nonresidue(&m1);
    assert_eq!(QuadraticNonResidue, m1.legendre());
}

#[test]
fn test_fq2_mul_nonresidue() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let nqr = Fq2::new(Fq::zero(), Fq::one());

    for _ in 0..1000 {
        let mut a = Fq2::rand(&mut rng);
        let mut b = a;
        a = Fq4Parameters::mul_fp2_by_nonresidue(&a);
        b.mul_assign(&nqr);

        assert_eq!(a, b);
    }
}

#[test]
fn test_fq4_mul_by_023() {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    for _ in 0..1000 {
        let c0 = Fq::rand(&mut rng);
        let c2 = Fq::rand(&mut rng);
        let c3 = Fq::rand(&mut rng);
        let to_mul = Fq4::new(Fq2::new(c0, Fq::zero()), Fq2::new(c2, c3));
        let a = Fq4::rand(&mut rng);
        let mut b = a;

        b.mul_assign(&to_mul);

        assert_eq!(a.mul_by_023(&to_mul), b);
    }
}
