use crate::{
    field_new,
    biginteger::BigInteger768 as BigInteger,
    fields::{
        fp6_2over3::{Fp6, Fp6Parameters},
        mnt6753::{
            fq753b::Fq,
            fq3::{Fq3, Fq3b7Parameters},
        },
    },
};

pub type Fq6 = Fp6<Fq6b7Parameters>;

pub struct Fq6b7Parameters;

impl Fp6Parameters for Fq6b7Parameters {
    type Fp3Params = Fq3b7Parameters;

// 15664CB51F657D6070F328598582FA2EC8A689E11FE486C43D93CD263C4661714868F67A7CF4BDB24398ACBB12E86EB5409814CC7808279A5240144C32F00B32E3562232D379FF502478C2D062EDF76AA224C87975C57CB1308DB6D9215A2
    const NONRESIDUE: Fq3 = field_new!(Fq3, 
        field_new!(Fq, BigInteger([
            0xCB1308DB6D9215A2,
            0x76AA224C87975C57,
            0xF502478C2D062EDF,
            0xB32E3562232D379F,
            0x79A5240144C32F00,
            0xEB5409814CC78082,
            0xDB24398ACBB12E86,
            0x1714868F67A7CF4B,
            0x6C43D93CD263C466,
            0xA2EC8A689E11FE48,
            0xD6070F328598582F,
            0x15664CB51F657,
        ])),
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    );

    const FROBENIUS_COEFF_FP6_C1: [Fq; 6] = [
        field_new!(Fq, BigInteger([
            0xb99680147fff6f42,
            0x4eb16817b589cea8,
            0xa1ebd2d90c79e179,
            0xf725caec549c0da,
            0xab0c4ee6d3e6dad4,
            0x9fbca908de0ccb62,
            0x320c3bb713338498,
            0x598b4302d2f00a62,
            0x4074c9cbfd8ca621,
            0xfa47edb3865e88c,
            0x95455fb31ff9a195,
            0x7b479ec8e242,
        ])),
// Fq6 Frobenius_coeffs_c1[1]  19897d1eb1ca04df5e7055d775a97c45deb6fc774821ee7cb22af3bb2a0123f3ec8159c1da9be186cec139ca9d216943d1b65fbd3ae53fd566bac405fc44b4ed706302c0e2b723da2c7145efc7815622b1a03ad8db3c724fd77ccbf95f5c2
        field_new!(Fq, BigInteger([
            0x24fd77ccbf95f5c2,
            0x622b1a03ad8db3c7,
            0x3da2c7145efc7815,
            0x4ed706302c0e2b72,
            0xfd566bac405fc44b,
            0x943d1b65fbd3ae53,
            0x186cec139ca9d216,
            0x3f3ec8159c1da9be,
            0xe7cb22af3bb2a012,
            0xc45deb6fc774821e,
            0x4df5e7055d775a97,
            0x19897d1eb1ca0,
        ])),
// Fq6 Frobenius_coeffs_c1[2]  11d5033223a5db8b087523d7db902b4b96c948f0e9992a75658e33e25f9f0e5b38512c92d9f5be660b05c89764d7df480725d1dc6e2f1524a1cc56c78e9773f64a98166c46a979bb6f43b5282969c1379b1ebf803e51e6b66f7b83f968680
        field_new!(Fq, BigInteger([
            0x6b66f7b83f968680,
            0x1379b1ebf803e51e,
            0x9bb6f43b5282969c,
            0x3f64a98166c46a97,
            0x524a1cc56c78e977,
            0xf480725d1dc6e2f1,
            0xe660b05c89764d7d,
            0xe5b38512c92d9f5b,
            0xa75658e33e25f9f0,
            0xb4b96c948f0e9992,
            0xb8b087523d7db902,
            0x11d5033223a5d,
        ])),
// Fq6 Frobenius_coeffs_c1[3]  1497e8ec9e1ce7add306fcee92c18a85518752329c7611e431f2d6f0b3251ae72762315b0e32b67c4e9228e27730512afb31fea4cde3893df7bad553a4b62aa6d9cc76f4f79ca34d7aee33286761dffef30ff5a176ba71f70f6cdc00090bf
        field_new!(Fq, BigInteger([
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
        ])),
// Fq6 Frobenius_coeffs_c1[4]  2c2e5ba7a770c22ca91d916b7315f39babe0941b2dce76ecc64a30e53860c8bef1104c8343cf816438c604b125871e2f40c2cc85fb4741955ee7e8c161eb6b08f346088b0f329920baa7e003df81ec757f1362138688b409ff15806a0a3f
        field_new!(Fq, BigInteger([
            0xb409ff15806a0a3f,
            0xec757f1362138688,
            0x9920baa7e003df81,
            0x6b08f346088b0f32,
            0x41955ee7e8c161eb,
            0x1e2f40c2cc85fb47,
            0x816438c604b12587,
            0xc8bef1104c8343cf,
            0x76ecc64a30e53860,
            0xf39babe0941b2dce,
            0xc22ca91d916b7315,
            0x2c2e5ba7a770,
        ])),
// Fq6 Frobenius_coeffs_c1[5]  a775fa7089b3577208d0b16514ab03402abbcc81165ab76190162e71de82224a34131f734e31b370747d17e4aa1fbdebe9cbaa92c6a9eca1adcebca83cbf7a7b4ff4cdd4d00d3b0c8d80ec7dc0fb3b26e72b179d55316da07f2a00697981
        field_new!(Fq, BigInteger([
            0x6da07f2a00697981,
            0x3b26e72b179d5531,
            0x3b0c8d80ec7dc0fb,
            0x7a7b4ff4cdd4d00d,
            0xeca1adcebca83cbf,
            0xbdebe9cbaa92c6a9,
            0xb370747d17e4aa1f,
            0x224a34131f734e31,
            0xb76190162e71de82,
            0x3402abbcc81165a,
            0x577208d0b16514ab,
            0xa775fa7089b3,
        ])),
    ];
}
