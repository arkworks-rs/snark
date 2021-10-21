use crate::{
    biginteger::BigInteger768 as BigInteger,
    field_new,
    fields::{
        fp6_2over3::{Fp6, Fp6Parameters},
        mnt6753::{
            fq::{Fq, FQ_ONE, FQ_ZERO},
            fq3::{Fq3, Fq3Parameters},
        },
    },
};

pub type Fq6 = Fp6<Fq6Parameters>;

pub struct Fq6Parameters;

impl Fp6Parameters for Fq6Parameters {
    type Fp3Params = Fq3Parameters;

    const NONRESIDUE: Fq3 = field_new!(Fq3, FQ_ZERO, FQ_ONE, FQ_ZERO);

    const FROBENIUS_COEFF_FP6_C1: &'static [Fq] = &[
        //alpha^((q^0 - 1)/ 6) = 1
        field_new!(
            Fq,
            BigInteger([
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
            ])
        ),
        //alpha^((q^1 - 1)/ 6)
        field_new!(
            Fq,
            BigInteger([
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
            ])
        ),
        //alpha^((q^2 - 1)/ 6)
        field_new!(
            Fq,
            BigInteger([
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
            ])
        ),
        //alpha^((q^3 - 1)/ 6)
        field_new!(
            Fq,
            BigInteger([
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
            ])
        ),
        //alpha^((q^4 - 1)/ 6)
        field_new!(
            Fq,
            BigInteger([
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
            ])
        ),
        //alpha^((q^5 - 1)/ 6)
        field_new!(
            Fq,
            BigInteger([
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
            ])
        ),
    ];
}
