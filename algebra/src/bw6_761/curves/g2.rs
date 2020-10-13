use crate::{
    biginteger::{BigInteger384, BigInteger768},
    bw6_761::{Fq, Fr},
    curves::{
        models::{ModelParameters, SWModelParameters},
        short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        GLVParameters,
    },
    field_new,
    fields::PrimeField,
};

pub type G2Affine = GroupAffine<Parameters>;
pub type G2Projective = GroupProjective<Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl ModelParameters for Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl GLVParameters for Parameters {
    type WideBigInt = BigInteger768;

    /// phi((x, y)) = (\omega x, y)
    /// \omega = 0x531dc16c6ecd27aa846c61024e4cca6c1f31e53bd9603c2d17be416c5e44
    /// 26ee4a737f73b6f952ab5e57926fa701848e0a235a0a398300c65759fc4518315
    /// 1f2f082d4dcb5e37cb6290012d96f8819c547ba8a4000002f962140000000002a
    const OMEGA: Fq = field_new!(
        Fq,
        BigInteger768([
            9193734820520314185,
            15390913228415833887,
            5309822015742495676,
            5431732283202763350,
            17252325881282386417,
            298854800984767943,
            15252629665615712253,
            11476276919959978448,
            6617989123466214626,
            293279592164056124,
            3271178847573361778,
            76563709148138387
        ])
    );

    /// lambda in Z s.t. phi(P) = lambda*P for all P
    /// \lambda = 0x9b3af05dd14f6ec619aaf7d34594aabc5ed1347970dec00452217cc900000008508c00000000001
    const LAMBDA: Self::ScalarField = field_new!(
        Fr,
        (BigInteger384([
            15766275933608376691,
            15635974902606112666,
            1934946774703877852,
            18129354943882397960,
            15437979634065614942,
            101285514078273488
        ]))
    );
    /// |round(B1 * R / n)|
    const Q2: <Self::ScalarField as PrimeField>::BigInt = BigInteger384([
        14430678704534329733,
        14479735877321354361,
        6958676793196883088,
        21,
        0,
        0,
    ]);
    const B1: <Self::ScalarField as PrimeField>::BigInt = BigInteger384([
        9586122913090633729,
        9963140610363752448,
        2588746559005780992,
        0,
        0,
        0,
    ]);
    const B1_IS_NEG: bool = false;
    /// |round(B2 * R / n)|
    const Q1: <Self::ScalarField as PrimeField>::BigInt = BigInteger384([
        11941976086484053770,
        4826578625773784813,
        2319558931065627696,
        7,
        0,
        0,
    ]);
    const B2: <Self::ScalarField as PrimeField>::BigInt = BigInteger384([
        6390748608727089153,
        3321046870121250816,
        862915519668593664,
        0,
        0,
        0,
    ]);
    const R_BITS: u32 = 384;
}

impl SWModelParameters for Parameters {
    /// COEFF_A = 0
    #[rustfmt::skip]

    const COEFF_A: Fq = field_new!(Fq, BigInteger768([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));

    /// COEFF_B = 4
    #[rustfmt::skip]
    const COEFF_B: Fq = field_new!(Fq, BigInteger768([
        0x136efffffffe16c9,
        0x82cf5a6dcffe3319,
        0x6458c05f1f0e0741,
        0xd10ae605e52a4eda,
        0x41ca591c0266e100,
        0x7d0fd59c3626929f,
        0x9967dc004d00c112,
        0x1ccff9c033379af5,
        0x9ad6ec10a23f63af,
        0x5cec11251a72c235,
        0x8d18b1ae789ba83e,
        10403402007434220,
    ]));

    /// COFACTOR =
    /// 26642435879335816683987677701488073867751118270052650655942102502312977592501693353047140953112195348280268661194869
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        0x3de5800000000075,
        0x832ba4061000003b,
        0xc61c554757551c0c,
        0xc856a0853c9db94c,
        0x2c77d5ac34cb12ef,
        0xad1972339049ce76,
    ];

    /// COFACTOR^(-1) mod r =
    /// 214911522365886453591244899095480747723790054550866810551297776298664428889000553861210287833206024638187939842124
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger384([
        14378295991815829998,
        14586153992421458638,
        9788477762582722914,
        12654821707953664524,
        15185631607604703397,
        26723985783783076,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G2_GENERATOR_X, G2_GENERATOR_Y);
    #[inline(always)]
    fn mul_by_a(_elem: &Self::BaseField) -> Self::BaseField {
        use crate::Zero;
        Self::BaseField::zero()
    }

    #[inline(always)]
    fn has_glv() -> bool {
        true
    }

    #[inline(always)]
    fn glv_endomorphism_in_place(elem: &mut Self::BaseField) {
        *elem *= &<Self as GLVParameters>::OMEGA;
    }

    #[inline]
    fn glv_scalar_decomposition(
        k: <Self::ScalarField as PrimeField>::BigInt,
    ) -> (
        (bool, <Self::ScalarField as PrimeField>::BigInt),
        (bool, <Self::ScalarField as PrimeField>::BigInt),
    ) {
        <Self as GLVParameters>::glv_scalar_decomposition_inner(k)
    }
}

/// G2_GENERATOR_X =
///  6445332910596979336035888152774071626898886139774101364933948236926875073754470830732273879639675437155036544153105017729592600560631678554299562762294743927912429096636156401171909259073181112518725201388196280039960074422214428
#[rustfmt::skip]
pub const G2_GENERATOR_X: Fq = field_new!(Fq, BigInteger768([
    0x3d902a84cd9f4f78,
    0x864e451b8a9c05dd,
    0xc2b3c0d6646c5673,
    0x17a7682def1ecb9d,
    0xbe31a1e0fb768fe3,
    0x4df125e09b92d1a6,
    0x0943fce635b02ee9,
    0xffc8e7ad0605e780,
    0x8165c00a39341e95,
    0x8ccc2ae90a0f094f,
    0x73a8b8cc0ad09e0c,
    0x11027e203edd9f4,
]));

/// G2_GENERATOR_Y =
/// 562923658089539719386922163444547387757586534741080263946953401595155211934630598999300396317104182598044793758153214972605680357108252243146746187917218885078195819486220416605630144001533548163105316661692978285266378674355041
#[rustfmt::skip]
pub const G2_GENERATOR_Y: Fq = field_new!(Fq, BigInteger768([
    0x9a159be4e773f67c,
    0x6b957244aa8f4e6b,
    0xa27b70c9c945a38c,
    0xacb6a09fda11d0ab,
    0x3abbdaa9bb6b1291,
    0xdbdf642af5694c36,
    0xb6360bb9560b369f,
    0xac0bd1e822b8d6da,
    0xfa355d17afe6945f,
    0x8d6a0fc1fbcad35e,
    0x72a63c7874409840,
    0x114976e5b0db280,
]));
