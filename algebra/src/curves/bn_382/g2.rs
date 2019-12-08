use crate::{
    biginteger::BigInteger384 as BigInteger,
    curves::models::{ModelParameters, SWModelParameters},
    field_new,
    fields::{
        bn_382::{Fp, Fq, Fq2},
        Field,
    },
};

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct Bn_382G2Parameters;

impl ModelParameters for Bn_382G2Parameters {
    type BaseField = Fq2;
    type ScalarField = Fp;
}

impl SWModelParameters for Bn_382G2Parameters {
    // y^2 = x^3 + 14 / ((2*sqrt(7))^5) = x^3 + sqrt(7)

    /// COEFF_A = [0, 0]
    const COEFF_A: Fq2 = field_new!(
        Fq2,
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0])),
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0])),
    );

    // (14 / (2*sqrt(7))^5
    // == 0 + sqrt7 *
    // 671741409037656549287655731709824109253980562797465531047568917158473772953357661607607074171171789249425365013734
    /// COEFF_B = [4, 4]
    const COEFF_B: Fq2 = field_new!(
        Fq2,
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0])),
        field_new!(
            Fq,
            BigInteger([
                0x1f58d0fac687d635,
                0x4924924a4690f812,
                0x7065db87ffa97412,
                0x63c12e0f38c8a8c2,
                0x849567165657c4df,
                0x1b860998b4eca50d
            ])
        ),
    );

    /// COFACTOR =
    /// 5543634365110765627805495722742127385843376434033820803594923240297849259333798279370015902197046673895926135783425
    const COFACTOR: &'static [u64] = &[
        1,
        6443243544,
        16149412065705181664,
        13009269933821593857,
        8165092549055070088,
        2595350192619816846,
    ];

    /// COFACTOR_INV = COFACTOR^{-1} mod r
    /// 2771817182555382813902747861510987483969025170954974114210393735761260242342976119835350290882788898215364334190594
    const COFACTOR_INV: Fp = field_new!(
        Fp,
        BigInteger([
            0x30f31f7dc6e0f46c,
            0x77ea8678846f2b23,
            0x1391dd7fdcdaa32,
            0xf3afdcd2c5fbcc02,
            0xf9573d5af51e891e,
            0x1580839fc48bbede
        ])
    );

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G2_GENERATOR_X, G2_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }
}

pub const G2_GENERATOR_X: Fq2 = field_new!(Fq2, G2_GENERATOR_X_C0, G2_GENERATOR_X_C1);
pub const G2_GENERATOR_Y: Fq2 = field_new!(Fq2, G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1);

// Generator:
// (3519382844713541579002133617775236000337302709092053889907196608497211512910083011998063983635946531824025900302318*sqrt7 + 5091479006341624589567896397635435258574014748076809289641574502625108749078943401554928186045022840715545119724980,
// 4780208203968490754926961189497985186234872265999339695883768193648752722495923801951744797689497641942946290071424*sqrt7 + 3934462613855637686263305666415197064493526818650772586512345121679314757894509046665527945441022114959626478116310)
/// G2_GENERATOR_X_C0 =
/// 5091479006341624589567896397635435258574014748076809289641574502625108749078943401554928186045022840715545119724980
pub const G2_GENERATOR_X_C0: Fq = field_new!(
    Fq,
    BigInteger([
        0xaa14afc8b159ead4,
        0x277e239ba378f9a2,
        0x53af3b32c1685d95,
        0x5ba7b0850c1229e6,
        0x9b26be1e9ea6773c,
        0x10b027422ad1ccfa
    ])
);

/// G2_GENERATOR_X_C1 =
/// 3519382844713541579002133617775236000337302709092053889907196608497211512910083011998063983635946531824025900302318
pub const G2_GENERATOR_X_C1: Fq = field_new!(
    Fq,
    BigInteger([
        0x54a865e3edc90bb9,
        0x7f5128acfea3ee45,
        0xf2a0887d3d43fc7,
        0x638886d3adde77dc,
        0x40357f01221b9804,
        0x1e9051d8c102908
    ])
);

/// G2_GENERATOR_Y_C0 =
/// 3934462613855637686263305666415197064493526818650772586512345121679314757894509046665527945441022114959626478116310
pub const G2_GENERATOR_Y_C0: Fq = field_new!(
    Fq,
    BigInteger([
        0x82c0caabec3e17d3,
        0x8c65a2ca682b24ae,
        0x44686fe645fd5f7d,
        0xfc711a7c8fdb18fa,
        0x7113b7187965708b,
        0xe837725a5a84e1d
    ])
);

/// G2_GENERATOR_Y_C1 =
/// 4780208203968490754926961189497985186234872265999339695883768193648752722495923801951744797689497641942946290071424
pub const G2_GENERATOR_Y_C1: Fq = field_new!(
    Fq,
    BigInteger([
        0xc97b8db7d421302e,
        0x2cfd45fb7e39854c,
        0x31bf077f06393b01,
        0x2e5cbd509a1ab02c,
        0xbbbde4f85ef3a1aa,
        0xe10c91e02074189
    ])
);
