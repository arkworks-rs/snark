use crate::{
    biginteger::BigInteger384 as BigInteger,
    curves::models::{ModelParameters, SWModelParameters},
    Field, field_new,
    fields::bn_382::*
};

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct Bn382G2Parameters;

impl ModelParameters for Bn382G2Parameters {
    type BaseField = Fq2;
    type ScalarField = Fr;
}

impl SWModelParameters for Bn382G2Parameters {
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
    const COEFF_B: Fq2 = field_new!(
        Fq2,
        field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0])),
        field_new!(
            Fq,
            BigInteger([
                0xaaaaaaaaaaaaaaa6,
                0xaaaaaaa3aa723a3a,
                0xaa39c8039c9ca3aa,
                0xb57b03a1ed2c80a3,
                0x99dee612c7ca822c,
                0x2957f80ad62dcbc
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
    const COFACTOR_INV: Fr = field_new!(
        Fr,
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
        0x2761a83ee6ccb6c6,
        0x681d8b2b656ce886,
        0x3540fb52bab89b4,
        0x81e6427c08680553,
        0xa9ccf8c26dcf6e1,
        0x1851476ea8077fc6
    ])
);

/// G2_GENERATOR_X_C1 =
/// 3519382844713541579002133617775236000337302709092053889907196608497211512910083011998063983635946531824025900302318
pub const G2_GENERATOR_X_C1: Fq = field_new!(
    Fq,
    BigInteger([
        0x14498e69e5b53113,
        0xee8cd774d8d88e77,
        0xc6f3b5ce2ace1aef,
        0x3502bb8b846944a9,
        0xc95e755dd7927cae,
        0x7c0beebd73ab8f5
    ])
);

/// G2_GENERATOR_Y_C0 =
/// 3934462613855637686263305666415197064493526818650772586512345121679314757894509046665527945441022114959626478116310
pub const G2_GENERATOR_Y_C0: Fq = field_new!(
    Fq,
    BigInteger([
        0xe8028a161e1bbc9a,
        0x266bd5d118d75d9b,
        0xacc76640f1e4baa9,
        0xa70bd81b6be756f4,
        0x5161ebf3eef9a86e,
        0x657061a71f10b07
    ])
);

/// G2_GENERATOR_Y_C1 =
/// 4780208203968490754926961189497985186234872265999339695883768193648752722495923801951744797689497641942946290071424
pub const G2_GENERATOR_Y_C1: Fq = field_new!(
    Fq,
    BigInteger([
        0xc981004fff54161b,
        0x12540e2eb9d55972,
        0xf89981d85302a29a,
        0xb69dcab62945321d,
        0x841115a42fc75e00,
        0x65c141e4b7455c1
    ])
);
