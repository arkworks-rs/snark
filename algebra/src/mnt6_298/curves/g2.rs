use crate::mnt6_298::{self, g1, Fq, Fq3, Fr};
use algebra_core::{
    biginteger::BigInteger320,
    curves::{
        mnt6,
        mnt6::MNT6Parameters,
        models::{ModelParameters, SWModelParameters},
    },
    field_new,
};

pub type G2Affine = mnt6::G2Affine<mnt6_298::Parameters>;
pub type G2Projective = mnt6::G2Projective<mnt6_298::Parameters>;
pub type G2Prepared = mnt6::G2Prepared<mnt6_298::Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl ModelParameters for Parameters {
    type BaseField = Fq3;
    type ScalarField = Fr;
}

/// MUL_BY_A_C0 = NONRESIDUE * COEFF_A
    #[rustfmt::skip]
pub const MUL_BY_A_C0: Fq = field_new!(Fq, BigInteger320([
    0xa07b458bf1496fab,
    0xde8254e6541f9fb4,
    0xb1b5cc7bf859c3ea,
    0xf83c4d58364645a9,
    0x30a29b55fa2,
]));

/// MUL_BY_A_C1 = NONRESIDUE * COEFF_A
    #[rustfmt::skip]
pub const MUL_BY_A_C1: Fq = field_new!(Fq, BigInteger320([
    0xa07b458bf1496fab,
    0xde8254e6541f9fb4,
    0xb1b5cc7bf859c3ea,
    0xf83c4d58364645a9,
    0x30a29b55fa2,
]));

/// MUL_BY_A_C2 = COEFF_A
pub const MUL_BY_A_C2: Fq = g1::Parameters::COEFF_A;

impl SWModelParameters for Parameters {
    const COEFF_A: Fq3 = mnt6_298::Parameters::TWIST_COEFF_A;
    #[rustfmt::skip]
    const COEFF_B: Fq3 = field_new!(Fq3,
        field_new!(Fq, BigInteger320([
            0x79a4c2cea3c84026,
            0x4b50cad0f3233baa,
            0x9ded82770e7a4410,
            0x5ade8b105838b95d,
            0xe4036e0a3a,
        ])),
        field_new!(Fq, BigInteger320([0, 0, 0, 0, 0])),
        field_new!(Fq, BigInteger320([0, 0, 0, 0, 0])),
    );

    /// COFACTOR =
    /// 226502022472576270196498690498308461791828762732602586162207535351960270082712694977333372361549082214519252261735048131889018501404377856786623430385820659037970876666767495659520
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        15308190245346869248,
        10669098443577192943,
        4561413759929581409,
        3680089780298582849,
        17336300687782721465,
        10745756320947240891,
        17479264233688728128,
        16828697388537672097,
        4184034152442024798,
        915787,
    ];

    /// COFACTOR^(-1) mod r =
    /// 79320381028210220958891541608841408590854146655427655872973753568875979721417185067925504
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger320([
        5837598184463018016,
        7845868194417674836,
        12170332588914158076,
        6950611683754678431,
        102280178745,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G2_GENERATOR_X, G2_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(elt: &Fq3) -> Fq3 {
        field_new!(
            Fq3,
            MUL_BY_A_C0 * &elt.c1,
            MUL_BY_A_C1 * &elt.c2,
            MUL_BY_A_C2 * &elt.c0,
        )
    }
}

const G2_GENERATOR_X: Fq3 =
    field_new!(Fq3, G2_GENERATOR_X_C0, G2_GENERATOR_X_C1, G2_GENERATOR_X_C2);
const G2_GENERATOR_Y: Fq3 =
    field_new!(Fq3, G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1, G2_GENERATOR_Y_C2);

#[rustfmt::skip]
pub const G2_GENERATOR_X_C0: Fq = field_new!(Fq, BigInteger320([
    0x15ca12fc5d551ea7,
    0x9e0b2b2b2bb8b979,
    0xe6e66283ad5a786a,
    0x46ba0aedcc383c07,
    0x243853463ed,
]));

#[rustfmt::skip]
pub const G2_GENERATOR_X_C1: Fq = field_new!(Fq, BigInteger320([
    0x2c0e3dd7be176130,
    0x27a15d879495904b,
    0x6f1f0d2dd1502a82,
    0x9782ee3c70834da,
    0x2c28bb71862,
]));

#[rustfmt::skip]
pub const G2_GENERATOR_X_C2: Fq = field_new!(Fq, BigInteger320([
    0xf3e5f4eb9631e1f1,
    0x657801e80c50778,
    0x2d2abb128fee90f3,
    0x72e58e4c3aa3598c,
    0x100b8026b9d,
]));

#[rustfmt::skip]
pub const G2_GENERATOR_Y_C0: Fq = field_new!(Fq, BigInteger320([
    0xb1cddd6c64a67c5f,
    0xa01e90d89aa5d2ba,
    0x39e9a733be49ed1,
    0x9438f46f63d3264f,
    0x12cc928ef10,
]));

#[rustfmt::skip]
pub const G2_GENERATOR_Y_C1: Fq = field_new!(Fq, BigInteger320([
    0xa1529b7265ad4be7,
    0x21c5e827cf309306,
    0x9b3d647bd8c70b22,
    0x42835bf373e4b213,
    0xd3c77c9ff9,
]));

#[rustfmt::skip]
pub const G2_GENERATOR_Y_C2: Fq = field_new!(Fq, BigInteger320([
    0x610557ec4b58b8df,
    0x51a23865b52045f1,
    0x9dcfd915a09da608,
    0x6d65c95f69adb700,
    0x2d3c3d195a1,
]));
