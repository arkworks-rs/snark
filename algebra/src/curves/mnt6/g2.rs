use crate::field_new;
use crate::{
    biginteger::BigInteger320,
    bytes::ToBytes,
    curves::{
        mnt6::{g1::MNT6G1Parameters, G1Affine, MNT6, TWIST_COEFF_A},
        models::{ModelParameters, SWModelParameters},
        short_weierstrass_projective::{GroupAffine, GroupProjective},
        AffineCurve, PairingCurve, PairingEngine,
    },
    fields::mnt6::{Fq, Fq3, Fq6, Fr},
};
use std::io::{Result as IoResult, Write};

pub type G2Affine = GroupAffine<MNT6G2Parameters>;
pub type G2Projective = GroupProjective<MNT6G2Parameters>;

impl PairingCurve for G2Affine {
    type Engine = MNT6;
    type Prepared = G2Prepared;
    type PairWith = G1Affine;
    type PairingResult = Fq6;

    #[inline(always)]
    fn prepare(&self) -> Self::Prepared {
        Self::Prepared::from_affine(self)
    }

    #[inline(always)]
    fn pairing_with(&self, other: &Self::PairWith) -> Self::PairingResult {
        MNT6::pairing(*other, *self)
    }
}

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct MNT6G2Parameters;

impl ModelParameters for MNT6G2Parameters {
    type BaseField = Fq3;
    type ScalarField = Fr;
}

/// MUL_BY_A_C0 = NONRESIDUE * COEFF_A
pub const MUL_BY_A_C0: Fq = field_new!(Fq, BigInteger320([
    0xa07b458bf1496fab,
    0xde8254e6541f9fb4,
    0xb1b5cc7bf859c3ea,
    0xf83c4d58364645a9,
    0x30a29b55fa2,
]));

/// MUL_BY_A_C1 = NONRESIDUE * COEFF_A
pub const MUL_BY_A_C1: Fq = field_new!(Fq, BigInteger320([
    0xa07b458bf1496fab,
    0xde8254e6541f9fb4,
    0xb1b5cc7bf859c3ea,
    0xf83c4d58364645a9,
    0x30a29b55fa2,
]));

/// MUL_BY_A_C2 = COEFF_A
pub const MUL_BY_A_C2: Fq = MNT6G1Parameters::COEFF_A;

impl SWModelParameters for MNT6G2Parameters {
    const COEFF_A: Fq3 = TWIST_COEFF_A;
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
        field_new!(Fq3, 
            MUL_BY_A_C0 * &elt.c1,
            MUL_BY_A_C1 * &elt.c2,
            MUL_BY_A_C2 * &elt.c0,
        )
    }
}

const G2_GENERATOR_X: Fq3 = field_new!(Fq3, G2_GENERATOR_X_C0, G2_GENERATOR_X_C1, G2_GENERATOR_X_C2);
const G2_GENERATOR_Y: Fq3 = field_new!(Fq3, G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1, G2_GENERATOR_Y_C2);

pub const G2_GENERATOR_X_C0: Fq = field_new!(Fq, BigInteger320([
    0x15ca12fc5d551ea7,
    0x9e0b2b2b2bb8b979,
    0xe6e66283ad5a786a,
    0x46ba0aedcc383c07,
    0x243853463ed,
]));

pub const G2_GENERATOR_X_C1: Fq = field_new!(Fq, BigInteger320([
    0x2c0e3dd7be176130,
    0x27a15d879495904b,
    0x6f1f0d2dd1502a82,
    0x9782ee3c70834da,
    0x2c28bb71862,
]));

pub const G2_GENERATOR_X_C2: Fq = field_new!(Fq, BigInteger320([
    0xf3e5f4eb9631e1f1,
    0x657801e80c50778,
    0x2d2abb128fee90f3,
    0x72e58e4c3aa3598c,
    0x100b8026b9d,
]));

pub const G2_GENERATOR_Y_C0: Fq = field_new!(Fq, BigInteger320([
    0xb1cddd6c64a67c5f,
    0xa01e90d89aa5d2ba,
    0x39e9a733be49ed1,
    0x9438f46f63d3264f,
    0x12cc928ef10,
]));

pub const G2_GENERATOR_Y_C1: Fq = field_new!(Fq, BigInteger320([
    0xa1529b7265ad4be7,
    0x21c5e827cf309306,
    0x9b3d647bd8c70b22,
    0x42835bf373e4b213,
    0xd3c77c9ff9,
]));

pub const G2_GENERATOR_Y_C2: Fq = field_new!(Fq, BigInteger320([
    0x610557ec4b58b8df,
    0x51a23865b52045f1,
    0x9dcfd915a09da608,
    0x6d65c95f69adb700,
    0x2d3c3d195a1,
]));

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct G2Prepared {
    pub x:                     Fq3,
    pub y:                     Fq3,
    pub x_over_twist:          Fq3,
    pub y_over_twist:          Fq3,
    pub double_coefficients:   Vec<AteDoubleCoefficients>,
    pub addition_coefficients: Vec<AteAdditionCoefficients>,
}

impl ToBytes for G2Prepared {
    fn write<W: Write>(&self, _writer: W) -> IoResult<()> {
        unimplemented!()
    }
}

impl G2Prepared {
    pub fn from_affine(point: &G2Affine) -> Self {
        MNT6::ate_precompute_g2(&point.into_projective())
    }
}

impl Default for G2Prepared {
    fn default() -> Self {
        Self::from_affine(&G2Affine::prime_subgroup_generator())
    }
}

pub(super) struct G2ProjectiveExtended {
    pub(crate) x: Fq3,
    pub(crate) y: Fq3,
    pub(crate) z: Fq3,
    pub(crate) t: Fq3,
}

#[derive(Eq, PartialEq, Copy, Clone, Debug)]
pub struct AteDoubleCoefficients {
    pub(crate) c_h:  Fq3,
    pub(crate) c_4c: Fq3,
    pub(crate) c_j:  Fq3,
    pub(crate) c_l:  Fq3,
}

#[derive(Eq, PartialEq, Copy, Clone, Debug)]
pub struct AteAdditionCoefficients {
    pub(crate) c_l1: Fq3,
    pub(crate) c_rz: Fq3,
}
