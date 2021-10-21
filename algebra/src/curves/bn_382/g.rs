use crate::{
    biginteger::BigInteger384,
    curves::{
        models::short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        ModelParameters, SWModelParameters,
    },
    field_new,
    fields::bn_382::*,
    Field,
};

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct Bn382GParameters;

impl ModelParameters for Bn382GParameters {
    type BaseField = Fr;
    type ScalarField = Fq;
}

pub type Affine = GroupAffine<Bn382GParameters>;
pub type Projective = GroupProjective<Bn382GParameters>;

impl SWModelParameters for Bn382GParameters {
    /// COEFF_A = 0
    const COEFF_A: Fr = field_new!(Fr, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0]));

    /// COEFF_B = 7
    const COEFF_B: Fr = field_new!(
        Fr,
        BigInteger384([
            0xffffffffffffffcf,
            0xffffffb67daf6367,
            0xdc87071c715188df,
            0x718ba6243a5346c8,
            0x4fa46fc531ce56d5,
            0x1b21bac71c8e0dbc
        ])
    );

    /// COFACTOR = 1
    const COFACTOR: &'static [u64] = &[0x1];

    /// COFACTOR_INV = 1
    const COFACTOR_INV: Fq = field_new!(
        Fq,
        BigInteger384([
            0xfffffffffffffff9,
            0xfffffff57fab5757,
            0x7f56ac056aeaf57f,
            0x10388572e3c2c0f5,
            0xe6ce591c2bafc343,
            0x3e03f4104144b1a
        ])
    );

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G_GENERATOR_X, G_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }
}

/// G_GENERATOR_X =
/// 1
pub const G_GENERATOR_X: Fr = field_new!(
    Fr,
    BigInteger384([
        0xfffffffffffffff9,
        0xfffffff57fab5757,
        0x1f8101041030381f,
        0x10388572e3c2c0f8,
        0xe6ce591c2bafc343,
        0x3e03f4104144b1a
    ])
);

/// G1_GENERATOR_Y =
/// 1587713460471950740217388326193312024737041813752165827005856534245539019723616944862168333942330219466268138558982
pub const G_GENERATOR_Y: Fr = field_new!(
    Fr,
    BigInteger384([
        0x7bbbac48dff48e8a,
        0x7f0b69a418192817,
        0x91be699f8043e89b,
        0xb9a47acffcccc09c,
        0xbd7a048e12f9984f,
        0x16e7846105853ac1
    ])
);
