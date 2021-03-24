use crate::{
    biginteger::BigInteger384,
    fields::bn_382::*,
    curves::{
        models::{ModelParameters, SWModelParameters},
    },
    Field, field_new,
};

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct Bn382G1Parameters;

impl ModelParameters for Bn382G1Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl SWModelParameters for Bn382G1Parameters {
    /// COEFF_A = 0
    const COEFF_A: Fq = field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0]));

    /// COEFF_B = 14
    const COEFF_B: Fq = field_new!(
        Fq,
        BigInteger384([
            0xffffffffffffff9d,
            0xffffff6b7b52aeb7,
            0x76a537ba55d66b7f,
            0x2e8d16344c0b846b,
            0x2df8a320b2feee22,
            0x123eec4e5e4393ea
        ])
    );

    /// COFACTOR = 1
    const COFACTOR: &'static [u64] = &[0x1];

    /// COFACTOR_INV = 1
    const COFACTOR_INV: Fr = field_new!(
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

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G1_GENERATOR_X, G1_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }
}

/// G1_GENERATOR_X =
/// 1
pub const G1_GENERATOR_X: Fq = field_new!(
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

/// G1_GENERATOR_Y =
/// 93360544046129830094757569027791679210844519762232758194920967606984287664392872848607365449491441272860487554919
pub const G1_GENERATOR_Y: Fq = field_new!(
    Fq,
    BigInteger384([
        0x9dfafa6eb6e5986a,
        0x320ae00a19eea8eb,
        0x740e245a3411fca8,
        0x7ad3304e255f5799,
        0x310b3464a5ff421d,
        0x12713e4c3440dde
    ])
);
