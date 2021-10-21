use crate::field_new;
use crate::{
    biginteger::BigInteger256,
    curves::{
        models::short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        ModelParameters, SWModelParameters,
    },
    fields::tweedle::*,
    Field,
};

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct TweedledeeParameters;

impl ModelParameters for TweedledeeParameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

pub type Affine = GroupAffine<TweedledeeParameters>;
pub type Projective = GroupProjective<TweedledeeParameters>;

impl SWModelParameters for TweedledeeParameters {
    /// COEFF_A = 0
    const COEFF_A: Fq = field_new!(Fq, BigInteger256([0x0, 0x0, 0x0, 0x0]));

    /// COEFF_B = 5
    const COEFF_B: Fq = field_new!(
        Fq,
        BigInteger256([
            0x30aef343ffffffed,
            0xbcb60a132dafff0b,
            0xffffffffffffffff,
            0x3fffffffffffffff
        ])
    );

    /// COFACTOR = 1
    const COFACTOR: &'static [u64] = &[0x1];

    /// COFACTOR_INV = 1
    const COFACTOR_INV: Fr = field_new!(
        Fr,
        BigInteger256([
            0x1c3ed159fffffffd,
            0xf5601c89bb41f2d3,
            0xffffffffffffffff,
            0x3fffffffffffffff
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
pub const G_GENERATOR_X: Fq = field_new!(
    Fq,
    BigInteger256([
        0x7379f083fffffffd,
        0xf5601c89c3d86ba3,
        0xffffffffffffffff,
        0x3fffffffffffffff
    ])
);

/// G1_GENERATOR_Y =
/// 14240188643175251183985684255458419213835105645119662786317263805424119994471
pub const G_GENERATOR_Y: Fq = field_new!(
    Fq,
    BigInteger256([
        0x1e28b41c4fc25056,
        0x544abc778a8b7ce5,
        0xe133d9afa567f32b,
        0x37705e17172ff461
    ])
);
