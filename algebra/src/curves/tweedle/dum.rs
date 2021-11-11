use crate::{
    biginteger::BigInteger256,
    curves::{
        models::short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        ModelParameters, SWModelParameters,
    },
    field_new,
    fields::tweedle::*,
    Field,
};

#[derive(Copy, Clone, Default, PartialEq, Eq)]
pub struct TweedledumParameters;

impl ModelParameters for TweedledumParameters {
    type BaseField = Fr;
    type ScalarField = Fq;
}

pub type Affine = GroupAffine<TweedledumParameters>;
pub type Projective = GroupProjective<TweedledumParameters>;

impl SWModelParameters for TweedledumParameters {
    /// COEFF_A = 0
    const COEFF_A: Fr = field_new!(Fr, BigInteger256([0x0, 0x0, 0x0, 0x0]));

    /// COEFF_B = 5
    const COEFF_B: Fr = field_new!(
        Fr,
        BigInteger256([
            0x8388339ffffffed,
            0xbcb60a12f74c5739,
            0xffffffffffffffff,
            0x3fffffffffffffff
        ])
    );

    /// COFACTOR = 1
    const COFACTOR: &'static [u64] = &[0x1];

    /// COFACTOR_INV = 1
    const COFACTOR_INV: Fq = field_new!(
        Fq,
        BigInteger256([
            0x7379f083fffffffd,
            0xf5601c89c3d86ba3,
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
pub const G_GENERATOR_X: Fr = field_new!(
    Fr,
    BigInteger256([
        0x1c3ed159fffffffd,
        0xf5601c89bb41f2d3,
        0xffffffffffffffff,
        0x3fffffffffffffff
    ])
);

/// G1_GENERATOR_Y =
/// 385654983219305453067387443941241858913435815837190103938162313975739315615
pub const G_GENERATOR_Y: Fr = field_new!(
    Fr,
    BigInteger256([
        0x7414a31870fe2315,
        0x5771cccafdb1a2b5,
        0x747fd502e877c849,
        0x3175a51e493b99fc
    ])
);
