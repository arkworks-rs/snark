use algebra_core::{
    biginteger::{BigInteger256, BigInteger512},
    curves::{
        bn,
        models::{ModelParameters, SWModelParameters},
    },
    field_new, impl_glv_for_sw, impl_scalar_mul_kernel, impl_scalar_mul_parameters, GLVParameters,
    PrimeField, Zero,
};

use crate::{bn254, bn254::*};

pub type G1Affine = bn::G1Affine<bn254::Parameters>;
pub type G1Projective = bn::G1Projective<bn254::Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl ModelParameters for Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl_scalar_mul_kernel!(bn254, "bn254", g1, G1Projective);

impl GLVParameters for Parameters {
    type WideBigInt = BigInteger512;
    const OMEGA: Self::BaseField = field_new!(
        Fq,
        BigInteger256([
            3697675806616062876,
            9065277094688085689,
            6918009208039626314,
            2775033306905974752
        ])
    );
    const LAMBDA: Self::ScalarField = field_new!(
        Fr,
        BigInteger256([
            244305545194690131,
            8351807910065594880,
            14266533074055306532,
            404339206190769364
        ])
    );
    /// |round(B1 * R / n)|
    /// |round(B1 * R / n)|
    const Q2: <Self::ScalarField as PrimeField>::BigInt =
        BigInteger256([6023842690951505253, 5534624963584316114, 2, 0]);
    const B1: <Self::ScalarField as PrimeField>::BigInt =
        BigInteger256([857057580304901387, 8020209761171036669, 0, 0]);
    const B1_IS_NEG: bool = false;
    /// |round(B2 * R / n)|
    const Q1: <Self::ScalarField as PrimeField>::BigInt =
        BigInteger256([15644699364383830999, 2, 0, 0]);
    const B2: <Self::ScalarField as PrimeField>::BigInt =
        BigInteger256([9931322734385697763, 0, 0, 0]);
    const R_BITS: u32 = 256;
}

impl SWModelParameters for Parameters {
    /// COEFF_A = 0
    const COEFF_A: Fq = field_new!(Fq, BigInteger256([0x0, 0x0, 0x0, 0x0]));

    /// COEFF_B = 3
    #[rustfmt::skip]
    const COEFF_B: Fq = field_new!(Fq, BigInteger256([
        0x7a17caa950ad28d7,
        0x1f6ac17ae15521b9,
        0x334bea4e696bd284,
        0x2a1f6744ce179d8e,
    ]));

    /// COFACTOR = 1
    const COFACTOR: &'static [u64] = &[0x1];

    /// COFACTOR_INV = COFACTOR^{-1} mod r = 1
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger256([
        0xac96341c4ffffffb,
        0x36fc76959f60cd29,
        0x666ea36f7879462e,
        0xe0a77c19a07df2f,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G1_GENERATOR_X, G1_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }

    impl_scalar_mul_parameters!(G1Projective);
    impl_glv_for_sw!();
}

/// G1_GENERATOR_X =
/// 1
#[rustfmt::skip]
pub const G1_GENERATOR_X: Fq = field_new!(Fq, BigInteger256([
    0xd35d438dc58f0d9d,
    0x0a78eb28f5c70b3d,
    0x666ea36f7879462c,
    0x0e0a77c19a07df2f,
]));

/// G1_GENERATOR_Y =
/// 2
#[rustfmt::skip]
pub const G1_GENERATOR_Y: Fq = field_new!(Fq, BigInteger256([
    0xa6ba871b8b1e1b3a,
    0x14f1d651eb8e167b,
    0xccdd46def0f28c58,
    0x1c14ef83340fbe5e,
]));
