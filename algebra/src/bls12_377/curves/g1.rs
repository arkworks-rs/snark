use algebra_core::{
    biginteger::{BigInteger256, BigInteger384, BigInteger512},
    curves::{
        bls12,
        models::{ModelParameters, SWModelParameters},
        GLVParameters,
    },
    field_new, impl_glv_for_sw, impl_scalar_mul_kernel, impl_scalar_mul_parameters, PrimeField,
    Zero,
};

use crate::{bls12_377, bls12_377::*};

pub type G1Affine = bls12::G1Affine<bls12_377::Parameters>;
pub type G1Projective = bls12::G1Projective<bls12_377::Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl ModelParameters for Parameters {
    type BaseField = Fq;
    type ScalarField = Fr;
}

impl_scalar_mul_kernel!(bls12_377, "bls12_377", g1, G1Projective);

impl GLVParameters for Parameters {
    type WideBigInt = BigInteger512;
    const OMEGA: Self::BaseField = field_new!(
        Fq,
        BigInteger384([
            15766275933608376691,
            15635974902606112666,
            1934946774703877852,
            18129354943882397960,
            15437979634065614942,
            101285514078273488
        ])
    );
    const LAMBDA: Self::ScalarField = field_new!(
        Fr,
        BigInteger256([
            12574070832645531618,
            10005695704657941814,
            1564543351912391449,
            657300228442948690
        ])
    );
    /// |round(B1 * R / n)|
    const Q2: <Self::ScalarField as PrimeField>::BigInt =
        BigInteger256([9183663392111466540, 12968021215939883360, 3, 0]);
    const B1: <Self::ScalarField as PrimeField>::BigInt =
        BigInteger256([725501752471715841, 4981570305181876225, 0, 0]);
    const B1_IS_NEG: bool = false;
    /// |round(B2 * R / n)|
    const Q1: <Self::ScalarField as PrimeField>::BigInt = BigInteger256([13, 0, 0, 0]);
    const B2: <Self::ScalarField as PrimeField>::BigInt = BigInteger256([1, 0, 0, 0]);
    const R_BITS: u32 = 256;
}

impl SWModelParameters for Parameters {
    /// COEFF_A = 0
    const COEFF_A: Fq = field_new!(Fq, BigInteger384([0x0, 0x0, 0x0, 0x0, 0x0, 0x0]));

    /// COEFF_B = 1
    #[rustfmt::skip]
    const COEFF_B: Fq = field_new!(Fq, BigInteger384([
        0x2cdffffffffff68,
        0x51409f837fffffb1,
        0x9f7db3a98a7d3ff2,
        0x7b4e97b76e7c6305,
        0x4cf495bf803c84e8,
        0x8d6661e2fdf49a,
    ]));

    /// COFACTOR = (x - 1)^2 / 3  = 30631250834960419227450344600217059328
    const COFACTOR: &'static [u64] = &[0x0, 0x170b5d4430000000];

    /// COFACTOR_INV = COFACTOR^{-1} mod r
    /// = 5285428838741532253824584287042945485047145357130994810877
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger256([
        2013239619100046060,
        4201184776506987597,
        2526766393982337036,
        1114629510922847535,
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
/// 81937999373150964239938255573465948239988671502647976594219695644855304257327692006745978603320413799295628339695
#[rustfmt::skip]
pub const G1_GENERATOR_X: Fq = field_new!(Fq, BigInteger384([
    0x260f33b9772451f4,
    0xc54dd773169d5658,
    0x5c1551c469a510dd,
    0x761662e4425e1698,
    0xc97d78cc6f065272,
    0xa41206b361fd4d,
]));

/// G1_GENERATOR_Y =
/// 241266749859715473739788878240585681733927191168601896383759122102112907357779751001206799952863815012735208165030
#[rustfmt::skip]
pub const G1_GENERATOR_Y: Fq = field_new!(Fq, BigInteger384([
    0x8193961fb8cb81f3,
    0x638d4c5f44adb8,
    0xfafaf3dad4daf54a,
    0xc27849e2d655cd18,
    0x2ec3ddb401d52814,
    0x7da93326303c71,
]));
