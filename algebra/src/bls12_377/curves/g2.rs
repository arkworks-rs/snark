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

pub type G2Affine = bls12::G2Affine<bls12_377::Parameters>;
pub type G2Projective = bls12::G2Projective<bls12_377::Parameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl ModelParameters for Parameters {
    type BaseField = Fq2;
    type ScalarField = Fr;
}

impl_scalar_mul_kernel!(bls12_377, "bls12_377", g2, G2Projective);

impl GLVParameters for Parameters {
    type WideBigInt = BigInteger512;
    const OMEGA: Self::BaseField = field_new!(
        Fq2,
        field_new!(
            Fq,
            BigInteger384([
                3203870859294639911,
                276961138506029237,
                9479726329337356593,
                13645541738420943632,
                7584832609311778094,
                101110569012358506
            ])
        ),
        field_new!(Fq, BigInteger384([0, 0, 0, 0, 0, 0]))
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
    /// COEFF_A = [0, 0]
    #[rustfmt::skip]
    const COEFF_A: Fq2 = field_new!(Fq2,
        g1::Parameters::COEFF_A,
        g1::Parameters::COEFF_A,
    );

    // As per https://eprint.iacr.org/2012/072.pdf,
    // this curve has b' = b/i, where b is the COEFF_B of G1, and x^6 -i is
    // the irreducible poly used to extend from Fp2 to Fp12.
    // In our case, i = u (App A.3, T_6).
    /// COEFF_B = [0,
    /// 155198655607781456406391640216936120121836107652948796323930557600032281009004493664981332883744016074664192874906]
    #[rustfmt::skip]
    const COEFF_B: Fq2 = field_new!(Fq2,
        field_new!(Fq, BigInteger384([0, 0, 0, 0, 0, 0])),
        field_new!(Fq, BigInteger384([
            9255502405446297221,
            10229180150694123945,
            9215585410771530959,
            13357015519562362907,
            5437107869987383107,
            16259554076827459,
        ])),
    );

    /// COFACTOR =
    /// 7923214915284317143930293550643874566881017850177945424769256759165301436616933228209277966774092486467289478618404761412630691835764674559376407658497
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        0x0000000000000001,
        0x452217cc90000000,
        0xa0f3622fba094800,
        0xd693e8c36676bd09,
        0x8c505634fae2e189,
        0xfbb36b00e1dcc40c,
        0xddd88d99a6f6a829,
        0x26ba558ae9562a,
    ];

    /// COFACTOR_INV = COFACTOR^{-1} mod r
    /// = 6764900296503390671038341982857278410319949526107311149686707033187604810669
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger256([
        15499857013495546999,
        4613531467548868169,
        14546778081091178013,
        549402535258503313,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G2_GENERATOR_X, G2_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }

    impl_scalar_mul_parameters!(G2Projective);
    impl_glv_for_sw!();
}

#[rustfmt::skip]
pub const G2_GENERATOR_X: Fq2 = field_new!(Fq2, G2_GENERATOR_X_C0, G2_GENERATOR_X_C1);
#[rustfmt::skip]
pub const G2_GENERATOR_Y: Fq2 = field_new!(Fq2, G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1);

/// G2_GENERATOR_X_C0 =
/// 233578398248691099356572568220835526895379068987715365179118596935057653620464273615301663571204657964920925606294
#[rustfmt::skip]
pub const G2_GENERATOR_X_C0: Fq = field_new!(Fq, BigInteger384([
    0x68904082f268725b,
    0x668f2ea74f45328b,
    0xebca7a65802be84f,
    0x1e1850f4c1ada3e6,
    0x830dc22d588ef1e9,
    0x1862a81767c0982,
]));

/// G2_GENERATOR_X_C1 =
/// 140913150380207355837477652521042157274541796891053068589147167627541651775299824604154852141315666357241556069118
#[rustfmt::skip]
pub const G2_GENERATOR_X_C1: Fq = field_new!(Fq, BigInteger384([
    0x5f02a915c91c7f39,
    0xf8c553ba388da2a7,
    0xd51a416dbd198850,
    0xe943c6f38ae3073a,
    0xffe24aa8259a4981,
    0x11853391e73dfdd,
]));

/// G2_GENERATOR_Y_C0 =
/// 63160294768292073209381361943935198908131692476676907196754037919244929611450776219210369229519898517858833747423
#[rustfmt::skip]
pub const G2_GENERATOR_Y_C0: Fq = field_new!(Fq, BigInteger384([
    0xd5b19b897881430f,
    0x5be9118a5b371ed,
    0x6063f91f86c131ee,
    0x3244a61be8f4ec19,
    0xa02e425b9f9a3a12,
    0x18af8c04f3360d2,
]));

/// G2_GENERATOR_Y_C1 =
/// 149157405641012693445398062341192467754805999074082136895788947234480009303640899064710353187729182149407503257491
#[rustfmt::skip]
pub const G2_GENERATOR_Y_C1: Fq = field_new!(Fq, BigInteger384([
    0x57601ac71a5b96f5,
    0xe99acc1714f2440e,
    0x2339612f10118ea9,
    0x8321e68a3b1cd722,
    0x2b543b050cc74917,
    0x590182b396c112,
]));
