use algebra_core::{
    biginteger::BigInteger256,
    curves::models::{ModelParameters, SWModelParameters},
    field_new, Zero,
};

use crate::bn254::{g1, Fq, Fq2, Fr};

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Parameters;

impl ModelParameters for Parameters {
    type BaseField = Fq2;
    type ScalarField = Fr;
}

impl SWModelParameters for Parameters {
    /// COEFF_A = [0, 0]
    #[rustfmt::skip]
    const COEFF_A: Fq2 = field_new!(Fq2,
        g1::Parameters::COEFF_A,
        g1::Parameters::COEFF_A,
    );

    /// COEFF_B = 3/(u+9)
    ///         = (19485874751759354771024239261021720505790618469301721065564631296452457478373, 266929791119991161246907387137283842545076965332900288569378510910307636690)
    #[rustfmt::skip]
    const COEFF_B: Fq2 = field_new!(Fq2,
        field_new!(Fq, BigInteger256([
            0x3bf938e377b802a8,
            0x020b1b273633535d,
            0x26b7edf049755260,
            0x2514c6324384a86d,
        ])),
        field_new!(Fq, BigInteger256([
            0x38e7ecccd1dcff67,
            0x65f0b37d93ce0d3e,
            0xd749d0dd22ac00aa,
            0x0141b9ce4a688d4d,
        ])),
    );

    /// COFACTOR = (36 * X^4) + (36 * X^3) + (30 * X^2) + 6*X + 1
    ///          = 21888242871839275222246405745257275088844257914179612981679871602714643921549
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        0x345f2299c0f9fa8d,
        0x06ceecda572a2489,
        0xb85045b68181585e,
        0x30644e72e131a029,
    ];

    /// COFACTOR_INV = COFACTOR^{-1} mod r
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = field_new!(Fr, BigInteger256([
        0x7fff17d53ff2895e,
        0xd0617390cf7919e5,
        0xb9af426b22d0eb61,
        0x270485e31bd72a4d,
    ]));

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G2_GENERATOR_X, G2_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }
}

#[rustfmt::skip]
pub const G2_GENERATOR_X: Fq2 = field_new!(Fq2, G2_GENERATOR_X_C0, G2_GENERATOR_X_C1);
#[rustfmt::skip]
pub const G2_GENERATOR_Y: Fq2 = field_new!(Fq2, G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1);

/// G2_GENERATOR_X_C0 =
/// 10857046999023057135944570762232829481370756359578518086990519993285655852781
#[rustfmt::skip]
pub const G2_GENERATOR_X_C0: Fq = field_new!(Fq, BigInteger256([
    0x8e83b5d102bc2026,
    0xdceb1935497b0172,
    0xfbb8264797811adf,
    0x19573841af96503b,
]));

/// G2_GENERATOR_X_C1 =
/// 11559732032986387107991004021392285783925812861821192530917403151452391805634
#[rustfmt::skip]
pub const G2_GENERATOR_X_C1: Fq = field_new!(Fq, BigInteger256([
    0xafb4737da84c6140,
    0x6043dd5a5802d8c4,
    0x09e950fc52a02f86,
    0x14fef0833aea7b6b,
]));

/// G2_GENERATOR_Y_C0 =
/// 8495653923123431417604973247489272438418190587263600148770280649306958101930
#[rustfmt::skip]
pub const G2_GENERATOR_Y_C0: Fq = field_new!(Fq, BigInteger256([
    0x619dfa9d886be9f6,
    0xfe7fd297f59e9b78,
    0xff9e1a62231b7dfe,
    0x28fd7eebae9e4206,
]));

/// G2_GENERATOR_Y_C1 =
/// 4082367875863433681332203403145435568316851327593401208105741076214120093531
#[rustfmt::skip]
pub const G2_GENERATOR_Y_C1: Fq = field_new!(Fq, BigInteger256([
    0x64095b56c71856ee,
    0xdc57f922327d3cbb,
    0x55f935be33351076,
    0x0da4a0e693fd6482,
]));
