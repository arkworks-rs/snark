use crate::{
    bw6_761::*,
    field_new,
    biginteger::BigInteger768 as BigInteger,
};

use algebra_core::curves::{
    bw6,
    bw6::{BW6, BW6Parameters, TwistType},
};

pub mod g1;
pub mod g2;

#[cfg(test)]
mod tests;

pub struct Parameters;

impl BW6Parameters for Parameters {
    const X: BigInteger = BigInteger([0x8508c00000000001, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]);
    /// `x` is positive.
    const X_IS_NEGATIVE: bool = false;
    // X+1
    const ATE_LOOP_COUNT_1: &'static [u64] = &[0x8508c00000000002];
    const ATE_LOOP_COUNT_1_IS_NEGATIVE: bool = false;
    // X^3-X^2-X
    const ATE_LOOP_COUNT_2: &'static [u64] = &[
        0xffffffffffffffff,
        0x8a442f991fffffff,
        0x23ed1347970dec00,
    ];
    const ATE_LOOP_COUNT_2_IS_NEGATIVE: bool = false;
    const TWIST_TYPE: TwistType = TwistType::D;
    const TWIST: Fq = field_new!(Fq, BigInteger([
        0x0405ffffffff0baa,
        0xb4b04c6b1fff19ce,
        0x3d32dc8704ff55bc,
        0xb4d5fe641dc8fbe9,
        0xd9d3967c3b297017,
        0x81cccf44a4904817,
        0x4e9b4b7fb95a720b,
        0x46a5cffc8c5e4207,
        0xf6acb100116390f8,
        0x8b0914c7ce22045e,
        0xaf503d773ecb53be,
        0xa3eefde24fd0fb,
    ]));
    type Fp = Fq;
    type Fp3Params = Fq3Parameters;
    type Fp6Params = Fq6Parameters;
    type G1Parameters = g1::Parameters;
    type G2Parameters = g2::Parameters;
}

pub type BW6_761 = BW6<Parameters>;

pub type G1Affine = bw6::G1Affine<Parameters>;
pub type G1Projective = bw6::G1Projective<Parameters>;
pub type G2Affine = bw6::G2Affine<Parameters>;
pub type G2Projective = bw6::G2Projective<Parameters>;
