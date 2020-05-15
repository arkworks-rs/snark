use crate::bw6_761::*;
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
    const X: &'static [u64] = &[0x8508c00000000001];
    /// `x` is positive.
    const X_IS_NEGATIVE: bool = false;
    const ATE_LOOP_COUNT_1: &'static [u64] = &[0x8508c00000000002];
    const ATE_LOOP_COUNT_1_IS_NEGATIVE: bool = false;
    const ATE_LOOP_COUNT_2: &'static [u64] = &[
        0xffffffffffffffff,
        0x8a442f991fffffff,
        0x23ed1347970dec00,
    ];
    const ATE_LOOP_COUNT_2_IS_NEGATIVE: bool = false;
    const TWIST_TYPE: TwistType = TwistType::D;
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
