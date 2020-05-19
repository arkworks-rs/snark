use crate::{
    bw6_761::*,
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
    const ATE_LOOP_COUNT: &'static [u64] = &[
        0x467a80000000000f,
        0x70b5d44300000007,
        0x58490fb409869401,
        0xb55fc0d440cb48f0,
        0x1ab2f9cb6145aeec,
        0x15d8f58f3501dbec,
    ];
    const ATE_LOOP_COUNT_IS_NEGATIVE: bool = false;
    const FINAL_EXPONENT_LAST_CHUNK_1: BigInteger = BigInteger([
        0x3de5800000000089,
        0x832ba4061000003b,
        0xc61c554757551c0c,
        0xc856a0853c9db94c,
        0x2c77d5ac34cb12ef,
        0xad1972339049ce76,
        0x0,
        0x0,
        0x0,
        0x0,
        0x0,
        0x0,
    ]);
    const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = false;
    const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger = BigInteger([
        0x6f9440000000008c,
        0x1aff40fcf0000082,
        0x9521646d73808c51,
        0x3ba806d298c79fc5,
        0xb521a3d9309c6dd0,
        0x824cd7cfb1e8685a,
        0xa7f6ef02c228c497,
        0xa311dc0a5ef6ff10,
        0x96a147eaf584608d,
        0x828e2c6f9f4f1494,
        0x68f6427062e1b0b,
        0x0,
    ]);
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
