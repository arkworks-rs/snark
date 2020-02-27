use crate::bls12_377::*;
use algebra_core::curves::{
    bls12,
    bls12::{Bls12, Bls12Parameters, TwistType},
};

pub mod g1;
pub mod g2;

#[cfg(test)]
mod tests;

pub struct Parameters;

impl Bls12Parameters for Parameters {
    const X: &'static [u64] = &[0x8508c00000000001];
    /// `x` is positive.
    const X_IS_NEGATIVE: bool = false;
    const TWIST_TYPE: TwistType = TwistType::D;
    type Fp = Fq;
    type Fp2Params = Fq2Parameters;
    type Fp6Params = Fq6Parameters;
    type Fp12Params = Fq12Parameters;
    type G1Parameters = g1::Parameters;
    type G2Parameters = g2::Parameters;
}

pub type Bls12_377 = Bls12<Parameters>;

pub type G1Affine = bls12::G1Affine<Parameters>;
pub type G1Projective = bls12::G1Projective<Parameters>;
pub type G2Affine = bls12::G2Affine<Parameters>;
pub type G2Projective = bls12::G2Projective<Parameters>;
