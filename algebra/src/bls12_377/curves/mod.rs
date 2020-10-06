use crate::bls12_377::*;
use algebra_core::curves::bls12::{Bls12, Bls12Parameters, TwistType};

pub mod g1;
pub use self::g1::{G1Affine, G1Projective};

pub mod g2;
pub use self::g2::{G2Affine, G2Projective};

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
