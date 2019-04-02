use crate::{
    curves::bls12::{Bls12, Bls12Parameters, TwistType},
    fields::bls12_381::{Fq, Fq12Parameters, Fq2Parameters, Fq6Parameters},
};

pub mod g1;
pub mod g2;
#[cfg(test)]
mod tests;

use self::{g1::Bls12_381G1Parameters, g2::Bls12_381G2Parameters};

pub use self::{
    g1::{G1Affine, G1Projective},
    g2::{G2Affine, G2Projective},
};

pub type Bls12_381 = Bls12<Bls12_381Parameters>;

pub struct Bls12_381Parameters;

impl Bls12Parameters for Bls12_381Parameters {
    const X: &'static [u64] = &[0xd201000000010000];
    const X_IS_NEGATIVE: bool = true;
    const TWIST_TYPE: TwistType = TwistType::M;
    type Fp = Fq;
    type Fp2Params = Fq2Parameters;
    type Fp6Params = Fq6Parameters;
    type Fp12Params = Fq12Parameters;
    type G1Parameters = Bls12_381G1Parameters;
    type G2Parameters = Bls12_381G2Parameters;
}
