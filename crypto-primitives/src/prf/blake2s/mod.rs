use blake2::Blake2s as b2s;
use digest::Digest;

use super::PRF;
use crate::CryptoError;

#[cfg(feature = "r1cs")]
pub mod constraints;

#[derive(Clone)]
pub struct Blake2s;

impl PRF for Blake2s {
    type Input = [u8; 32];
    type Output = [u8; 32];
    type Seed = [u8; 32];

    fn evaluate(seed: &Self::Seed, input: &Self::Input) -> Result<Self::Output, CryptoError> {
        let eval_time = start_timer!(|| "Blake2s::Eval");
        let mut h = b2s::new();
        h.input(seed.as_ref());
        h.input(input.as_ref());
        let mut result = [0u8; 32];
        result.copy_from_slice(&h.result());
        end_timer!(eval_time);
        Ok(result)
    }
}
