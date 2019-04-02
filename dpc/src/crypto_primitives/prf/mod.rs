use algebra::bytes::{FromBytes, ToBytes};
use std::{fmt::Debug, hash::Hash};

use crate::crypto_primitives::CryptoError;

pub mod blake2s;
pub use self::blake2s::*;

pub trait PRF {
    type Input: FromBytes + Default;
    type Output: ToBytes + Eq + Clone + Default + Hash;
    type Seed: FromBytes + ToBytes + Clone + Default + Debug;

    fn evaluate(seed: &Self::Seed, input: &Self::Input) -> Result<Self::Output, CryptoError>;
}
