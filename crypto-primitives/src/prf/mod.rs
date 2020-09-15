use algebra_core::bytes::{FromBytes, ToBytes};
use core::{fmt::Debug, hash::Hash};

use crate::CryptoError;

#[cfg(feature = "r1cs")]
pub mod constraints;
#[cfg(feature = "r1cs")]
pub use constraints::*;

pub mod blake2s;
pub use self::blake2s::*;

pub trait PRF {
    type Input: FromBytes + Default;
    type Output: ToBytes + Eq + Clone + Debug + Default + Hash;
    type Seed: FromBytes + ToBytes + Clone + Default + Debug;

    fn evaluate(seed: &Self::Seed, input: &Self::Input) -> Result<Self::Output, CryptoError>;
}
