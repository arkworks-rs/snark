use algebra::bytes::{FromBytes, ToBytes};
use std::{fmt::Debug, hash::Hash};

use crate::CryptoError;

pub mod blake2s;
pub use self::blake2s::*;
use serde::{Deserialize, Serialize};

pub trait PRF {
    type Input: FromBytes + Serialize + for<'a> Deserialize<'a> + Default;
    type Output: ToBytes + Serialize + for<'a> Deserialize<'a> + Eq + Clone + Default + Hash;
    type Seed: FromBytes + ToBytes + Serialize + for<'a> Deserialize<'a> + Clone + Default + Debug;

    fn evaluate(seed: &Self::Seed, input: &Self::Input) -> Result<Self::Output, CryptoError>;
}
