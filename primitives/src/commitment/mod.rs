use algebra::bytes::ToBytes;
use algebra::UniformRand;
use rand::Rng;
use std::{fmt::Debug, hash::Hash};

use serde::{Deserialize, Serialize};

pub mod blake2s;
pub mod injective_map;
pub mod pedersen;

use crate::Error;

pub trait CommitmentScheme {
    type Output: ToBytes + Serialize + for<'a> Deserialize<'a> + Clone + Default + Eq + Hash + Debug;
    type Parameters: Clone + Serialize + for<'a> Deserialize<'a>;
    type Randomness: Clone
        + ToBytes
        + Serialize
        + for<'a> Deserialize<'a>
        + Default
        + Eq
        + UniformRand
        + Debug;

    fn setup<R: Rng>(r: &mut R) -> Result<Self::Parameters, Error>;

    fn commit(
        parameters: &Self::Parameters,
        input: &[u8],
        r: &Self::Randomness,
    ) -> Result<Self::Output, Error>;
}
