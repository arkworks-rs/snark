use rand::{Rand, Rng};
use std::{fmt::Debug, hash::Hash};

use algebra::bytes::ToBytes;

pub mod blake2s;
pub mod injective_map;
pub mod pedersen;

use crate::Error;

pub trait CommitmentScheme {
    type Output: ToBytes + Clone + Default + Eq + Hash + Debug;
    type Parameters: Clone;
    type Randomness: Clone + ToBytes + Default + Eq + Rand + Debug;

    fn setup<R: Rng>(r: &mut R) -> Result<Self::Parameters, Error>;
    fn commit(
        parameters: &Self::Parameters,
        input: &[u8],
        r: &Self::Randomness,
    ) -> Result<Self::Output, Error>;
}
