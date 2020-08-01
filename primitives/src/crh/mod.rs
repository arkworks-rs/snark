use algebra::{
    Field, bytes::ToBytes
};
use rand::Rng;
use std::hash::Hash;

pub mod bowe_hopwood;
pub mod injective_map;
pub mod pedersen;

pub mod poseidon;
pub use self::poseidon::*;

use crate::Error;


pub trait FixedLengthCRH {
    const INPUT_SIZE_BITS: usize;
    type Output: ToBytes + Clone + Eq + Hash + Default;
    type Parameters: Clone + Default;

    fn setup<R: Rng>(r: &mut R) -> Result<Self::Parameters, Error>;
    fn evaluate(parameters: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error>;
}

pub trait FieldBasedHashParameters{
    type Fr: Field;
}

pub trait FieldBasedHash {
    type Data: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Data>;

    // Initialize the hash to a null state, or with `personalization` if specified.
    fn init(personalization: Option<&[Self::Data]>) -> Self;

    // Update the hash with `input`.
    fn update(&mut self, input: Self::Data) -> &mut Self;

    // Return the hash. This method is idempotent, and calling it multiple times will
    // give the same result. It's also possible to `update` with more inputs in between.
    fn finalize(&self) -> Self::Data;

    // Reset self to its initial state, allowing to change `personalization` too if needed.
    fn reset(&mut self, personalization: Option<&[Self::Data]>) -> &mut Self;
}

pub trait BatchFieldBasedHash {
    type Data: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Data>;

    // Given an `input_array` of size n * hash_rate, batches the computation of the n hashes
    // and outputs the n hash results.
    // NOTE: The hashes are independent from each other, therefore the output is not some sort
    // of aggregated hash but it's actually the hash result of each of the inputs, grouped in
    // hash_rate chunks.
    fn batch_evaluate(input_array: &[Self::Data]) -> Result<Vec<Self::Data>, Error>;

    // Given an `input_array` of size n * hash_rate, batches the computation of the n hashes
    // and outputs the n hash results.
    // Avoids a copy by requiring to pass the `output_array` already as input to the
    // function.
    // NOTE: The hashes are independent from each other, therefore the output is not some sort
    // of aggregated hash but it's actually the hash result of each of the inputs, grouped in
    // hash_rate chunks.
    fn batch_evaluate_in_place(input_array: &mut[Self::Data], output_array: &mut[Self::Data]);
}