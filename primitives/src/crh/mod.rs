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

    fn evaluate(input: &[Self::Data]) -> Result<Self::Data, Error>;
}

pub trait BatchFieldBasedHash {
    type Data: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Data>;

    fn batch_evaluate(input_array: &[Self::Data]) -> Result<Vec<Self::Data>, Error>;
    fn batch_evaluate_in_place(input_array: &mut[Self::Data], output_array: &mut[Self::Data]);

}

pub trait UpdatableFieldBasedHash {
    type Data: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Data>;

    // Updates the hash with `input`
    fn update(&mut self, input: &Self::Data) -> &mut Self;

    // Returns the hash. This method is idempotent, and calling it multiple times will
    // give the same result. It's also possible to `update` with more inputs in between.
    fn finalize(&self) -> Self::Data;
}

pub trait UpdatableBatchFieldBasedHash {
    type Data: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Data>;

    // Updates the batches hash outputs with `input` batch
    fn update(&mut self, input: &[Self::Data]) -> &mut Self;

    // Returns the hash output for each batch. This method is idempotent, and calling it multiple
    // times will give the same result. It's also possible to `update` with more input batches
    // in between.
    fn finalize(&self) -> Vec<Self::Data>;

    // Will free memory from the batches hash outputs computed until now.
    fn clear(&mut self);
}