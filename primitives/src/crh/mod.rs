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

    fn batch_evaluate(input_array: &mut[Self::Data], output_array: &mut[Self::Data]);

    fn merkle_tree(input_vec: &mut[Self::Data], output_vec: &mut[Self::Data], input_size: usize);
    fn merkle_tree_2_1(array: &mut [Self::Data], input_size: usize);

}


