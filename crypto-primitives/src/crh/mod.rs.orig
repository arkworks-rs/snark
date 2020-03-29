use algebra::bytes::ToBytes;
use rand::Rng;
use std::hash::Hash;

pub mod bowe_hopwood;
pub mod injective_map;
pub mod pedersen;
pub mod poseidon;

use crate::Error;

#[cfg(feature = "r1cs")]
pub mod constraints;
#[cfg(feature = "r1cs")]
pub use constraints::*;

pub trait FixedLengthCRH {
    const INPUT_SIZE_BITS: usize;
    type Output: ToBytes + Clone + Eq + Hash + Default;
    type Parameters: Clone + Default;

    fn setup<R: Rng>(r: &mut R) -> Result<Self::Parameters, Error>;
    fn evaluate(parameters: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error>;
}

// pub trait Batched2to1CRH {
//     const INPUT_NUM_PAIRS: usize;
//     type Output: ToBytes + Clone + Eq + Hash + Default;
//     type Parameters: Clone + 'static;
//
//     fn evaluate(input: &[u8]) -> Result<Self::Output, Error>;
// }

use algebra::Field;

pub trait FieldBasedHashParameters: Sized{
    type Fr: Field;
}

pub trait FieldBasedHash {
    type Data: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Data>;

    fn evaluate(input: &[Self::Data]) -> Result<Self::Data, Error>;
    //fn batch_evaluate(input: Vec<Vec<Self::Data>>) -> Result<Vec<Self::Data>, Error>;
}

pub trait BatchFieldBasedHash {
    type Data: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Data>;

    fn batch_evaluate_2_1(input_array: &mut[Self::Data]);
}

