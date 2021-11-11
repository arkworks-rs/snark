use algebra::{bytes::ToBytes, Field};
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::hash::Hash;

pub mod bowe_hopwood;
pub mod injective_map;
pub mod pedersen;

pub mod sbox;
pub use self::sbox::*;

pub mod poseidon;
pub use self::poseidon::*;

use crate::{CryptoError, Error};
use rayon::prelude::*;

pub trait FixedLengthCRH {
    const INPUT_SIZE_BITS: usize;
    type Output: ToBytes + Serialize + for<'a> Deserialize<'a> + Clone + Eq + Hash + Default;
    type Parameters: Clone + Serialize + for<'a> Deserialize<'a> + Default;

    fn setup<R: Rng>(r: &mut R) -> Result<Self::Parameters, Error>;
    fn evaluate(parameters: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error>;
}

/// Define parameters required to implement a hash function working with field arithmetics.
/// TODO: Depending on the hash construction some parameters may be present and others not
///       we should think about particularizing or generalizing this trait definition.
pub trait FieldBasedHashParameters: Clone {
    type Fr: Field;

    /// The rate of the hash function
    const R: usize;
}

pub trait FieldBasedHash {
    type Data: Field;
    type Parameters: FieldBasedHashParameters<Fr = Self::Data>;

    /// Initialize a Field Hash accepting inputs of constant length `input_size`:
    /// any attempt to finalize the hash after having updated the Self instance
    /// with a number of inputs not equal to `input_size` should result in an error.
    /// Initialize the hash to a null state, or with `personalization` if specified.
    fn init_constant_length(input_size: usize, personalization: Option<&[Self::Data]>) -> Self;

    /// Initialize a Field Hash accepting inputs of variable length.
    /// It is able to serve two different modes, selected by the boolean `mod_rate`:
    /// - `mod_rate` = False is for the ususal variable length hash, whereas
    /// - `mod_rate` = True allows the input only to be a multiple of the rate (and hence
    /// should throw an error when trying to finalize with a non-multiple of rate input).
    /// This mode allows an optimized handling of padding, saving constraints in SNARK applications;
    fn init_variable_length(mod_rate: bool, personalization: Option<&[Self::Data]>) -> Self;

    /// Update the hash with `input`.
    fn update(&mut self, input: Self::Data) -> &mut Self;

    /// Return the hash. This method is idempotent, and calling it multiple times will
    /// give the same result.
    fn finalize(&self) -> Result<Self::Data, Error>;

    /// Reset self to its initial state, allowing to change `personalization` too if needed.
    fn reset(&mut self, personalization: Option<&[Self::Data]>) -> &mut Self;
}

/// Helper allowing to hash the implementor of this trait into a Field
pub trait FieldHasher<F: Field, H: FieldBasedHash<Data = F>> {
    /// Hash `self`, given some optional `personalization` into a Field
    fn hash(&self, personalization: Option<&[H::Data]>) -> Result<H::Data, Error>;
}

pub trait BatchFieldBasedHash {
    type Data: Field;

    /// Specification of this type allows to provide a default implementation and more flexibility
    /// when included in other traits/struct (i.e. a FieldBasedMerkleTree using both simple and
    /// batch hashes can only specify this trait, simplifying its design and usage).
    /// Still, it's a reasonable addition for a trait like this.
    type BaseHash: FieldBasedHash<Data = Self::Data>;

    /// Given an `input_array` of size n * hash_rate, batches the computation of the n hashes
    /// and outputs the n hash results.
    /// NOTE: The hashes are independent from each other, therefore the output is not some sort
    /// of aggregated hash but it's actually the hash result of each of the inputs, grouped in
    /// hash_rate chunks.
    fn batch_evaluate(input_array: &[Self::Data]) -> Result<Vec<Self::Data>, Error> {
        let rate = <<Self::BaseHash as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;
        if input_array.len() % rate != 0 {
            return Err(Box::new(CryptoError::Other(
                "The length of the input data array is not a multiple of the rate".to_owned(),
            )));
        }
        if input_array.len() == 0 {
            return Err(Box::new(CryptoError::Other(
                "Input data array does not contain any data".to_owned(),
            )));
        }

        Ok(input_array
            .par_chunks(rate)
            .map(|chunk| {
                let mut digest =
                    <Self::BaseHash as FieldBasedHash>::init_constant_length(rate, None);
                chunk.iter().for_each(|input| {
                    digest.update(input.clone());
                });
                digest.finalize().unwrap()
            })
            .collect::<Vec<_>>())
    }

    /// Given an `input_array` of size n * hash_rate, batches the computation of the n hashes
    /// and outputs the n hash results.
    /// Avoids a copy by requiring to pass the `output_array` already as input to the
    /// function.
    /// NOTE: The hashes are independent from each other, therefore the output is not some sort
    /// of aggregated hash but it's actually the hash result of each of the inputs, grouped in
    /// hash_rate chunks.
    fn batch_evaluate_in_place(
        input_array: &mut [Self::Data],
        output_array: &mut [Self::Data],
    ) -> Result<(), Error> {
        let rate = <<Self::BaseHash as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;
        if input_array.len() % rate != 0 {
            return Err(Box::new(CryptoError::Other(
                "The length of the input data array is not a multiple of the rate".to_owned(),
            )));
        }
        if input_array.len() == 0 {
            return Err(Box::new(CryptoError::Other(
                "Input data array does not contain any data".to_owned(),
            )));
        }
        if output_array.len() != input_array.len() / rate {
            return Err(Box::new(CryptoError::Other(format!(
                "Output array size must be equal to input_array_size/rate. Output array size: {}, Input array size: {}, Rate: {}",
                output_array.len(),
                input_array.len(),
                rate
            ))));
        }
        input_array
            .par_chunks(rate)
            .zip(output_array.par_iter_mut())
            .for_each(|(inputs, output)| {
                let mut digest =
                    <Self::BaseHash as FieldBasedHash>::init_constant_length(rate, None);
                inputs.iter().for_each(|input| {
                    digest.update(input.clone());
                });
                *output = digest.finalize().unwrap();
            });
        Ok(())
    }
}

#[cfg(test)]
mod test {

    use algebra::{fields::mnt4753::Fr as MNT4753Fr, Field, UniformRand};

    use super::BatchFieldBasedHash;
    use crate::crh::poseidon::{MNT4BatchPoseidonHash, MNT4PoseidonHash};

    use crate::{FieldBasedHash, FieldBasedHashParameters};
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    struct DummyMNT4BatchPoseidonHash;

    impl BatchFieldBasedHash for DummyMNT4BatchPoseidonHash {
        type Data = MNT4753Fr;
        type BaseHash = MNT4PoseidonHash;
    }

    pub(crate) fn constant_length_field_based_hash_test<H: FieldBasedHash>(
        digest: &mut H,
        inputs: Vec<H::Data>,
    ) {
        let inputs_len = inputs.len();

        let final_elem = inputs[inputs_len - 1].clone();

        digest.reset(None);
        inputs.into_iter().take(inputs_len - 1).for_each(|fe| {
            digest.update(fe);
        });

        // Test call to finalize() with too few inputs with respect to the declared size
        // results in an error.
        assert!(
            digest.finalize().is_err(),
            "Success call to finalize despite smaller number of inputs"
        );

        //Test finalize() being idempotent
        digest.update(final_elem);
        let output = digest.finalize().unwrap();
        assert_eq!(
            output,
            digest.finalize().unwrap(),
            "Two subsequent calls to finalize gave different results"
        );

        // Test call to finalize() with too much inputs with respect to the declared size
        // results in an error.
        digest.update(final_elem);
        assert!(
            digest.finalize().is_err(),
            "Success call to finalize despite bigger number of inputs"
        );
    }

    pub(crate) fn variable_length_field_based_hash_test<H: FieldBasedHash>(
        digest: &mut H,
        inputs: Vec<H::Data>,
        mod_rate: bool,
    ) {
        let rate = <H::Parameters as FieldBasedHashParameters>::R;

        let pad_inputs = |mut inputs: Vec<H::Data>| -> Vec<H::Data> {
            inputs.push(H::Data::one());
            inputs.append(&mut vec![H::Data::zero(); rate - (inputs.len() % rate)]);
            inputs
        };

        if mod_rate {
            constant_length_field_based_hash_test(digest, inputs);
        } else {
            // Check padding is added correctly and that the hash is collision free when input
            // is not modulus rate
            let output = digest.finalize().unwrap();
            let padded_inputs = pad_inputs(inputs.clone());
            digest.reset(None);
            padded_inputs.iter().for_each(|fe| {
                digest.update(fe.clone());
            });
            assert_ne!(
                output,
                digest.finalize().unwrap(),
                "Incorrect padding: collision detected"
            );

            // Check padding is added correctly and that the hash is collision free also when input
            // happens to be modulus rate
            let output = digest.finalize().unwrap();
            let padded_inputs = pad_inputs(padded_inputs);
            digest.reset(None);
            padded_inputs.into_iter().for_each(|fe| {
                digest.update(fe);
            });
            assert_ne!(
                output,
                digest.finalize().unwrap(),
                "Incorrect padding: collision detected"
            );
        }
    }

    #[ignore]
    #[test]
    fn test_default_batch_hash_implementation() {
        let rate = 2;
        let num_inputs = 100;
        let mut inputs = Vec::with_capacity(num_inputs);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..num_inputs {
            inputs.push(MNT4753Fr::rand(&mut rng))
        }

        let batch_hash_output = MNT4BatchPoseidonHash::batch_evaluate(inputs.as_slice()).unwrap();
        let dummy_batch_hash_output =
            DummyMNT4BatchPoseidonHash::batch_evaluate(inputs.as_slice()).unwrap();
        assert_eq!(batch_hash_output, dummy_batch_hash_output);

        let mut batch_hash_output_new = vec![MNT4753Fr::zero(); num_inputs / rate];
        let mut dummy_batch_hash_output_new = vec![MNT4753Fr::zero(); num_inputs / rate];

        MNT4BatchPoseidonHash::batch_evaluate_in_place(
            inputs.as_mut_slice(),
            batch_hash_output_new.as_mut_slice(),
        )
        .unwrap();
        DummyMNT4BatchPoseidonHash::batch_evaluate_in_place(
            inputs.as_mut_slice(),
            dummy_batch_hash_output_new.as_mut_slice(),
        )
        .unwrap();

        assert_eq!(batch_hash_output_new, dummy_batch_hash_output_new);
        assert_eq!(batch_hash_output, batch_hash_output_new);
        assert_eq!(dummy_batch_hash_output, dummy_batch_hash_output_new);
    }
}
