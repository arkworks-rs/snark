use algebra::{
    Field, PrimeField, bytes::ToBytes
};
use rand::Rng;
use std::hash::Hash;
use serde::{Serialize, Deserialize};

pub mod bowe_hopwood;
pub mod injective_map;
pub mod pedersen;

pub mod sbox;
pub use self::sbox::*;

pub mod poseidon;
pub use self::poseidon::*;

use crate::Error;
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

    /// Initialize the hash to a null state, or with `personalization` if specified.
    fn init(personalization: Option<&[Self::Data]>) -> Self;

    /// Update the hash with `input`.
    fn update(&mut self, input: Self::Data) -> &mut Self;

    /// Return the hash. This method is idempotent, and calling it multiple times will
    /// give the same result. It's also possible to `update` with more inputs in between.
    fn finalize(&self) -> Self::Data;

    /// Reset self to its initial state, allowing to change `personalization` too if needed.
    fn reset(&mut self, personalization: Option<&[Self::Data]>) -> &mut Self;
}

/// Helper allowing to hash the implementor of this trait into a Field
pub trait FieldHasher<F: Field, H: FieldBasedHash<Data = F>> {

    /// Hash `self`, given some optional `personalization` into a Field
    fn hash(
        &self,
        personalization: Option<&[H::Data]>
    ) -> Result<H::Data, Error>;
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
        assert_eq!(input_array.len() % rate, 0, "The length of the input data array is not a multiple of the rate.");
        assert_ne!(input_array.len(), 0, "Input data array does not contain any data.");

        Ok(input_array.par_chunks(rate).map(|chunk| {
            let mut digest = <Self::BaseHash as FieldBasedHash>::init(None);
            chunk.iter().for_each(|input| { digest.update(input.clone()); } );
            digest.finalize()
        }).collect::<Vec<_>>())
    }

    /// Given an `input_array` of size n * hash_rate, batches the computation of the n hashes
    /// and outputs the n hash results.
    /// Avoids a copy by requiring to pass the `output_array` already as input to the
    /// function.
    /// NOTE: The hashes are independent from each other, therefore the output is not some sort
    /// of aggregated hash but it's actually the hash result of each of the inputs, grouped in
    /// hash_rate chunks.
    fn batch_evaluate_in_place(input_array: &mut[Self::Data], output_array: &mut[Self::Data]) {
        let rate = <<Self::BaseHash as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;
        input_array.par_chunks(rate).zip(output_array.par_iter_mut())
            .for_each(|(inputs, output)| {
                let mut digest = <Self::BaseHash as FieldBasedHash>::init(None);
                inputs.iter().for_each(|input| { digest.update(input.clone()); } );
                *output = digest.finalize();
            });
    }
}

#[derive(Clone, Debug)]
pub enum SpongeMode {
    Absorbing,
    Squeezing,
}

/// the trait for algebraic sponge
pub trait AlgebraicSponge<F: PrimeField>: Clone {
    /// Initialize the sponge
    fn new() -> Self;
    /// Initialize the sponge given a state
    fn from_state(state: Vec<F>) -> Self;
    /// Get the sponge internal state
    fn get_state(&self) -> &[F];
    /// Update the sponge with `elems`
    fn absorb(&mut self, elems: Vec<F>);
    /// Output `num` field elements from the sponge.
    fn squeeze(&mut self, num: usize) -> Vec<F>;
    /// Reset the sponge to its initial state
    fn reset(&mut self) {
        *self = Self::new();
    }
}

#[cfg(test)]
mod test {

    use algebra::{
        fields::mnt4753::Fr as MNT4753Fr, Field, PrimeField, UniformRand
    };

    use super::BatchFieldBasedHash;
    use crate::crh::poseidon::{
        MNT4PoseidonHash, MNT4BatchPoseidonHash
    };

    use rand_xorshift::XorShiftRng;
    use rand::{SeedableRng, thread_rng};
    use crate::{FieldBasedHash, AlgebraicSponge};
    use std::collections::HashSet;

    struct DummyMNT4BatchPoseidonHash;

    impl BatchFieldBasedHash for DummyMNT4BatchPoseidonHash {
        type Data = MNT4753Fr;
        type BaseHash = MNT4PoseidonHash;
    }

    pub(crate) fn field_based_hash_test<H: FieldBasedHash>(
        personalization: Option<&[H::Data]>,
        inputs: Vec<H::Data>,
        expected_output: H::Data
    )
    {
        // Test H(inputs) == expected_output
        let mut digest = H::init(personalization);
        inputs.iter().for_each(|fe| { digest.update(fe.clone()); });
        let output = digest.finalize();
        assert_eq!(output, expected_output, "Outputs do not match");

        let inputs_len = inputs.len();
        if inputs_len > 1 {
            let final_elem = inputs[inputs_len - 1].clone();
            // Test finalize() holding the state and allowing updates in between different calls to it
            digest.reset(None);
            inputs.into_iter().take(inputs_len - 1).for_each(|fe| { digest.update(fe); });

            digest.finalize();
            digest.update(final_elem);
            assert_eq!(output, digest.finalize());

            //Test finalize() being idempotent
            assert_eq!(output, digest.finalize());
        }
    }

    pub(crate) fn algebraic_sponge_test<H: AlgebraicSponge<F>, F: PrimeField>(
        to_absorb: Vec<F>,
        expected_squeeze: F
    )
    {
        let mut sponge = H::new();

        // Absorb all field elements
        sponge.absorb(to_absorb);

        // Squeeze and check the output
        assert_eq!(expected_squeeze, sponge.squeeze(1)[0]);

        // Check that calling squeeze() multiple times without absorbing
        // changes the output
        let mut prev = expected_squeeze;
        for _ in 0..100 {
            let curr = sponge.squeeze(1)[0];
            assert!(prev != curr);
            prev = curr;
        }

        let rng = &mut thread_rng();

        // Check squeeze() outputs the correct number of field elements
        // all different from each others
        let mut set = HashSet::new();
        for i in 0..=10 {
            sponge.absorb(vec![F::rand(rng); i]);
            let outs = sponge.squeeze(i);
            assert_eq!(i, outs.len());

            // HashSet::insert(val) returns false if val was already present, so to check
            // that all the elements output by the sponge are different, we assert insert()
            // returning always true
            outs.into_iter().for_each(|f| assert!(set.insert(f)));
        }

        //Test edge cases. Assumption: R = 2
        sponge.reset();

        // Absorb nothing. Check that the internal state is not changed.
        let prev_state = sponge.get_state().to_vec();
        sponge.absorb(vec![]);
        assert_eq!(prev_state, sponge.get_state());

        // Squeeze nothing. Check that the internal state is not changed.
        let prev_state = sponge.get_state().to_vec();
        sponge.squeeze(0);
        assert_eq!(prev_state, sponge.get_state());

        // Absorb up to rate elements and trigger a permutation. Assert that calling squeeze()
        // afterwards won't trigger another permutation.
        sponge.absorb(vec![F::rand(rng); 2]);
        let prev_state = sponge.get_state().to_vec();
        sponge.squeeze(1);
        let curr_state = sponge.get_state().to_vec();
        assert_eq!(prev_state, curr_state);

        // The next squeeze() should instead change the state
        sponge.squeeze(1);
        assert!(curr_state != sponge.get_state());
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
        let dummy_batch_hash_output = DummyMNT4BatchPoseidonHash::batch_evaluate(inputs.as_slice()).unwrap();
        assert_eq!(batch_hash_output, dummy_batch_hash_output);

        let mut batch_hash_output_new = vec![MNT4753Fr::zero(); num_inputs/rate];
        let mut dummy_batch_hash_output_new = vec![MNT4753Fr::zero(); num_inputs/rate];

        MNT4BatchPoseidonHash::batch_evaluate_in_place(inputs.as_mut_slice(), batch_hash_output_new.as_mut_slice());
        DummyMNT4BatchPoseidonHash::batch_evaluate_in_place(inputs.as_mut_slice(), dummy_batch_hash_output_new.as_mut_slice());

        assert_eq!(batch_hash_output_new, dummy_batch_hash_output_new);
        assert_eq!(batch_hash_output, batch_hash_output_new);
        assert_eq!(dummy_batch_hash_output, dummy_batch_hash_output_new);
    }
}