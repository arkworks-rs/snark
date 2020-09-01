use std::{clone::Clone, fmt::Debug};

pub mod naive;
pub use self::naive::*;

pub mod poseidon;
pub use self::poseidon::*;

#[cfg(feature = "smt")]
pub mod smt;
#[cfg(feature = "smt")]
pub use self::smt::*;

use algebra::{
    biginteger::BigInteger768, fields::mnt4753::Fr as MNT4753Fr, field_new, Field};

pub trait FieldBasedMerkleTreeParameters: 'static + Clone {
    type Data: Field;

    // The arity of the Sparse Merkle Tree
    const MERKLE_ARITY: usize;
    // The pre-computed hashes of the empty nodes for the different levels of the SMT
    const EMPTY_HASH_CST: &'static [Self::Data];
}

pub trait FieldBasedMerkleTree: Clone {
    type Data: Field;
    type MerklePath: Clone + Debug;
    type Parameters: FieldBasedMerkleTreeParameters;

    // Initialize this tree. The user must pass the maximum number of leaves the tree must
    // support.
    fn init(num_leaves: usize) -> Self;

    // This function reset the pointing indices for the new element positions and processed positions
    // Mainly used for benchmarking
    fn reset(&mut self) -> &mut Self;

    // This function appends a new leaf to the Merkle tree
    fn append(&mut self, leaf: Self::Data) -> &mut Self;

    // This function finalizes the computation of the Merkle tree and returns an updated
    // copy of it. This method is idempotent, and calling it multiple times will
    // give the same result. It's also possible to `update` with more inputs in between.
    fn finalize(&self) -> Self;

    // This function finalizes the computation of the Merkle tree
    // Once this function is called, it is not possible to further update the tree.
    fn finalize_in_place(&mut self) -> &mut Self;

    // Returns the root of the Merkle Tree. Returns None if the tree has not been
    // finalized before calling this function.
    fn root(&self) -> Option<Self::Data>;

    // Given an `index` returns the MerklePath of the leaf at that index up until the root of the
    // `self` Merkle Tree. Returns None if the tree has not been finalized before
    // calling this function.
    fn get_merkle_path(&self, leaf_index: usize) -> Option<Self::MerklePath>;

    // Verify `path` given `leaf`. Returns None if the tree has not been
    // finalized before calling this function.
    fn verify_merkle_path(
        &self,
        leaf: &Self::Data,
        path: &Self::MerklePath
    ) -> Option<bool>;
}

// PoseidonHash("This represents an empty Merkle Root for a MNT4753PoseidonHash based Merkle Tree.") padded with 0s
pub const MNT4753_PHANTOM_MERKLE_ROOT: MNT4753Fr =
    field_new!(MNT4753Fr, BigInteger768([
        8534937690304963668,
        5486292534803213323,
        1720870611961422927,
        11405914840660719672,
        7162329517212056783,
        11658292353137306079,
        17490588101047840223,
        12735752881395833110,
        11735157047413601083,
        6658060155531600932,
        1470933043432945054,
        312822709740712
    ]));


    #[ignore]
    #[test]
    fn generate_mnt4753_phantom_merkle_root(){
        use crate::crh::{FieldBasedHash, poseidon::MNT4PoseidonHash};
        use algebra::{FromBytes, PrimeField, FpParameters};

        let field_size_in_bytes = (MNT4753Fr::size_in_bits() + (<MNT4753Fr as PrimeField>::Params::REPR_SHAVE_BITS as usize))/8;
        let magic_string = b"This represents an empty Merkle Root for a MNT4753PoseidonHash based Merkle Tree.";

        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(magic_string);
        for _ in magic_string.len()..field_size_in_bytes { hash_input.push(0u8) }
        let hash_input_f = MNT4753Fr::read(hash_input.as_slice()).unwrap();

        let hash = MNT4PoseidonHash::init(None).update(hash_input_f).finalize();
        assert_eq!(hash, MNT4753_PHANTOM_MERKLE_ROOT);
    }