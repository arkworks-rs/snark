pub mod path;
pub use self::path::*;

pub mod naive;
pub use self::naive::*;

pub mod optimized;
pub use self::optimized::*;

#[cfg(feature = "smt")]
pub mod smt;
#[cfg(feature = "smt")]
pub use self::smt::*;

pub mod parameters;
pub use self::parameters::*;

use algebra::{
    Field, ToBytes, FromBytes,
};
use crate::{FieldBasedHash, BatchFieldBasedHash, Error, FieldBasedHashParameters, MerkleTreeError};
use std::{clone::Clone, fmt::Debug};
use serde::{Serialize, Deserialize};

/// Definition of parameters needed to implement and optimize a Merkle Tree whose nodes and leaves
/// are Field elements. The trait is generic with respect to the arity of the Merkle Tree.
pub trait FieldBasedMerkleTreeParameters: 'static + Clone {
    type Data: Field; /// Actually unnecessary, but simplifies the overall design
    type H: FieldBasedHash<Data = Self::Data>;
    /// The arity of the Merkle Tree
    const MERKLE_ARITY: usize;
    /// The pre-computed hashes of the empty nodes for the different levels of the Merkle Tree
    const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>>;
}

/// Pre-computed hashes of the empty nodes for the different levels of the Merkle Tree
pub struct FieldBasedMerkleTreePrecomputedEmptyConstants<'a, H: FieldBasedHash> {
    pub nodes: &'a [H::Data],
    pub merkle_arity: usize,
    pub max_height: usize,
}

/// For optimized Merkle Tree implementations, it provides the possibility to specify
/// a hash function able to efficiently compute many hashes at the same time
pub trait BatchFieldBasedMerkleTreeParameters: FieldBasedMerkleTreeParameters {
    type BH: BatchFieldBasedHash<
        Data = <Self as FieldBasedMerkleTreeParameters>::Data,
        BaseHash = <Self as FieldBasedMerkleTreeParameters>::H
    >;
}

pub(crate) fn check_precomputed_parameters<T: FieldBasedMerkleTreeParameters>(tree_height: usize) -> bool
{
    match T::EMPTY_HASH_CST {
        Some(supported_params) => {
            tree_height <= supported_params.max_height &&
                T::MERKLE_ARITY == supported_params.merkle_arity &&
                T::MERKLE_ARITY == <<T::H as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R
        }
        None => false
    }
}

/// Definition of a Merkle Tree whose leaves and nodes are Field Elements. The trait is
/// designed to be efficient in memory by providing a digest-like interface that allows
/// to update the tree one leaf at a time, without the necessity to keep all the leaves
/// in memory in order to create the tree. The trait is generic with respect to the arity
/// of the Merkle Tree and to the hash function used.
pub trait FieldBasedMerkleTree: Clone {
    type Parameters: FieldBasedMerkleTreeParameters;
    type MerklePath: FieldBasedMerkleTreePath<
        H = <Self::Parameters as FieldBasedMerkleTreeParameters>::H,
        Parameters = Self::Parameters
    >;

    /// Append a new leaf to the Merkle Tree. The moment in which the root will be computed
    /// is transparent to the user and obeys to pre-defined internal policies.
    fn append(&mut self, leaf: <Self::Parameters as FieldBasedMerkleTreeParameters>::Data) -> &mut Self;

    /// Force the computation of the root whatever its internal state and return an updated copy
    /// of the Merkle Tree. This function is idempotent, i.e. calling it multiple times will give
    /// the same result. It's also possible to `update` with more inputs in between.
    fn finalize(&self) -> Self;

    /// Force the computation of the root whatever its internal state and return an updated copy
    /// of the Merkle Tree. It's more efficient than `finalize` because avoids a copy; however,
    /// once this function is called, it is not possible to further `update` the tree.
    fn finalize_in_place(&mut self) -> &mut Self;

    /// Resets the internal state of the tree, bringing it back to the initial one.
    fn reset(&mut self) -> &mut Self;

    /// Return the root of the Merkle Tree. Return None if the tree has not been
    /// finalized before calling this function.
    fn root(&self) -> Option<<Self::Parameters as FieldBasedMerkleTreeParameters>::Data>;

    /// Given an `index` returns the MerklePath of the leaf at that index up until the root of the
    /// Merkle Tree. Returns None if the tree has not been finalized before calling this function.
    fn get_merkle_path(&self, leaf_index: usize) -> Option<Self::MerklePath>;

    /// Returns the height of this Merkle Tree.
    fn height(&self) -> usize;
}

/// Definition of a Merkle Path for a Merkle Tree whose leaves and nodes are field elements. The
/// trait is generic with respect to the arity of the Merkle Tree and to the hash function used.
pub trait FieldBasedMerkleTreePath {
    type H: FieldBasedHash;
    type Path: Clone + Debug + Serialize + for<'a> Deserialize<'a>;
    type Parameters: FieldBasedMerkleTreeParameters<
        Data = <Self::H as FieldBasedHash>::Data,
        H = Self::H
    >;

    /// Return a new instance of the struct implementing this trait given the raw `path`
    fn new(path: Self::Path) -> Self;

    /// Verify the Merkle Path for `leaf` given the `root` of a Merkle Tree with height `height`.
    fn verify(
        &self,
        height: usize,
        leaf: &<Self::H as FieldBasedHash>::Data,
        root: &<Self::H as FieldBasedHash>::Data
    ) -> Result<bool, Error> {
        let path_len = self.get_length();
        if path_len != height {
            Err(MerkleTreeError::IncorrectPathLength(path_len, height))?
        }
        self.verify_without_length_check(leaf, root)
    }

    /// Verify the Merkle Path for `leaf` given the `root` of a Merkle Tree. Doesn't check if the
    /// length of the Merkle Path is consistent with the height of the corresponding Merkle Tree,
    /// therefore it is advisable to use it when it's certain that `leaf` is actually a leaf.
    fn verify_without_length_check(
        &self,
        leaf: &<Self::H as FieldBasedHash>::Data,
        root: &<Self::H as FieldBasedHash>::Data
    ) -> Result<bool, Error>;

    /// Returns the underlying raw path
    fn get_raw_path(&self) -> Self::Path;

    /// Returns the length of the underlying raw path
    fn get_length(&self) -> usize;
}