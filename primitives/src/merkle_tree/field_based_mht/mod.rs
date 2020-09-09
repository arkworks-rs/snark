use std::{clone::Clone, fmt::Debug};

pub mod naive;
pub use self::naive::*;

pub mod optimized;
pub use self::optimized::*;

pub mod poseidon;
pub use self::poseidon::*;

#[cfg(feature = "smt")]
pub mod smt;
#[cfg(feature = "smt")]
pub use self::smt::*;

use algebra::Field;
use crate::{FieldBasedHash, BatchFieldBasedHash, Error, FieldBasedHashParameters, MerkleTreeError};

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
    type Path: Clone + Debug;
    type Parameters: FieldBasedMerkleTreeParameters<
        Data = <Self::H as FieldBasedHash>::Data,
        H = Self::H
    >;

    fn new(path: Self::Path) -> Self;

    /// Verify the Merkle Path for `leaf` given the `root` of the Merkle Tree and its `height`.
    fn verify(
        &self,
        height: usize,
        leaf: &<Self::H as FieldBasedHash>::Data,
        root: &<Self::H as FieldBasedHash>::Data
    ) -> Result<bool, Error>;

    /// Returns the underlying raw path
    fn get_raw_path(&self) -> Self::Path;
}

/// An implementation of the FieldBasedMerkleTreePath trait, for a given FieldBasedHash and
/// FieldBasedMerkleTree with arbitrary arity.
/// TODO: Test for arity > 2
#[derive(Clone, Debug)]
pub struct FieldBasedMHTPath<T: FieldBasedMerkleTreeParameters>{
    path: Vec<(Vec<<T::H as FieldBasedHash>::Data>, usize)>,
}

impl<T: FieldBasedMerkleTreeParameters> PartialEq for FieldBasedMHTPath<T> {
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path
    }
}

impl<T: FieldBasedMerkleTreeParameters> FieldBasedMerkleTreePath for FieldBasedMHTPath<T> {
    type H = T::H;

    /// A Merkle Path for a leaf of a Merkle Tree with arity >= 2 will be made up of couples of nodes
    /// and an integer. The nodes are all the siblings of the leaf ( in number MERKLE_ARITY - 1 )
    /// and the integer is the position of the leaf, among its siblings, in the input of the hash
    /// function (valid values are from 0 to MERKLE_ARITY - 1).
    type Path = Vec<(Vec<<Self::H as FieldBasedHash>::Data>, usize)>;

    type Parameters = T;

    fn new(path: Self::Path) -> Self {
        Self { path }
    }

    fn verify(
        &self,
        height: usize,
        leaf: &<Self::H as FieldBasedHash>::Data,
        root: &<Self::H as FieldBasedHash>::Data
    ) -> Result<bool, Error> {
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design. Should be also enforced by the
        // MerkleTree that creates this instance, but let's do it again.
        let height = height + 1;
        assert_eq!(<<Self::H as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R, T::MERKLE_ARITY);
        if self.path.len() != height - 1 {
            Err(MerkleTreeError::IncorrectPathLength(self.path.len(), height - 1))?
        }
        let mut digest = <Self::H as FieldBasedHash>::init(None);
        let mut prev_node = leaf.clone();
        for (sibling_nodes, position) in self.path.as_slice() {
            assert_eq!(sibling_nodes.len(), T::MERKLE_ARITY - 1);

            // Get the position of the node among its siblings
            let prev_node_position = *position % T::MERKLE_ARITY;

            // Update the digest respecting the position of each sibling
            for i in 0..T::MERKLE_ARITY {
                if i == prev_node_position {
                    digest.update(prev_node.clone());
                } else {
                    let index = i % (T::MERKLE_ARITY - 1); // Be sure to not overflow the siblings vector
                    digest.update(sibling_nodes[index]);
                }
            }

            // Compute the parent node
            prev_node = digest.finalize();
            digest.reset(None);
        }

        // Check final computed node is equal to the root
        if prev_node == root.clone() {
            Ok(true)
        } else {
            Ok(false)
        }
    }

    fn get_raw_path(&self) -> Self::Path {
        self.path.clone()
    }
}

/// A wrapper around a Merkle Path for a FieldBasedMerkleTree of arity 2. Merkle Trees of arity
/// 2 are the most common and it's worth to explicitly create a separate struct
#[derive(Clone, Debug)]
pub struct FieldBasedBinaryMHTPath<T: FieldBasedMerkleTreeParameters>{
    path: Vec<(<T::H as FieldBasedHash>::Data, bool)>,
}

impl<T: FieldBasedMerkleTreeParameters> FieldBasedMerkleTreePath for FieldBasedBinaryMHTPath<T> {
    type H = T::H;
    type Path = Vec<(<T::H as FieldBasedHash>::Data, bool)>;
    type Parameters = T;

    fn new(path: Self::Path) -> Self {
        Self { path }
    }

    fn verify(
        &self,
        height: usize,
        leaf: &<Self::H as FieldBasedHash>::Data,
        root: &<Self::H as FieldBasedHash>::Data
    ) -> Result<bool, Error> {
        let mut v = Vec::with_capacity(self.path.len());
        for &(node, direction) in &self.path {
            v.push((vec![node], if !direction {0} else {1}));
        }
        FieldBasedMHTPath::<T>::new(v).verify(height, leaf, root)
    }

    fn get_raw_path(&self) -> Self::Path {
        self.path.clone()
    }
}

impl<T: FieldBasedMerkleTreeParameters> PartialEq for FieldBasedBinaryMHTPath<T> {
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path
    }
}

impl<T: FieldBasedMerkleTreeParameters> PartialEq<FieldBasedBinaryMHTPath<T>> for FieldBasedMHTPath<T> {
    fn eq(&self, other: &FieldBasedBinaryMHTPath<T>) -> bool {
        let binary_path = other.get_raw_path();
        if self.path.len() != binary_path.len() { return false };

        for ((p1_node, p1_pos), (p2_node, p2_pos)) in
            self.path.iter().zip(binary_path) {
            if  p1_node.len() != 1 || // In a binary Merkle Tree, there is only one sibling for each node
                p1_node[0] != p2_node || // The two nodes along the same leves of the paths must be equal
                (*p1_pos != 0 && *p1_pos != 1) || // In a binary Merkle Tree, position of the acutal node is either 0 or 1
                (*p1_pos == 0 && p2_pos) || (*p1_pos == 1 && !p2_pos) // Second path position is expressed as a boolean flag
            {
                return false
            }
        }
        return true;
    }
}