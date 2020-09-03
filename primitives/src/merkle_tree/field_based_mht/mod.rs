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
use std::marker::PhantomData;

pub trait FieldBasedMerkleTreeParameters: 'static + Clone {
    type Data: Field;
    // The arity of the Sparse Merkle Tree
    const MERKLE_ARITY: usize;
    // The pre-computed hashes of the empty nodes for the different levels of the SMT
    const EMPTY_HASH_CST: &'static [Self::Data];
}

pub trait BaseFieldBasedMerkleTreeParameters: FieldBasedMerkleTreeParameters {
    type H: FieldBasedHash<Data = <Self as FieldBasedMerkleTreeParameters>::Data>;
}

pub trait BatchFieldBasedMerkleTreeParameters: FieldBasedMerkleTreeParameters {
    type H: BatchFieldBasedHash<Data = <Self as FieldBasedMerkleTreeParameters>::Data>;
}

pub trait FieldBasedMerkleTree: Clone {
    type Parameters: FieldBasedMerkleTreeParameters;
    type MerklePath: FieldBasedMerkleTreePath<Data = <Self::Parameters as FieldBasedMerkleTreeParameters>::Data>;

    // Initialize this tree. The user must pass the maximum number of leaves the tree must
    // support.
    fn init(num_leaves: usize) -> Self;

    // This function reset the pointing indices for the new element positions and processed positions
    // Mainly used for benchmarking
    fn reset(&mut self) -> &mut Self;

    // This function appends a new leaf to the Merkle tree
    fn append(&mut self, leaf: <Self::Parameters as FieldBasedMerkleTreeParameters>::Data) -> &mut Self;

    // This function finalizes the computation of the Merkle tree and returns an updated
    // copy of it. This method is idempotent, and calling it multiple times will
    // give the same result. It's also possible to `update` with more inputs in between.
    fn finalize(&self) -> Self;

    // This function finalizes the computation of the Merkle tree
    // Once this function is called, it is not possible to further update the tree.
    fn finalize_in_place(&mut self) -> &mut Self;

    // Returns the root of the Merkle Tree. Returns None if the tree has not been
    // finalized before calling this function.
    fn root(&self) -> Option<<Self::Parameters as FieldBasedMerkleTreeParameters>::Data>;

    // Given an `index` returns the MerklePath of the leaf at that index up until the root of the
    // `self` Merkle Tree. Returns None if the tree has not been finalized before
    // calling this function.
    fn get_merkle_path(&self, leaf_index: usize) -> Option<Self::MerklePath>;
}

pub trait FieldBasedMerkleTreePath {
    type Data: Field;

    // A Merkle Path for a leaf of a Merkle Tree with arity >= 2 will be made up of couples of nodes
    // and an integer. The nodes are all the siblings of the leaf (MERKLE_ARITY - 1) and the integer
    // is the position of the leaf, among its siblings, in the input to the hash function (valid
    // values are from 0 to MERKLE_ARITY - 1).
    fn new(path: &[(Vec<Self::Data>, usize)]) -> Self;

    // Verify the Merkle Path for `leaf` given the `root` of the Merkle Tree and its `height`.
    fn verify(&self, height: usize, leaf: &Self::Data, root: &Self::Data) -> Result<bool, Error>;

    fn get_raw_path(&self) -> Vec<(Vec<Self::Data>, usize)>;
}

#[derive(Clone, Debug)]
pub struct FieldBasedMHTPath<H: FieldBasedHash, T: FieldBasedMerkleTreeParameters<Data = H::Data>>{
    pub path: Vec<(Vec<H::Data>, usize)>,
    _tree_params: PhantomData<T>
}

impl<H: FieldBasedHash, T: FieldBasedMerkleTreeParameters<Data = H::Data>> FieldBasedMHTPath<H, T> {

    pub fn compare_with_binary(&self, binary_path: &[(H::Data, bool)]) -> bool {
        if self.path.len() != binary_path.len() { return false };

        for ((p1_node, p1_pos), &(p2_node, p2_pos)) in self.path.iter().zip(binary_path.iter()) {
            if  p1_node.len() != 1 || p1_node[0] != p2_node || (*p1_pos != 0 && *p1_pos != 1) ||
                (*p1_pos == 0 && p2_pos) || (*p1_pos == 1 && !p2_pos)
            {
                return false
            }
        }
        return true;
    }
}

impl<H: FieldBasedHash, T: FieldBasedMerkleTreeParameters<Data = H::Data>> FieldBasedMerkleTreePath for FieldBasedMHTPath<H, T> {

    type Data = H::Data;

    fn new(path: &[(Vec<H::Data>, usize)]) -> Self {
        Self {path: path.to_vec(), _tree_params: PhantomData}
    }

    fn verify(&self, height: usize, leaf: &Self::Data, root: &Self::Data) -> Result<bool, Error>{
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design.
        assert_eq!(<H::Parameters as FieldBasedHashParameters>::R , T::MERKLE_ARITY);
        if self.path.len() != height - 1 {
            Err(MerkleTreeError::IncorrectPathLength(self.path.len(), height - 1))?
        }
        let mut digest = H::init(None);
        let mut prev_node = leaf.clone();
        for (sibling_nodes, position) in self.path.as_slice() {
            assert_eq!(sibling_nodes.len(), T::MERKLE_ARITY - 1);
            let prev_node_position = *position % T::MERKLE_ARITY;
            for i in 0..T::MERKLE_ARITY {
                if i == prev_node_position {
                    digest.update(prev_node.clone());
                } else {
                    digest.update(sibling_nodes[i % (T::MERKLE_ARITY - 1)]);
                }
            }
            prev_node = digest.finalize();
            digest.reset(None);
        }
        if prev_node == root.clone() {
            Ok(true)
        } else {
            Ok(false)
        }
    }

    fn get_raw_path(&self) -> Vec<(Vec<Self::Data>, usize)> {
        self.path.clone()
    }
}