use algebra::Field;
use std::fmt::Debug;

pub mod poseidon;

pub trait RandomAccessMerkleTree: Clone {
    type Data: Field;
    type MerklePath: Clone + Debug;

    // Initialize this tree. The user must pass the maximum number of leaves the tree must
    // support.
    fn init(num_leaves: usize) -> Self;

    // This function reset the pointing indices for the new element positions and processed positions
    // Mainly used for benchmarking
    fn reset(&mut self) -> &mut Self;

    // This function appends a new leaf to the Merkle tree
    fn append(&mut self, leaf: Self::Data) -> &mut Self;

    // This function sets the leaf at `index` to `new_leaf`.
    fn set(&mut self, index: usize, new_leaf: Self::Data) -> &mut Self;

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

    // Given a `leaf` returns the MerklePath of that leaf up until the root of the
    // `self` Merkle Tree. Returns None if the tree has not been finalized before
    // calling this function.
    fn get_merkle_path(&self, leaf: &Self::Data) -> Option<Self::MerklePath>;

    // Verify `path` given `leaf` and `root`. Returns None if the tree has not been
    // finalized before calling this function.
    fn verify_merkle_path(
        &self,
        leaf: &Self::Data,
        root: &Self::Data,
        path: &Self::MerklePath
    ) -> Option<bool>;
}

pub trait NullableLeaf {
    type Data: Field;

    fn get_null_value() -> Self::Data;
}

pub trait RandomAccessDeleteMerkleTree: RandomAccessMerkleTree {
    type Leaf: NullableLeaf<Data = <Self as RandomAccessMerkleTree>::Data>;

    // Remove the leaf at `index`
    fn remove(&mut self, index: usize) -> &mut Self {
        self.set(index, <Self::Leaf as NullableLeaf>::get_null_value())
    }
}