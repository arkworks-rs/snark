use crate::{BatchFieldBasedHash};

pub mod poseidon;

pub trait BatchMerkleTree:Clone {
    type Hash: BatchFieldBasedHash;

    // This function adds a leaf to the Merkle tree
    // Whenever the number of added leafs reaches the block size passed as parameter during creation of the
    // Merkle tree, the parent nodes are updated with the corresponding hashes
    fn update(&mut self, leaf: <Self::Hash as BatchFieldBasedHash>::Data);

    // This function finalizes the computation of the Merkle tree
    // If the number of added leafs is smaller than the total number of leafs of the tree,
    // the remaining nodes are considered as zero elements
    // This function return a copy of the Merkle tree
    fn finalize(&self) -> Self;

    // This function finalizes the computation of the Merkle tree
    // If the number of added leafs is smaller than the total number of leafs of the tree,
    // the remaining nodes are considered as zero elements
    // Once this function is called, it is not possible to further update the tree.
    fn finalize_in_place(&mut self);

}

