use crate::{crh::FixedLengthCRH, Error};
use algebra::bytes::ToBytes;

use serde::{Deserialize, Serialize};
use std::{fmt, rc::Rc};

pub mod field_based_mht;
pub use self::field_based_mht::*;

pub trait MerkleTreeConfig {
    const HEIGHT: usize;
    type H: FixedLengthCRH;
}

/// Stores the hashes of a particular path (in order) from leaf to root.
/// Our path `is_left_child()` if the boolean in `path` is false.
#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: MerkleTreeConfig"),
    Debug(bound = "P: MerkleTreeConfig, <P::H as FixedLengthCRH>::Output: fmt::Debug")
)]
#[derive(Serialize, Deserialize)]
pub struct MerkleTreePath<P: MerkleTreeConfig> {
    pub path: Vec<(<P::H as FixedLengthCRH>::Output, bool)>,
}

pub type MerkleTreeParams<P> = <<P as MerkleTreeConfig>::H as FixedLengthCRH>::Parameters;
pub type MerkleTreeDigest<P> = <<P as MerkleTreeConfig>::H as FixedLengthCRH>::Output;

impl<P: MerkleTreeConfig> Default for MerkleTreePath<P> {
    fn default() -> Self {
        let mut path = Vec::with_capacity(P::HEIGHT as usize);
        for _i in 1..P::HEIGHT as usize {
            path.push((<P::H as FixedLengthCRH>::Output::default(), false));
        }
        Self { path }
    }
}

impl<P: MerkleTreeConfig> MerkleTreePath<P> {
    pub fn verify<L: ToBytes>(
        &self,
        parameters: &<P::H as FixedLengthCRH>::Parameters,
        root_hash: &<P::H as FixedLengthCRH>::Output,
        leaf: &L,
    ) -> Result<bool, Error> {
        if P::HEIGHT == 0 {
            Err(MerkleTreeError::Other(
                "Unable to verify: no existence proof defined for Merkle Tree of trivial height"
                    .to_owned(),
            ))?
        }

        if self.path.len() != P::HEIGHT as usize {
            return Err(MerkleTreeError::IncorrectPathLength(
                self.path.len(),
                P::HEIGHT as usize,
            ))?;
        }

        // Check that the given leaf matches the leaf in the membership proof.
        let mut prev = hash_leaf::<P::H, L>(parameters, leaf)?;

        // Check levels between leaf level and root.
        for &(ref sibling_hash, direction) in &self.path {
            // Check if the previous hash matches the correct current hash.
            prev = {
                if direction {
                    hash_inner_node::<P::H>(parameters, sibling_hash, &prev)
                } else {
                    hash_inner_node::<P::H>(parameters, &prev, sibling_hash)
                }
            }?;
        }

        if root_hash != &prev {
            return Ok(false);
        }
        Ok(true)
    }
}

/// WARNING. This Merkle Tree implementation:
/// 1) Stores all the nodes in memory, so please retain from using it if
///    the available amount of memory is limited compared to the number
///    of leaves to be stored;
/// 2) Leaves and nodes are hashed without using any kind of domain separation:
///    while this is ok for use cases where the Merkle Trees have always the
///    same height, it's not for all the others.
pub struct MerkleHashTree<P: MerkleTreeConfig> {
    tree: Vec<<P::H as FixedLengthCRH>::Output>,
    padding_tree: Vec<(
        <P::H as FixedLengthCRH>::Output,
        <P::H as FixedLengthCRH>::Output,
    )>,
    parameters: Rc<<P::H as FixedLengthCRH>::Parameters>,
    root: Option<<P::H as FixedLengthCRH>::Output>,
}

impl<P: MerkleTreeConfig> MerkleHashTree<P> {
    pub const HEIGHT: u8 = P::HEIGHT as u8;

    pub fn blank(parameters: Rc<<P::H as FixedLengthCRH>::Parameters>) -> Self {
        MerkleHashTree {
            tree: Vec::new(),
            padding_tree: Vec::new(),
            root: None,
            parameters,
        }
    }

    pub fn new<L: ToBytes>(
        parameters: Rc<<P::H as FixedLengthCRH>::Parameters>,
        leaves: &[L],
    ) -> Result<Self, Error> {
        // Deal with edge cases
        if Self::HEIGHT == 0 {
            // If height = 0, return tree with only the root
            // set to be hash_empty if leaves.is_empty() has been passed
            // or to the hash of the first (and only) leaf if !leaves.is_empty().
            let root = if leaves.is_empty() {
                hash_empty::<P::H>(&parameters)
            } else {
                if leaves.len() > 1 {
                    Err(MerkleTreeError::TooManyLeaves(Self::HEIGHT as usize))?
                }
                hash_leaf::<P::H, L>(&parameters, &leaves[0])
            }?;

            Ok(Self {
                tree: vec![root.clone()],
                padding_tree: vec![],
                parameters,
                root: Some(root),
            })
        } else {
            // Otherwise, compute root normally
            let new_time = start_timer!(|| "MerkleTree::New");
            let num_leaves = leaves.len();
            let last_level_size = if num_leaves <= 1 {
                2usize
            } else {
                num_leaves.next_power_of_two()
            };
            let tree_size = 2 * last_level_size - 1;
            let tree_height = tree_height(tree_size);
            if tree_height as u8 > Self::HEIGHT {
                Err(MerkleTreeError::TooManyLeaves(Self::HEIGHT as usize))?
            }

            // Initialize the merkle tree.
            let mut tree = Vec::with_capacity(tree_size);
            let empty_hash = hash_empty::<P::H>(&parameters)?;
            for _ in 0..tree_size {
                tree.push(empty_hash.clone());
            }

            // Compute the starting indices for each level of the tree.
            let mut index = 0;
            let mut level_indices = Vec::with_capacity(tree_height);
            for _ in 0..=tree_height {
                level_indices.push(index);
                index = left_child(index);
            }

            // Compute and store the hash values for each leaf.
            let last_level_index = level_indices.pop().unwrap();
            for (i, leaf) in leaves.iter().enumerate() {
                tree[last_level_index + i] = hash_leaf::<P::H, _>(&parameters, leaf)?;
            }

            // Compute the hash values for every node in the tree.
            let mut upper_bound = last_level_index;
            level_indices.reverse();
            for &start_index in &level_indices {
                // Iterate over the current level.
                for current_index in start_index..upper_bound {
                    let left_index = left_child(current_index);
                    let right_index = right_child(current_index);

                    // Compute Hash(left || right).
                    tree[current_index] = hash_inner_node::<P::H>(
                        &parameters,
                        &tree[left_index],
                        &tree[right_index],
                    )?;
                }
                upper_bound = start_index;
            }
            // Finished computing actual tree.
            // Now, we compute the dummy nodes until we hit our HEIGHT goal.
            let mut cur_height = tree_height + 1;
            let mut padding_tree = Vec::new();
            let mut cur_hash = tree[0].clone();
            while cur_height <= Self::HEIGHT as usize {
                cur_hash = hash_inner_node::<P::H>(&parameters, &cur_hash, &empty_hash)?;
                padding_tree.push((cur_hash.clone(), empty_hash.clone()));
                cur_height += 1;
            }

            let root_hash = cur_hash;

            end_timer!(new_time);

            Ok(MerkleHashTree {
                tree,
                padding_tree,
                parameters,
                root: Some(root_hash),
            })
        }
    }

    #[inline]
    pub fn root(&self) -> Option<<P::H as FixedLengthCRH>::Output> {
        self.root.clone()
    }

    pub fn generate_proof<L: ToBytes>(
        &self,
        index: usize,
        leaf: &L,
    ) -> Result<MerkleTreePath<P>, Error> {
        // Check that height is bigger than zero
        if P::HEIGHT == 0 {
            Err(MerkleTreeError::Other(
                "Unable to prove: no existence proof defined for Merkle Tree of trivial height"
                    .to_owned(),
            ))?
        }

        // Check that index is not bigger than num_leaves
        if index >= 1 << P::HEIGHT {
            Err(MerkleTreeError::IncorrectLeafIndex(index))?
        }

        let prove_time = start_timer!(|| "MerkleTree::GenProof");
        let mut path = Vec::new();

        let leaf_hash = hash_leaf::<P::H, _>(&self.parameters, leaf)?;
        let tree_height = tree_height(self.tree.len());
        let tree_index = convert_index_to_last_level(index, tree_height);

        // Check that the given index corresponds to the correct leaf.
        if leaf_hash != self.tree[tree_index] {
            Err(MerkleTreeError::IncorrectLeafIndex(tree_index))?
        }

        // Iterate from the leaf up to the root, storing all intermediate hash values.
        let mut current_node = tree_index;
        while !is_root(current_node) {
            let sibling_node = sibling(current_node).unwrap();
            let sibling_hash = self.tree[sibling_node].clone();
            if is_left_child(current_node) {
                path.push((sibling_hash, false));
            } else {
                path.push((sibling_hash, true));
            }
            current_node = parent(current_node).unwrap();
        }

        if path.len() > Self::HEIGHT as usize {
            Err(MerkleTreeError::IncorrectPathLength(
                path.len(),
                Self::HEIGHT as usize,
            ))?
        }

        //Push the other elements of the padding tree
        for &(_, ref sibling_hash) in &self.padding_tree {
            path.push((sibling_hash.clone(), false));
        }

        end_timer!(prove_time);
        if path.len() != Self::HEIGHT as usize {
            Err(MerkleTreeError::IncorrectPathLength(
                path.len(),
                Self::HEIGHT as usize,
            ))?
        } else {
            Ok(MerkleTreePath { path })
        }
    }
}

#[derive(Debug)]
pub enum MerkleTreeError {
    TooManyLeaves(usize),
    IncorrectLeafIndex(usize),
    IncorrectPathLength(usize, usize),
    Other(String),
}

impl std::fmt::Display for MerkleTreeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let msg = match self {
            MerkleTreeError::TooManyLeaves(height) => format!(
                "Reached maximum number of leaves for a tree of height {}",
                height
            ),
            MerkleTreeError::IncorrectLeafIndex(index) => {
                format!("incorrect leaf index: {}", index)
            }
            MerkleTreeError::IncorrectPathLength(actual_len, expected_len) => {
                format!(
                    "Incorrect path length. Expected {}, found {}",
                    expected_len, actual_len
                )
            }
            MerkleTreeError::Other(err_str) => format!("{}", err_str),
        };
        write!(f, "{}", msg)
    }
}

impl std::error::Error for MerkleTreeError {
    #[inline]
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

/// Returns the height of the tree, given the size of the tree.
#[inline]
fn tree_height(tree_size: usize) -> usize {
    if tree_size == 1 {
        return 1;
    }

    algebra::log2_floor(tree_size) as usize
}

/// Returns true iff the index represents the root.
#[inline]
fn is_root(index: usize) -> bool {
    index == 0
}

/// Returns the index of the left child, given an index.
#[inline]
fn left_child(index: usize) -> usize {
    2 * index + 1
}

/// Returns the index of the right child, given an index.
#[inline]
fn right_child(index: usize) -> usize {
    2 * index + 2
}

/// Returns the index of the sibling, given an index.
#[inline]
fn sibling(index: usize) -> Option<usize> {
    if index == 0 {
        None
    } else if is_left_child(index) {
        Some(index + 1)
    } else {
        Some(index - 1)
    }
}

/// Returns true iff the given index represents a left child.
#[inline]
fn is_left_child(index: usize) -> bool {
    index % 2 == 1
}

/// Returns the index of the parent, given an index.
#[inline]
fn parent(index: usize) -> Option<usize> {
    if index > 0 {
        Some((index - 1) >> 1)
    } else {
        None
    }
}

#[inline]
fn convert_index_to_last_level(index: usize, tree_height: usize) -> usize {
    index + (1 << tree_height) - 1
}

/// Returns the output hash, given a left and right hash value.
pub(crate) fn hash_inner_node<H: FixedLengthCRH>(
    parameters: &H::Parameters,
    left: &H::Output,
    right: &H::Output,
) -> Result<H::Output, Error> {
    use std::io::Cursor;

    let buffer = vec![0u8; H::INPUT_SIZE_BITS / 8];
    let mut writer = Cursor::new(buffer);
    // Construct left input.
    left.write(&mut writer)?;

    // Construct right input.
    right.write(&mut writer)?;

    let buffer = writer.into_inner();
    H::evaluate(parameters, &buffer[..])
}

/// Returns the hash of a leaf.
pub(crate) fn hash_leaf<H: FixedLengthCRH, L: ToBytes>(
    parameters: &H::Parameters,
    leaf: &L,
) -> Result<H::Output, Error> {
    use std::io::Cursor;
    let buffer = vec![0u8; H::INPUT_SIZE_BITS / 8];

    let mut writer = Cursor::new(buffer);
    leaf.write(&mut writer)?;

    let buffer = writer.into_inner();
    H::evaluate(parameters, &buffer[..])
}

pub(crate) fn hash_empty<H: FixedLengthCRH>(
    parameters: &H::Parameters,
) -> Result<H::Output, Error> {
    let empty_buffer = vec![0u8; H::INPUT_SIZE_BITS / 8];
    H::evaluate(parameters, &empty_buffer)
}

#[cfg(test)]
mod test {
    use crate::{
        crh::{pedersen::*, *},
        merkle_tree::*,
    };
    use algebra::curves::jubjub::JubJubAffine as JubJub;
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    #[derive(Clone)]
    pub(super) struct Window4x256;
    impl PedersenWindow for Window4x256 {
        const WINDOW_SIZE: usize = 4;
        const NUM_WINDOWS: usize = 256;
    }

    type H = PedersenCRH<JubJub, Window4x256>;

    struct JubJubMerkleTreeParams;

    impl MerkleTreeConfig for JubJubMerkleTreeParams {
        const HEIGHT: usize = 5;
        type H = H;
    }

    fn generate_merkle_tree<
        L: ToBytes + Clone + Eq,
        H: FixedLengthCRH,
        P: MerkleTreeConfig<H = H>,
    >(
        leaves: &[L],
    ) -> ()
    where
        H::Output: std::fmt::Debug,
    {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = MerkleHashTree::<P>::new(crh_parameters.clone(), &leaves).unwrap();
        let root = tree.root().unwrap();
        if P::HEIGHT > 0 {
            for (i, leaf) in leaves.iter().enumerate() {
                // Positive test
                let proof = tree.generate_proof(i, &leaf).unwrap();
                assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());

                // Negative test: pass non existing index
                assert!(tree.generate_proof(i + 1 << P::HEIGHT, &leaf).is_err());
            }
        } else {
            // Assert error when trying to create merkle proofs for trees of trivial height
            for (i, leaf) in leaves.iter().enumerate() {
                assert!(tree.generate_proof(i, &leaf).is_err());
            }

            // Check correct root value
            if leaves.len() == 0 {
                assert_eq!(root, hash_empty::<H>(&crh_parameters).unwrap());
            } else {
                assert_eq!(
                    root,
                    hash_leaf::<H, L>(&crh_parameters, &leaves[0]).unwrap()
                );
            }
        }
    }

    #[test]
    fn good_root_test() {
        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for i in 0..4u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        generate_merkle_tree::<_, _, JubJubMerkleTreeParams>(&leaves);

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for i in 0..16u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        generate_merkle_tree::<_, _, JubJubMerkleTreeParams>(&leaves);

        //Test #leaves = 2^HEIGHT
        let mut leaves = Vec::new();
        for i in 0..32u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        generate_merkle_tree::<_, _, JubJubMerkleTreeParams>(&leaves);
    }

    fn bad_merkle_tree_verify<
        L: ToBytes + Clone + Eq,
        H: FixedLengthCRH,
        P: MerkleTreeConfig<H = H>,
    >(
        leaves: &[L],
    ) -> () {
        let mut rng = XorShiftRng::seed_from_u64(13423423u64);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = MerkleHashTree::<P>::new(crh_parameters.clone(), &leaves).unwrap();
        let root = hash_empty::<H>(&crh_parameters).unwrap();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, &leaf).unwrap();
            assert!(!proof.verify(&crh_parameters, &root, &leaf).unwrap());
        }
    }

    #[test]
    fn bad_root_test() {
        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for i in 0..4u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        bad_merkle_tree_verify::<_, _, JubJubMerkleTreeParams>(&leaves);

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for i in 0..16u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        bad_merkle_tree_verify::<_, _, JubJubMerkleTreeParams>(&leaves);

        //Test #leaves = 2^HEIGHT
        let mut leaves = Vec::new();
        for i in 0..32u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        bad_merkle_tree_verify::<_, _, JubJubMerkleTreeParams>(&leaves);
    }

    // Params for Merkle Tree of height 0
    struct JubJubOnlyRootMerkleTreeParams;

    impl MerkleTreeConfig for JubJubOnlyRootMerkleTreeParams {
        const HEIGHT: usize = 0;
        type H = H;
    }

    // Params for Merkle Tree of height 1
    struct JubJubHeightOneMerkleTreeParams;

    impl MerkleTreeConfig for JubJubHeightOneMerkleTreeParams {
        const HEIGHT: usize = 1;
        type H = H;
    }

    #[test]
    fn test_edge_cases() {
        // HEIGHT > 1 with 0 or 1 leaves
        {
            // Generate empty Merkle Tree
            let mut leaves: Vec<Vec<u8>> = vec![];
            generate_merkle_tree::<_, _, JubJubMerkleTreeParams>(&leaves);

            // Generate Merkle Tree with only 1 leaf
            leaves.push(vec![1u8; 8]);
            generate_merkle_tree::<_, _, JubJubMerkleTreeParams>(&leaves);

            // Generate Merkle Tree with more leaves with respect to the specified height
            for i in 0..32u8 {
                leaves.push(vec![i, i, i, i, i, i, i, i]);
            }
            assert!(std::panic::catch_unwind(|| generate_merkle_tree::<
                _,
                _,
                JubJubMerkleTreeParams,
            >(&leaves))
            .is_err());
        }

        // HEIGHT == 1
        {
            // Generate empty Merkle Tree
            let mut leaves: Vec<Vec<u8>> = vec![];
            generate_merkle_tree::<_, _, JubJubHeightOneMerkleTreeParams>(&leaves);

            // Generate Merkle Tree with only 1 leaf
            leaves.push(vec![1u8; 8]);
            generate_merkle_tree::<_, _, JubJubHeightOneMerkleTreeParams>(&leaves);

            // Generate Merkle Tree with more leaves with respect to the specified height
            for i in 0..2u8 {
                leaves.push(vec![i, i, i, i, i, i, i, i]);
            }
            assert!(std::panic::catch_unwind(|| generate_merkle_tree::<
                _,
                _,
                JubJubHeightOneMerkleTreeParams,
            >(&leaves))
            .is_err());
        }

        // HEIGHT == 0
        {
            let mut leaves: Vec<Vec<u8>> = vec![];

            // Generate Merkle Tree with only the root, but without passing any leaf
            generate_merkle_tree::<_, _, JubJubOnlyRootMerkleTreeParams>(&leaves);

            // Generate Merkle Tree with only the root, and passing one leaf
            leaves.push(vec![1u8; 8]);
            generate_merkle_tree::<_, _, JubJubOnlyRootMerkleTreeParams>(&leaves);

            // Generate Merkle Tree with only the root, passing more than one leaf. Assert error
            leaves.push(vec![2u8; 8]);
            assert!(std::panic::catch_unwind(|| generate_merkle_tree::<
                _,
                _,
                JubJubOnlyRootMerkleTreeParams,
            >(&leaves))
            .is_err());
        }
    }
}
