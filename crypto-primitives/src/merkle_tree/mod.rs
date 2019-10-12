use crate::crh::FixedLengthCRH;
use algebra::bytes::ToBytes;
use crate::Error;
use std::{fmt, rc::Rc};


#[cfg(feature = "r1cs")]
pub mod constraints;

pub trait MerkleTreeConfig {
    const HEIGHT: usize;
    type H: FixedLengthCRH;
}

/// Stores the hashes of a particular path (in order) from leaf to root.
/// Our path `is_left_child()` if the boolean in `path` is true.
#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: MerkleTreeConfig"),
    Debug(bound = "P: MerkleTreeConfig, <P::H as FixedLengthCRH>::Output: fmt::Debug")
)]
pub struct MerkleTreePath<P: MerkleTreeConfig> {
    pub(crate) path: Vec<(<P::H as FixedLengthCRH>::Output, <P::H as FixedLengthCRH>::Output)>,
}

pub type MerkleTreeParams<P> = <<P as MerkleTreeConfig>::H as FixedLengthCRH>::Parameters;
pub type MerkleTreeDigest<P> = <<P as MerkleTreeConfig>::H as FixedLengthCRH>::Output;


impl<P: MerkleTreeConfig> Default for MerkleTreePath<P> {
    fn default() -> Self {
        let mut path = Vec::with_capacity(P::HEIGHT as usize);
        for _i in 1..P::HEIGHT as usize {
            path.push((<P::H as FixedLengthCRH>::Output::default(), <P::H as FixedLengthCRH>::Output::default()));
        }
        Self {
            path,
        }
    }
}

impl<P: MerkleTreeConfig> MerkleTreePath<P> {
    pub fn verify<L: ToBytes>(
        &self,
        parameters: &<P::H as FixedLengthCRH>::Parameters,
        root_hash: &<P::H as FixedLengthCRH>::Output,
        leaf: &L,
    ) -> Result<bool, Error> {
        if self.path.len() != (P::HEIGHT - 1) as usize {
            return Ok(false);
        }
        // Check that the given leaf matches the leaf in the membership proof.
        let mut buffer = [0u8; 128];

        if !self.path.is_empty() {
            let claimed_leaf_hash = hash_leaf::<P::H, L>(parameters, leaf, &mut buffer)?;

            // Check if leaf is one of the bottom-most siblings.
            if claimed_leaf_hash != self.path[0].0 && claimed_leaf_hash != self.path[0].1 {
                return Ok(false);
            };

            let mut prev = claimed_leaf_hash;
            // Check levels between leaf level and root.
            for &(ref hash, ref sibling_hash) in &self.path {
                // Check if the previous hash matches the correct current hash.
                if &prev != hash && &prev != sibling_hash {
                    return Ok(false);
                };
                prev = hash_inner_node::<P::H>(parameters, hash, sibling_hash, &mut buffer)?;
            }

            if root_hash != &prev {
                return Ok(false);
            }
            Ok(true)
        } else {
            Ok(false)
        }
    }
}

pub struct MerkleHashTree<P: MerkleTreeConfig> {
    tree:         Vec<<P::H as FixedLengthCRH>::Output>,
    padding_tree: Vec<(<P::H as FixedLengthCRH>::Output, <P::H as FixedLengthCRH>::Output)>,
    parameters:   Rc<<P::H as FixedLengthCRH>::Parameters>,
    root:         Option<<P::H as FixedLengthCRH>::Output>,
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

    pub fn new<L: ToBytes>(parameters: Rc<<P::H as FixedLengthCRH>::Parameters>, leaves: &[L]) -> Result<Self, Error> {
        let new_time = start_timer!(|| "MerkleTree::New");

        let last_level_size = leaves.len().next_power_of_two();
        let tree_size = 2 * last_level_size - 1;
        let tree_height = tree_height(tree_size);
        assert!(tree_height as u8 <= Self::HEIGHT);

        // Initialize the merkle tree.
        let mut tree = Vec::with_capacity(tree_size);
        let empty_hash = hash_empty::<P::H>(&parameters)?;
        for _ in 0..tree_size {
            tree.push(empty_hash.clone());
        }

        // Compute the starting indices for each level of the tree.
        let mut index = 0;
        let mut level_indices = Vec::with_capacity(tree_height);
        for _ in 0..tree_height {
            level_indices.push(index);
            index = left_child(index);
        }

        // Compute and store the hash values for each leaf.
        let last_level_index = level_indices.pop().unwrap();
        let mut buffer = [0u8; 128];
        for (i, leaf) in leaves.iter().enumerate() {
            tree[last_level_index + i] = hash_leaf::<P::H, _>(&parameters, leaf, &mut buffer)?;
        }

        // Compute the hash values for every node in the tree.
        let mut upper_bound = last_level_index;
        let mut buffer = [0u8; 128];
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
                    &mut buffer,
                )?;
            }
            upper_bound = start_index;
        }
        // Finished computing actual tree.
        // Now, we compute the dummy nodes until we hit our HEIGHT goal.
        let mut cur_height = tree_height;
        let mut padding_tree = Vec::new();
        let mut cur_hash = tree[0].clone();
        while cur_height < (Self::HEIGHT - 1) as usize {
            cur_hash = hash_inner_node::<P::H>(&parameters, &cur_hash, &empty_hash, &mut buffer)?;
            padding_tree.push((cur_hash.clone(), empty_hash.clone()));
            cur_height += 1;
        }

        let root_hash = hash_inner_node::<P::H>(&parameters, &cur_hash, &empty_hash, &mut buffer)?;

        end_timer!(new_time);

        Ok(MerkleHashTree {
            tree,
            padding_tree,
            parameters,
            root: Some(root_hash),
        })
    }

    #[inline]
    pub fn root(&self) -> <P::H as FixedLengthCRH>::Output {
        self.root.clone().unwrap()
    }

    pub fn generate_proof<L: ToBytes>(
        &self,
        index: usize,
        leaf: &L,
    ) -> Result<MerkleTreePath<P>, Error> {
        let prove_time = start_timer!(|| "MerkleTree::GenProof");
        let mut path = Vec::new();

        let mut buffer = [0u8; 128];
        let leaf_hash = hash_leaf::<P::H, _>(&self.parameters, leaf, &mut buffer)?;
        let tree_height = tree_height(self.tree.len());
        let tree_index = convert_index_to_last_level(index, tree_height);
        let empty_hash = hash_empty::<P::H>(&self.parameters)?;

        // Check that the given index corresponds to the correct leaf.
        if leaf_hash != self.tree[tree_index] {
            Err(MerkleTreeError::IncorrectLeafIndex(tree_index))?
        }

        // Iterate from the leaf up to the root, storing all intermediate hash values.
        let mut current_node = tree_index;
        while !is_root(current_node) {
            let sibling_node = sibling(current_node).unwrap();
            let (curr_hash, sibling_hash) = (
                self.tree[current_node].clone(),
                self.tree[sibling_node].clone(),
            );
            if is_left_child(current_node) {
                path.push((curr_hash, sibling_hash));
            } else {
                path.push((sibling_hash, curr_hash));
            }
            current_node = parent(current_node).unwrap();
        }

        // Store the root node. Set boolean as true for consistency with digest
        // location.
        assert!(path.len() < Self::HEIGHT as usize);
        if path.len() != (Self::HEIGHT - 1) as usize {
            path.push((self.tree[0].clone(), empty_hash));
            for &(ref hash, ref sibling_hash) in &self.padding_tree {
                path.push((hash.clone(), sibling_hash.clone()));
            }
        }
        end_timer!(prove_time);
        if path.len() != (Self::HEIGHT - 1) as usize {
            Err(MerkleTreeError::IncorrectPathLength(path.len()))?
        } else {
            Ok(MerkleTreePath {
                path,
            })
        }
    }
}

#[derive(Debug)]
pub enum MerkleTreeError {
    IncorrectLeafIndex(usize),
    IncorrectPathLength(usize),
}

impl std::fmt::Display for MerkleTreeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let msg = match self {
            MerkleTreeError::IncorrectLeafIndex(index) => format!("incorrect leaf index: {}", index),
            MerkleTreeError::IncorrectPathLength(len) => format!("incorrect path length: {}", len),
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

/// Returns the log2 value of the given number.
#[inline]
fn log2(number: usize) -> usize {
    (number as f64).log2() as usize
}

/// Returns the height of the tree, given the size of the tree.
#[inline]
fn tree_height(tree_size: usize) -> usize {
    log2(tree_size + 1)
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
    index + (1 << (tree_height - 1)) - 1
}

/// Returns the output hash, given a left and right hash value.
pub(crate) fn hash_inner_node<H: FixedLengthCRH>(
    parameters: &H::Parameters,
    left: &H::Output,
    right: &H::Output,
    buffer: &mut [u8],
) -> Result<H::Output, Error> {
    use std::io::Cursor;
    let mut writer = Cursor::new(buffer);
    // Construct left input.
    left.write(&mut writer)?;

    // Construct right input.
    right.write(&mut writer)?;

    let buffer = writer.into_inner();
    H::evaluate(parameters, &buffer[..(H::INPUT_SIZE_BITS / 8)])
}

/// Returns the hash of a leaf.
pub(crate) fn hash_leaf<H: FixedLengthCRH, L: ToBytes>(
    parameters: &H::Parameters,
    leaf: &L,
    buffer: &mut [u8],
) -> Result<H::Output, Error> {
    use std::io::Cursor;
    let mut writer = Cursor::new(buffer);
    leaf.write(&mut writer)?;

    let buffer = writer.into_inner();
    H::evaluate(parameters, &buffer[..(H::INPUT_SIZE_BITS / 8)])
}

pub(crate) fn hash_empty<H: FixedLengthCRH>(
    parameters: &H::Parameters,
) -> Result<H::Output, Error> {
    let empty_buffer = vec![0u8; H::INPUT_SIZE_BITS / 8];
    H::evaluate(parameters, &empty_buffer)
}


#[cfg(test)]
mod test {
    use crate::{crh::{pedersen::*, *}, merkle_tree::*};
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
        const HEIGHT: usize = 32;
        type H = H;
    }
    type JubJubMerkleTree = MerkleHashTree<JubJubMerkleTreeParams>;


    fn generate_merkle_tree<L: ToBytes + Clone + Eq>(leaves: &[L]) -> () {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = JubJubMerkleTree::new(crh_parameters.clone(), &leaves).unwrap();
        let root = tree.root();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());
        }
    }

    #[test]
    fn good_root_test() {
        let mut leaves = Vec::new();
        for i in 0..4u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        generate_merkle_tree(&leaves);
        let mut leaves = Vec::new();
        for i in 0..100u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        generate_merkle_tree(&leaves);
    }

    fn bad_merkle_tree_verify<L: ToBytes + Clone + Eq>(leaves: &[L]) -> () {
        use algebra::groups::Group;
        let mut rng = XorShiftRng::seed_from_u64(13423423u64);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = JubJubMerkleTree::new(crh_parameters.clone(), &leaves).unwrap();
        let root = JubJub::zero();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());
        }
    }

    #[should_panic]
    #[test]
    fn bad_root_test() {
        let mut leaves = Vec::new();
        for i in 0..4u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        generate_merkle_tree(&leaves);
        let mut leaves = Vec::new();
        for i in 0..100u8 {
            leaves.push([i, i, i, i, i, i, i, i]);
        }
        bad_merkle_tree_verify(&leaves);
    }
}
