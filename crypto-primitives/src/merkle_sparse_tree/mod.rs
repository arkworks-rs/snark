use crate::{crh::FixedLengthCRH, Error, Vec};
use algebra_core::{bytes::ToBytes, io::Cursor};
use core::fmt;

#[cfg(not(feature = "std"))]
pub(crate) use alloc::collections::{HashMap, HashSet};

#[cfg(feature = "std")]
pub(crate) use std::collections::{HashMap, HashSet};

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
    pub(crate) path: Vec<(
        <P::H as FixedLengthCRH>::Output,
        <P::H as FixedLengthCRH>::Output,
    )>,
}

// The two-path variant is to represent Merkle tree update.
pub struct MerkleTreeTwoPaths<P: MerkleTreeConfig> {
    pub(crate) old_path: MerkleTreePath<P>,
    pub(crate) new_path: MerkleTreePath<P>,
}

pub type MerkleTreeParams<P> = <<P as MerkleTreeConfig>::H as FixedLengthCRH>::Parameters;
pub type MerkleTreeDigest<P> = <<P as MerkleTreeConfig>::H as FixedLengthCRH>::Output;

impl<P: MerkleTreeConfig> Default for MerkleTreePath<P> {
    fn default() -> Self {
        let mut path = Vec::with_capacity(P::HEIGHT as usize);
        for _i in 1..P::HEIGHT as usize {
            path.push((
                <P::H as FixedLengthCRH>::Output::default(),
                <P::H as FixedLengthCRH>::Output::default(),
            ));
        }
        Self { path }
    }
}
impl<P: MerkleTreeConfig> Default for MerkleTreeTwoPaths<P> {
    fn default() -> Self {
        let old_path: MerkleTreePath<P> = MerkleTreePath::default();
        let new_path: MerkleTreePath<P> = MerkleTreePath::default();
        Self { old_path, new_path }
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

            if claimed_leaf_hash != self.path[0].0 && claimed_leaf_hash != self.path[0].1 {
                return Ok(false);
            }

            let mut prev = claimed_leaf_hash;
            // Check levels between leaf level and root.
            for &(ref left_hash, ref right_hash) in &self.path {
                // Check if the previous hash matches the correct current hash.
                if &prev != left_hash && &prev != right_hash {
                    return Ok(false);
                }
                prev = hash_inner_node::<P::H>(parameters, left_hash, right_hash, &mut buffer)?;
            }

            if root_hash != &prev {
                return Ok(false);
            }
            Ok(true)
        } else {
            Ok(false)
        }
    }

    pub fn verify_with_index<L: ToBytes>(
        &self,
        parameters: &<P::H as FixedLengthCRH>::Parameters,
        root_hash: &<P::H as FixedLengthCRH>::Output,
        leaf: &L,
        index: usize,
    ) -> Result<bool, Error> {
        if self.path.len() != (P::HEIGHT - 1) as usize {
            return Ok(false);
        }
        // Check that the given leaf matches the leaf in the membership proof.
        let mut buffer = [0u8; 128];

        let last_level_index: usize = (1usize << (P::HEIGHT - 1) as usize) - 1;
        let tree_index: usize = last_level_index + index;

        let mut index_from_path: usize = last_level_index;
        let mut index_offset: usize = 1;

        if !self.path.is_empty() {
            let claimed_leaf_hash = hash_leaf::<P::H, L>(parameters, leaf, &mut buffer)?;

            if tree_index % 2 == 1 {
                if claimed_leaf_hash != self.path[0].0 {
                    return Ok(false);
                }
            } else {
                if claimed_leaf_hash != self.path[0].1 {
                    return Ok(false);
                }
            }

            let mut prev = claimed_leaf_hash;
            let mut prev_index = tree_index;
            // Check levels between leaf level and root.
            for &(ref left_hash, ref right_hash) in &self.path {
                // Check if the previous hash matches the correct current hash.
                if prev_index % 2 == 1 {
                    if &prev != left_hash {
                        return Ok(false);
                    }
                    index_from_path += index_offset * 0;
                } else {
                    if &prev != right_hash {
                        return Ok(false);
                    }
                    index_from_path += index_offset * 1;
                }
                index_offset *= 2;
                prev_index = (prev_index - 1) / 2;
                prev = hash_inner_node::<P::H>(parameters, left_hash, right_hash, &mut buffer)?;
            }

            if root_hash != &prev {
                return Ok(false);
            }

            if index_from_path != tree_index {
                return Ok(false);
            }

            Ok(true)
        } else {
            Ok(false)
        }
    }
}

impl<P: MerkleTreeConfig> MerkleTreeTwoPaths<P> {
    pub fn verify<L: ToBytes>(
        &self,
        parameters: &<P::H as FixedLengthCRH>::Parameters,
        old_root_hash: &<P::H as FixedLengthCRH>::Output,
        new_root_hash: &<P::H as FixedLengthCRH>::Output,
        leaf: &L,
        index: usize,
    ) -> Result<bool, Error> {
        if self.old_path.path.len() != (P::HEIGHT - 1) as usize
            || self.new_path.path.len() != (P::HEIGHT - 1) as usize
        {
            return Ok(false);
        }
        // Check that the given leaf matches the leaf in the membership proof.
        let mut buffer = [0u8; 128];

        let last_level_index: usize = (1usize << (P::HEIGHT - 1) as usize) - 1;
        let tree_index: usize = last_level_index + index;

        let mut index_from_path: usize = last_level_index;
        let mut index_offset: usize = 1;

        if !self.old_path.path.is_empty() && !self.new_path.path.is_empty() {
            // Check the new path first
            let claimed_leaf_hash = hash_leaf::<P::H, L>(parameters, leaf, &mut buffer)?;

            if tree_index % 2 == 1 {
                if claimed_leaf_hash != self.new_path.path[0].0 {
                    return Ok(false);
                }
            } else {
                if claimed_leaf_hash != self.new_path.path[0].1 {
                    return Ok(false);
                }
            }

            let mut prev = claimed_leaf_hash;
            let mut prev_index = tree_index;

            // Check levels between leaf level and root.
            for &(ref left_hash, ref right_hash) in &self.new_path.path {
                // Check if the previous hash matches the correct current hash.
                if prev_index % 2 == 1 {
                    if &prev != left_hash {
                        return Ok(false);
                    }
                    index_from_path += index_offset * 0;
                } else {
                    if &prev != right_hash {
                        return Ok(false);
                    }
                    index_from_path += index_offset * 1;
                }
                index_offset *= 2;
                prev_index = (prev_index - 1) / 2;
                prev = hash_inner_node::<P::H>(parameters, left_hash, right_hash, &mut buffer)?;
            }

            if new_root_hash != &prev {
                return Ok(false);
            }

            if index_from_path != tree_index {
                return Ok(false);
            }

            if tree_index % 2 == 1 {
                prev = self.old_path.path[0].0.clone();
            } else {
                prev = self.old_path.path[0].1.clone();
            }

            prev_index = tree_index;
            let mut new_path_iter = self.new_path.path.iter();
            for &(ref left_hash, ref right_hash) in &self.old_path.path {
                // Check if the previous hash matches the correct current hash.
                if prev_index % 2 == 1 {
                    if &prev != left_hash {
                        return Ok(false);
                    }
                } else {
                    if &prev != right_hash {
                        return Ok(false);
                    }
                }

                let new_path_corresponding_entry = new_path_iter.next();

                // Check the co-path is unchanged
                match new_path_corresponding_entry {
                    Some(x) => {
                        if prev_index % 2 == 1 {
                            if *right_hash != x.1 {
                                return Ok(false);
                            }
                        } else {
                            if *left_hash != x.0 {
                                return Ok(false);
                            }
                        }
                    },
                    None => return Ok(false),
                }

                prev_index = (prev_index - 1) / 2;
                prev = hash_inner_node::<P::H>(parameters, left_hash, right_hash, &mut buffer)?;
            }

            if old_root_hash != &prev {
                return Ok(false);
            }

            Ok(true)
        } else {
            Ok(false)
        }
    }
}

pub struct MerkleHashTree<P: MerkleTreeConfig> {
    tree:         HashMap<usize, <P::H as FixedLengthCRH>::Output>,
    parameters:   <P::H as FixedLengthCRH>::Parameters,
    root:         Option<<P::H as FixedLengthCRH>::Output>,
    empty_hashes: Vec<<P::H as FixedLengthCRH>::Output>,
}

impl<P: MerkleTreeConfig> MerkleHashTree<P> {
    pub const HEIGHT: u8 = P::HEIGHT as u8;

    pub fn blank(parameters: <P::H as FixedLengthCRH>::Parameters) -> Self {
        let empty_hashes = gen_empty_hashes::<P>(&parameters).unwrap();

        MerkleHashTree {
            tree: HashMap::new(),
            parameters,
            root: Some(empty_hashes[P::HEIGHT - 1].clone()),
            empty_hashes,
        }
    }

    pub fn new<L: ToBytes>(
        parameters: <P::H as FixedLengthCRH>::Parameters,
        leaves: &HashMap<usize, L>,
    ) -> Result<Self, Error> {
        let new_time = start_timer!(|| "MerkleTree::New");

        let last_level_size = leaves.len().next_power_of_two();
        let tree_size = 2 * last_level_size - 1;
        let tree_height = tree_height(tree_size);
        assert!(tree_height as u8 <= Self::HEIGHT);

        // Initialize the merkle tree.
        let mut tree: HashMap<usize, <P::H as FixedLengthCRH>::Output> = HashMap::new();
        let empty_hashes = gen_empty_hashes::<P>(&parameters)?;

        // Compute and store the hash values for each leaf.
        let last_level_index: usize = (1usize << (Self::HEIGHT - 1) as usize) - 1;
        let mut buffer = [0u8; 128];
        for (i, leaf) in leaves.iter() {
            tree.insert(
                last_level_index + *i,
                hash_leaf::<P::H, _>(&parameters, leaf, &mut buffer)?,
            );
        }

        let mut middle_nodes: HashSet<usize> = HashSet::new();
        for i in leaves.keys() {
            middle_nodes.insert(parent(last_level_index + *i).unwrap());
        }

        // Compute the hash values for every node in parts of the tree.
        let mut buffer = [0u8; 128];
        for level in 0..Self::HEIGHT {
            // Iterate over the current level.
            for current_index in &middle_nodes {
                let left_index = left_child(*current_index);
                let right_index = right_child(*current_index);

                let mut left_hash = empty_hashes[level as usize].clone();
                let mut right_hash = empty_hashes[level as usize].clone();

                if tree.contains_key(&left_index) {
                    match tree.get(&left_index) {
                        Some(x) => left_hash = x.clone(),
                        _ => return Err(MerkleTreeError::IncorrectTreeStructure.into()),
                    }
                }

                if tree.contains_key(&right_index) {
                    match tree.get(&right_index) {
                        Some(x) => right_hash = x.clone(),
                        _ => return Err(MerkleTreeError::IncorrectTreeStructure.into()),
                    }
                }

                // Compute Hash(left || right).
                tree.insert(
                    *current_index,
                    hash_inner_node::<P::H>(&parameters, &left_hash, &right_hash, &mut buffer)?,
                );
            }

            let tmp_middle_nodes = middle_nodes.clone();
            middle_nodes.clear();
            for i in tmp_middle_nodes {
                if !is_root(i) {
                    middle_nodes.insert(parent(i).unwrap());
                }
            }
        }

        let root_hash;
        match tree.get(&0) {
            Some(x) => root_hash = (*x).clone(),
            _ => return Err(MerkleTreeError::IncorrectTreeStructure.into()),
        }

        end_timer!(new_time);

        Ok(MerkleHashTree {
            tree,
            parameters,
            root: Some(root_hash),
            empty_hashes,
        })
    }

    #[inline]
    pub fn root(&self) -> <P::H as FixedLengthCRH>::Output {
        self.root.clone().unwrap()
    }

    pub fn generate_membership_proof(&self, index: usize) -> Result<MerkleTreePath<P>, Error> {
        let prove_time = start_timer!(|| "MerkleTree::GenProof");
        let mut path = Vec::new();

        let tree_height = Self::HEIGHT;
        let tree_index = convert_index_to_last_level(index, tree_height as usize);

        // Iterate from the leaf up to the root, storing all intermediate hash values.
        let mut current_node = tree_index;
        let mut empty_hashes_iter = self.empty_hashes.iter();
        while !is_root(current_node) {
            let sibling_node = sibling(current_node).unwrap();

            let mut current_hash = empty_hashes_iter.next().unwrap().clone();
            let mut sibling_hash = current_hash.clone();

            if self.tree.contains_key(&current_node) {
                match self.tree.get(&current_node) {
                    Some(x) => current_hash = x.clone(),
                    _ => return Err(MerkleTreeError::IncorrectTreeStructure.into()),
                }
            }

            if self.tree.contains_key(&sibling_node) {
                match self.tree.get(&sibling_node) {
                    Some(x) => sibling_hash = x.clone(),
                    _ => return Err(MerkleTreeError::IncorrectTreeStructure.into()),
                }
            }

            if is_left_child(current_node) {
                path.push((current_hash, sibling_hash));
            } else {
                path.push((sibling_hash, current_hash));
            }
            current_node = parent(current_node).unwrap();
        }

        end_timer!(prove_time);
        if path.len() != (Self::HEIGHT - 1) as usize {
            return Err(MerkleTreeError::IncorrectPathLength(path.len()).into());
        } else {
            Ok(MerkleTreePath { path })
        }
    }

    pub fn generate_proof<L: ToBytes>(
        &self,
        index: usize,
        leaf: &L,
    ) -> Result<MerkleTreePath<P>, Error> {
        let mut buffer = [0u8; 128];
        let leaf_hash = hash_leaf::<P::H, _>(&self.parameters, leaf, &mut buffer)?;
        let tree_height = Self::HEIGHT;
        let tree_index = convert_index_to_last_level(index, tree_height as usize);

        // Check that the given index corresponds to the correct leaf.
        match self.tree.get(&tree_index) {
            Some(x) => {
                if leaf_hash != *x {
                    return Err(MerkleTreeError::IncorrectTreeStructure.into());
                }
            },
            _ => (),
            // _= > allow leaf_hash to be mismatched
        };

        self.generate_membership_proof(index)
    }

    pub fn update_and_prove<L: ToBytes>(
        &mut self,
        index: usize,
        new_leaf: &L,
    ) -> Result<MerkleTreeTwoPaths<P>, Error> {
        let old_path = self.generate_membership_proof(index)?;

        let mut buffer = [0u8; 128];
        let new_leaf_hash = hash_leaf::<P::H, _>(&self.parameters, new_leaf, &mut buffer)?;

        let tree_height = Self::HEIGHT;
        let tree_index = convert_index_to_last_level(index, tree_height as usize);

        // Update the leaf and update the parents
        self.tree.insert(tree_index, new_leaf_hash);

        // Iterate from the leaf up to the root, storing all intermediate hash values.
        let mut current_node = tree_index;
        current_node = parent(current_node).unwrap();

        let mut empty_hashes_iter = self.empty_hashes.iter();
        loop {
            let left_node = left_child(current_node);
            let right_node = right_child(current_node);

            let mut left_hash = empty_hashes_iter.next().unwrap().clone();
            let mut right_hash = left_hash.clone();

            if self.tree.contains_key(&left_node) {
                match self.tree.get(&left_node) {
                    Some(x) => left_hash = x.clone(),
                    _ => return Err(MerkleTreeError::IncorrectTreeStructure.into()),
                }
            }

            if self.tree.contains_key(&right_node) {
                match self.tree.get(&right_node) {
                    Some(x) => right_hash = x.clone(),
                    _ => return Err(MerkleTreeError::IncorrectTreeStructure.into()),
                }
            }

            self.tree.insert(
                current_node,
                hash_inner_node::<P::H>(&self.parameters, &left_hash, &right_hash, &mut buffer)?,
            );

            if is_root(current_node) {
                break;
            }

            current_node = parent(current_node).unwrap();
        }

        match self.tree.get(&0) {
            Some(x) => self.root = Some((*x).clone()),
            None => return Err(MerkleTreeError::IncorrectTreeStructure.into()),
        }

        let new_path = self.generate_proof(index, new_leaf)?;

        Ok(MerkleTreeTwoPaths { old_path, new_path })
    }
}

#[derive(Debug)]
pub enum MerkleTreeError {
    IncorrectLeafIndex(usize),
    IncorrectPathLength(usize),
    IncorrectTreeStructure,
}

impl core::fmt::Display for MerkleTreeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let msg = match self {
            MerkleTreeError::IncorrectLeafIndex(index) => {
                format!("incorrect leaf index: {}", index)
            },
            MerkleTreeError::IncorrectPathLength(len) => format!("incorrect path length: {}", len),
            MerkleTreeError::IncorrectTreeStructure => format!("incorrect tree structure"),
        };
        write!(f, "{}", msg)
    }
}

#[cfg(feature = "std")]
impl std::error::Error for MerkleTreeError {
    #[inline]
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

#[cfg(not(feature = "std"))]
impl algebra_core::Error for MerkleTreeError {}

/// Returns the log2 value of the given number.
#[inline]
fn log2(number: usize) -> usize {
    algebra_core::log2(number) as usize
}

/// Returns the height of the tree, given the size of the tree.
#[inline]
fn tree_height(tree_size: usize) -> usize {
    log2(tree_size)
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
    let mut writer = Cursor::new(&mut *buffer);
    // Construct left input.
    left.write(&mut writer)?;

    // Construct right input.
    right.write(&mut writer)?;

    H::evaluate(parameters, &buffer[..(H::INPUT_SIZE_BITS / 8)])
}

/// Returns the hash of a leaf.
pub(crate) fn hash_leaf<H: FixedLengthCRH, L: ToBytes>(
    parameters: &H::Parameters,
    leaf: &L,
    buffer: &mut [u8],
) -> Result<H::Output, Error> {
    let mut writer = Cursor::new(&mut *buffer);
    leaf.write(&mut writer)?;

    H::evaluate(parameters, &buffer[..(H::INPUT_SIZE_BITS / 8)])
}

pub(crate) fn hash_empty<H: FixedLengthCRH>(
    parameters: &H::Parameters,
) -> Result<H::Output, Error> {
    let empty_buffer = vec![0u8; H::INPUT_SIZE_BITS / 8];
    H::evaluate(parameters, &empty_buffer)
}

pub(crate) fn gen_empty_hashes<P: MerkleTreeConfig>(
    parameters: &<P::H as FixedLengthCRH>::Parameters,
) -> Result<Vec<<P::H as FixedLengthCRH>::Output>, Error> {
    let mut buffer = [0u8; 128];
    let mut empty_hashes = Vec::with_capacity(P::HEIGHT);

    let mut empty_hash = hash_empty::<P::H>(&parameters)?;
    empty_hashes.push(empty_hash.clone());

    for _ in 1..=P::HEIGHT {
        empty_hash = hash_inner_node::<P::H>(&parameters, &empty_hash, &empty_hash, &mut buffer)?;
        empty_hashes.push(empty_hash.clone());
    }

    Ok(empty_hashes)
}

#[cfg(test)]
mod test {
    use crate::{
        crh::{pedersen::*, *},
        merkle_sparse_tree::*,
    };
    use algebra::{jubjub::JubJubAffine as JubJub, ToBytes, Zero};
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    #[cfg(not(feature = "std"))]
    pub(crate) use alloc::collections::HashMap;

    #[cfg(feature = "std")]
    pub(crate) use std::collections::HashMap;

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

    fn generate_merkle_tree_and_test_membership<L: ToBytes + Clone + Eq>(
        leaves: &HashMap<usize, L>,
    ) -> () {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        let crh_parameters = H::setup(&mut rng).unwrap();
        let tree = JubJubMerkleTree::new(crh_parameters.clone(), leaves).unwrap();
        let root = tree.root();
        for (i, leaf) in leaves.iter() {
            let proof = tree.generate_proof(*i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());
            assert!(proof
                .verify_with_index(&crh_parameters, &root, &leaf, *i as usize)
                .unwrap());
        }
    }

    #[test]
    fn good_root_membership_test() {
        let mut leaves: HashMap<usize, u8> = HashMap::new();
        for i in 0..4u8 {
            leaves.insert(i as usize, i);
        }
        generate_merkle_tree_and_test_membership(&leaves);
        let mut leaves: HashMap<usize, u8> = HashMap::new();
        for i in 0..100u8 {
            leaves.insert(i as usize, i);
        }
        generate_merkle_tree_and_test_membership(&leaves);
    }

    fn generate_merkle_tree_with_bad_root_and_test_membership<L: ToBytes + Clone + Eq>(
        leaves: &HashMap<usize, L>,
    ) -> () {
        let mut rng = XorShiftRng::seed_from_u64(13423423u64);

        let crh_parameters = H::setup(&mut rng).unwrap();
        let tree = JubJubMerkleTree::new(crh_parameters.clone(), leaves).unwrap();
        let root = JubJub::zero();
        for (i, leaf) in leaves.iter() {
            let proof = tree.generate_proof(*i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());
            assert!(proof
                .verify_with_index(&crh_parameters, &root, &leaf, *i as usize)
                .unwrap());
        }
    }

    #[should_panic]
    #[test]
    fn bad_root_membership_test() {
        let mut leaves: HashMap<usize, u8> = HashMap::new();
        for i in 0..100u8 {
            leaves.insert(i as usize, i);
        }
        generate_merkle_tree_with_bad_root_and_test_membership(&leaves);
    }

    fn generate_merkle_tree_and_test_update<L: ToBytes + Clone + Eq>(
        old_leaves: &HashMap<usize, L>,
        new_leaves: &HashMap<usize, L>,
    ) -> () {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        let crh_parameters = H::setup(&mut rng).unwrap();
        let mut tree = JubJubMerkleTree::new(crh_parameters.clone(), old_leaves).unwrap();
        for (i, new_leaf) in new_leaves.iter() {
            let old_root = tree.root.unwrap();
            let old_leaf_option = old_leaves.get(i);

            match old_leaf_option {
                Some(old_leaf) => {
                    let old_leaf_membership_proof = tree.generate_proof(*i, &old_leaf).unwrap();
                    let update_proof = tree.update_and_prove(*i, &new_leaf).unwrap();
                    let new_leaf_membership_proof = tree.generate_proof(*i, &new_leaf).unwrap();
                    let new_root = tree.root.unwrap();

                    assert!(old_leaf_membership_proof
                        .verify_with_index(&crh_parameters, &old_root, &old_leaf, *i as usize)
                        .unwrap());
                    assert!(
                        !(old_leaf_membership_proof
                            .verify_with_index(&crh_parameters, &new_root, &old_leaf, *i as usize)
                            .unwrap())
                    );
                    assert!(new_leaf_membership_proof
                        .verify_with_index(&crh_parameters, &new_root, &new_leaf, *i as usize)
                        .unwrap());
                    assert!(
                        !(new_leaf_membership_proof
                            .verify_with_index(&crh_parameters, &new_root, &old_leaf, *i as usize)
                            .unwrap())
                    );

                    assert!(update_proof
                        .verify(
                            &crh_parameters,
                            &old_root,
                            &new_root,
                            &new_leaf,
                            *i as usize
                        )
                        .unwrap());
                },
                None => {
                    let update_proof = tree.update_and_prove(*i, &new_leaf).unwrap();
                    let new_leaf_membership_proof = tree.generate_proof(*i, &new_leaf).unwrap();
                    let new_root = tree.root.unwrap();

                    assert!(new_leaf_membership_proof
                        .verify_with_index(&crh_parameters, &new_root, &new_leaf, *i as usize)
                        .unwrap());
                    assert!(update_proof
                        .verify(
                            &crh_parameters,
                            &old_root,
                            &new_root,
                            &new_leaf,
                            *i as usize
                        )
                        .unwrap());
                },
            }
        }
    }

    #[test]
    fn good_root_update_test() {
        let mut old_leaves: HashMap<usize, u8> = HashMap::new();
        for i in 0..10u8 {
            old_leaves.insert(i as usize, i);
        }
        let mut new_leaves: HashMap<usize, u8> = HashMap::new();
        for i in 0..20u8 {
            new_leaves.insert(i as usize, i + 1);
        }
        generate_merkle_tree_and_test_update(&old_leaves, &new_leaves);
    }
}
