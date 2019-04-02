use crate::{config::MAX_MERKLE_TREE_HEIGHT, crypto_primitives::crh::FixedLengthCRH};
use algebra::bytes::ToBytes;
use failure::{format_err, Error};
use std::{fmt, marker::PhantomData, rc::Rc};

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

/// Stores the hashes of a particular path (in order) from leaf to root.
/// Our path `is_left_child()` if the boolean in `path` is true.
#[derive(Derivative)]
#[derivative(
    Clone(bound = "H: FixedLengthCRH, L: ToBytes + Eq"),
    Debug(bound = "H: FixedLengthCRH, L: ToBytes + Eq, H::Output: fmt::Debug")
)]
pub struct HashMembershipProof<H: FixedLengthCRH, L: ToBytes + Eq> {
    pub(crate) path: Vec<(H::Output, H::Output)>,
    _leaf:           PhantomData<L>,
}

impl<H: FixedLengthCRH, L: ToBytes + Eq> Default for HashMembershipProof<H, L> {
    fn default() -> Self {
        let mut path = Vec::with_capacity(Self::MAX_HEIGHT as usize);
        for _i in 1..(Self::MAX_HEIGHT as usize) {
            path.push((H::Output::default(), H::Output::default()));
        }
        Self {
            path,
            _leaf: PhantomData,
        }
    }
}

impl<H: FixedLengthCRH, L: ToBytes + Eq> HashMembershipProof<H, L> {
    pub const MAX_HEIGHT: u8 = MAX_MERKLE_TREE_HEIGHT;

    pub fn verify(
        &self,
        parameters: &H::Parameters,
        root_hash: &H::Output,
        leaf: &L,
    ) -> Result<bool, Error> {
        if self.path.len() != (Self::MAX_HEIGHT - 1) as usize {
            return Ok(false);
        }
        // Check that the given leaf matches the leaf in the membership proof.
        let mut buffer = [0u8; 128];

        if !self.path.is_empty() {
            let claimed_leaf_hash = hash_leaf::<H, _>(parameters, leaf, &mut buffer)?;

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
                prev = hash_inner_node::<H>(parameters, hash, sibling_hash, &mut buffer)?;
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

pub struct MerkleHashTree<H: FixedLengthCRH, L: ToBytes + Eq> {
    tree:         Vec<H::Output>,
    padding_tree: Vec<(H::Output, H::Output)>,
    parameters:   Rc<H::Parameters>,
    _leaf:        PhantomData<L>,
    root:         Option<H::Output>,
}

impl<H: FixedLengthCRH, L: ToBytes + Eq + Clone> MerkleHashTree<H, L> {
    pub const MAX_HEIGHT: u8 = MAX_MERKLE_TREE_HEIGHT;

    pub fn blank(parameters: Rc<H::Parameters>) -> Self {
        MerkleHashTree {
            tree: Vec::new(),
            padding_tree: Vec::new(),
            root: None,
            _leaf: PhantomData,
            parameters,
        }
    }

    pub fn new(parameters: Rc<H::Parameters>, leaves: &[L]) -> Result<Self, Error> {
        let new_time = timer_start!(|| "MHT::New");

        let last_level_size = leaves.len().next_power_of_two();
        let tree_size = 2 * last_level_size - 1;
        let tree_height = tree_height(tree_size);
        assert!(tree_height as u8 <= Self::MAX_HEIGHT);

        // Initialize the merkle tree.
        let mut tree = Vec::with_capacity(tree_size);
        let empty_hash = hash_empty::<H>(&parameters)?;
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
            tree[last_level_index + i] = hash_leaf::<H, L>(&parameters, leaf, &mut buffer)?;
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
                tree[current_index] = hash_inner_node::<H>(
                    &parameters,
                    &tree[left_index],
                    &tree[right_index],
                    &mut buffer,
                )?;
            }
            upper_bound = start_index;
        }
        // Finished computing actual tree.
        // Now, we compute the dummy nodes until we hit our MAX_HEIGHT goal.
        let mut cur_height = tree_height;
        let mut padding_tree = Vec::new();
        let mut cur_hash = tree[0].clone();
        while cur_height < (Self::MAX_HEIGHT - 1) as usize {
            cur_hash = hash_inner_node::<H>(&parameters, &cur_hash, &empty_hash, &mut buffer)?;
            padding_tree.push((cur_hash.clone(), empty_hash.clone()));
            cur_height += 1;
        }

        let root_hash = hash_inner_node::<H>(&parameters, &cur_hash, &empty_hash, &mut buffer)?;

        timer_end!(new_time);

        Ok(MerkleHashTree {
            tree,
            padding_tree,
            parameters,
            _leaf: PhantomData,
            root: Some(root_hash),
        })
    }

    #[inline]
    pub fn root(&self) -> H::Output {
        self.root.clone().unwrap()
    }

    pub fn generate_proof(
        &self,
        index: usize,
        leaf: &L,
    ) -> Result<HashMembershipProof<H, L>, Error> {
        let prove_time = timer_start!(|| "MHT::GenProof");
        let mut path = Vec::new();

        let mut buffer = [0u8; 128];
        let leaf_hash = hash_leaf::<H, L>(&self.parameters, leaf, &mut buffer)?;
        let tree_height = tree_height(self.tree.len());
        let tree_index = convert_index_to_last_level(index, tree_height);
        let empty_hash = hash_empty::<H>(&self.parameters)?;

        // Check that the given index corresponds to the correct leaf.
        if leaf_hash != self.tree[tree_index] {
            return Err(format_err!("incorrect leaf index."));
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
        assert!(path.len() < Self::MAX_HEIGHT as usize);
        if path.len() != (Self::MAX_HEIGHT - 1) as usize {
            path.push((self.tree[0].clone(), empty_hash));
            for &(ref hash, ref sibling_hash) in &self.padding_tree {
                path.push((hash.clone(), sibling_hash.clone()));
            }
        }
        timer_end!(prove_time);
        if path.len() != (Self::MAX_HEIGHT - 1) as usize {
            return Err(format_err!("incorrect path length."));
        } else {
            Ok(HashMembershipProof {
                path,
                _leaf: PhantomData,
            })
        }
    }
}

#[cfg(test)]
mod test {
    use crate::crypto_primitives::{
        crh::{pedersen::*, *},
        mht::*,
    };
    use algebra::curves::jubjub::JubJubAffine as JubJub;
    use rand::{ChaChaRng, SeedableRng};

    #[derive(Clone)]
    pub(super) struct Window4x256;
    impl PedersenWindow for Window4x256 {
        const WINDOW_SIZE: usize = 4;
        const NUM_WINDOWS: usize = 256;
    }

    type H = PedersenCRH<JubJub, Window4x256>;
    type JubJubMHT<L> = MerkleHashTree<H, L>;

    fn generate_merkle_tree<L: ToBytes + Clone + Eq>(leaves: &[L]) -> () {
        let seed: [u32; 8] = [
            2053759276, 152413135, 1690980041, 4293109333, 2390175708, 686052238, 1844363894,
            1379683288,
        ];
        let mut rng = ChaChaRng::from_seed(&seed);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = JubJubMHT::<L>::new(crh_parameters.clone(), &leaves).unwrap();
        let root = tree.root();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());
        }
    }

    #[test]
    fn mht_test() {
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
        let seed: [u32; 8] = [
            2053759276, 152413135, 1690980041, 4293109333, 2390175708, 686052238, 1844363894,
            1379683288,
        ];
        let mut rng = ChaChaRng::from_seed(&seed);

        let crh_parameters = Rc::new(H::setup(&mut rng).unwrap());
        let tree = JubJubMHT::<L>::new(crh_parameters.clone(), &leaves).unwrap();
        let root = JubJub::zero();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, &leaf).unwrap();
            assert!(proof.verify(&crh_parameters, &root, &leaf).unwrap());
        }
    }

    #[should_panic]
    #[test]
    fn bad_root_mht_test() {
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
