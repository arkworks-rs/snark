use crate::{
    crh::FieldBasedHash, Error
};
use super::*;

pub trait FieldBasedMerkleTreeConfig {
    const HEIGHT: usize;
    type H: FieldBasedHash;
}

use algebra::{
    biginteger::BigInteger768,
    fields::mnt4753::Fr as MNT4753Fr,
    field_new,
};

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

/// Stores the hashes of a particular path (in order) from leaf to root.
/// Our path `is_left_child()` if the boolean in `path` is false.
#[derive(Derivative)]
#[derivative(
Clone(bound = "P: FieldBasedMerkleTreeConfig"),
Debug(bound = "P: FieldBasedMerkleTreeConfig, <P::H as FieldBasedHash>::Data: fmt::Debug")
)]
pub struct FieldBasedMerkleTreePath<P: FieldBasedMerkleTreeConfig> {
    pub path: Vec<(
        <P::H as FieldBasedHash>::Data,
        bool,
    )>,
}

pub type FieldBasedMerkleTreeDigest<P> = <<P as FieldBasedMerkleTreeConfig>::H as FieldBasedHash>::Data;

impl<P: FieldBasedMerkleTreeConfig> Default for FieldBasedMerkleTreePath<P> {
    fn default() -> Self {
        let mut path = Vec::with_capacity(P::HEIGHT as usize);
        for _i in 1..P::HEIGHT as usize {
            path.push((
                <P::H as FieldBasedHash>::Data::default(),
                false,
            ));
        }
        Self { path }
    }
}

impl<P: FieldBasedMerkleTreeConfig> FieldBasedMerkleTreePath<P> {
    pub fn verify(
        &self,
        root_hash: &<P::H as FieldBasedHash>::Data,
        leaf: &<P::H as FieldBasedHash>::Data,
    ) -> Result<bool, Error>
    where
    {
        if self.path.len() != (P::HEIGHT - 1) as usize {
            return Err(MerkleTreeError::IncorrectPathLength(self.path.len()))?
        }

        if !self.path.is_empty() {
            let mut prev = *leaf;

            // Check levels between leaf level and root.
            for &(sibling_hash, direction) in &self.path {

                // Check if the previous hash matches the correct current hash.
                prev = {
                    if direction {
                        hash_inner_node::<P::H>(sibling_hash, prev)
                    } else {
                        hash_inner_node::<P::H>(prev, sibling_hash)
                    }
                }?;
            }

            if root_hash != &prev {
                return Ok(false);
            }
            Ok(true)
        } else {
            return Err(MerkleTreeError::IncorrectPathLength(0))?
        }
    }
}

/// Merkle Tree whose leaves are field elements, best with hash functions
/// that works with field elements, such as Poseidon. This implementation
/// works with leaves of size 1 field element.
/// Leaves passed when creating a MerkleTree/MerklePath proof won't be
/// hashed, it's responsibility of the caller to do it, if desired.
pub struct FieldBasedMerkleHashTree<P: FieldBasedMerkleTreeConfig> {
    tree:         Vec<<P::H as FieldBasedHash>::Data>,
    padding_tree: Vec<(
        <P::H as FieldBasedHash>::Data,
        <P::H as FieldBasedHash>::Data,
    )>,
    root:         Option<<P::H as FieldBasedHash>::Data>,
}

impl<P: FieldBasedMerkleTreeConfig> FieldBasedMerkleHashTree<P> {
    pub const HEIGHT: u8 = P::HEIGHT as u8;

    pub fn blank() -> Self {
        FieldBasedMerkleHashTree {
            tree: Vec::new(),
            padding_tree: Vec::new(),
            root: None,
        }
    }

    pub fn new(
        leaves: &[<P::H as FieldBasedHash>::Data],
    ) -> Result<Self, Error>
    {
        let new_time = start_timer!(|| "MerkleTree::New");

        let last_level_size = leaves.len().next_power_of_two();
        let tree_size = 2 * last_level_size - 1;
        let tree_height = tree_height(tree_size);
        assert!(tree_height as u8 <= Self::HEIGHT);

        // Initialize the merkle tree.
        let mut tree = Vec::with_capacity(tree_size);
        let empty_hash = hash_empty::<P::H>()?;
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

        // Compute and store the values for each leaf.
        let last_level_index = level_indices.pop().unwrap();
        for (i, leaf) in leaves.iter().enumerate() {
            tree[last_level_index + i] = *leaf;
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
                    tree[left_index],
                    tree[right_index],
                )?;
            }
            upper_bound = start_index;
        }
        // Finished computing actual tree.
        // Now, we compute the dummy nodes until we hit our HEIGHT goal.
        let mut cur_height = tree_height;
        let mut padding_tree = Vec::new();
        let mut cur_hash = tree[0].clone();
        while cur_height < Self::HEIGHT as usize {
            cur_hash = hash_inner_node::<P::H>(cur_hash, empty_hash)?;
            padding_tree.push((cur_hash.clone(), empty_hash.clone()));
            cur_height += 1;
        }

        let root_hash = cur_hash;

        end_timer!(new_time);

        Ok(FieldBasedMerkleHashTree {
            tree,
            padding_tree,
            root: Some(root_hash),
        })
    }

    #[inline]
    pub fn root(&self) -> <P::H as FieldBasedHash>::Data {
        self.root.clone().unwrap()
    }

    #[inline]
    pub fn leaves(&self) -> &[<P::H as FieldBasedHash>::Data] {
        let leaf_index = convert_index_to_last_level(0, tree_height(self.tree.len()));
        &self.tree[leaf_index..]
    }

    pub fn generate_proof(
        &self,
        index: usize,
        leaf: &<P::H as FieldBasedHash>::Data,
    ) -> Result<FieldBasedMerkleTreePath<P>, Error>
    {
        let prove_time = start_timer!(|| "MerkleTree::GenProof");
        let mut path = Vec::new();

        let tree_height = tree_height(self.tree.len());
        let tree_index = convert_index_to_last_level(index, tree_height);

        // Check that the given index corresponds to the correct leaf.
        if *leaf != self.tree[tree_index] {
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

        assert!(path.len() < Self::HEIGHT as usize);

        //Push the other elements of the padding tree
        for &(_, ref sibling_hash) in &self.padding_tree {
            path.push((sibling_hash.clone(), false));
        }

        end_timer!(prove_time);
        if path.len() != (Self::HEIGHT - 1) as usize {
            Err(MerkleTreeError::IncorrectPathLength(path.len()))?
        } else {
            Ok(FieldBasedMerkleTreePath { path })
        }
    }
}

/// Returns the output hash, given a left and right hash value.
pub(crate) fn hash_inner_node<H: FieldBasedHash>(
    left: H::Data,
    right: H::Data,
) -> Result<H::Data, Error> {
    H::evaluate(&[left, right])
}

pub(crate) fn hash_empty<H: FieldBasedHash>() -> Result<H::Data, Error> {
    use algebra::Field;
    let dummy = <H::Data as Field>::one();
    H::evaluate(&[dummy])
}

#[cfg(test)]
mod test {
    use crate::{
        crh::MNT4PoseidonHash,
        merkle_tree::field_based_mht::*,
    };
    use algebra::{
        fields::mnt4753::Fr, Field,
        UniformRand
    };
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    struct MNT4753FieldBasedMerkleTreeParams;

    impl FieldBasedMerkleTreeConfig for MNT4753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 6;
        type H = MNT4PoseidonHash;
    }

    type MNT4753FieldBasedMerkleTree = FieldBasedMerkleHashTree<MNT4753FieldBasedMerkleTreeParams>;

    fn generate_merkle_tree(leaves: &[Fr])
    {
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();
        let root = tree.root();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, leaf).unwrap();
            assert!(proof.verify(&root, leaf).unwrap());
        }
    }

    #[ignore]
    #[test]
    fn generate_mnt4753_phantom_merkle_root(){
        use algebra::{FromBytes, PrimeField, FpParameters};

        let field_size_in_bytes = (Fr::size_in_bits() + (<Fr as PrimeField>::Params::REPR_SHAVE_BITS as usize))/8;
        let magic_string = b"This represents an empty Merkle Root for a MNT4753PoseidonHash based Merkle Tree.";

        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(magic_string);
        for _ in magic_string.len()..field_size_in_bytes { hash_input.push(0u8) }
        let hash_input_f = Fr::read(hash_input.as_slice()).unwrap();

        let hash = MNT4PoseidonHash::evaluate(&[hash_input_f]).unwrap();
        assert_eq!(hash, MNT4753_PHANTOM_MERKLE_ROOT);
    }

    #[test]
    fn good_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..4 {
            let f = Fr::rand(&mut rng);
            leaves.push(MNT4PoseidonHash::evaluate(&[f]).unwrap());
        }
        generate_merkle_tree(&leaves);

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for _ in 0..16 {
            let f = Fr::rand(&mut rng);
            leaves.push(f);
        }
        generate_merkle_tree(&leaves);

        //Test #leaves == 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..32 {
            let f = Fr::rand(&mut rng);
            leaves.push(f);
        }
        generate_merkle_tree(&leaves);
    }

    fn bad_merkle_tree_verify(leaves: &[Fr])
    {
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();
        let root = Fr::zero();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, leaf).unwrap();
            assert!(!proof.verify(&root, leaf).unwrap());
        }
    }

    #[test]
    fn bad_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..4 {
            let f = Fr::rand(&mut rng);
            leaves.push(MNT4PoseidonHash::evaluate(&[f]).unwrap());
        }
        bad_merkle_tree_verify(&leaves);

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for _ in 0..16 {
            let f = Fr::rand(&mut rng);
            leaves.push(f);
        }
        bad_merkle_tree_verify(&leaves);

        //Test #leaves == 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..32 {
            let f = Fr::rand(&mut rng);
            leaves.push(f);
        }
        bad_merkle_tree_verify(&leaves);
    }
}