use algebra::Field;
use crate::crh::FieldBasedHash;
use crate::merkle_tree::*;
use crate::Error;

/// Merkle Tree whose leaves are field elements, best with hash functions
/// that works with field elements, such as Poseidon. This implementation
/// works with leaves of size 1 field element.
/// Leaves passed when creating a MerkleTree/MerklePath proof won't be
/// hashed, it's responsibility of the caller to do it, if desired.
/// WARNING. This Merkle Tree implementation:
/// 1) Stores all the nodes in memory, so please retain from using it if
///    the available amount of memory is limited compared to the number
///    of leaves to be stored;
/// 2) Leaves and nodes are hashed without using any kind of domain separation:
///    while this is ok for use cases where the Merkle Trees have always the
///    same height, it's not for all the others.
pub struct NaiveMerkleTree<P: FieldBasedMerkleTreeParameters> {
    height:       usize,
    tree:         Vec<<P::H as FieldBasedHash>::Data>,
    padding_tree: Vec<(
        <P::H as FieldBasedHash>::Data,
        <P::H as FieldBasedHash>::Data,
    )>,
    root:         Option<<P::H as FieldBasedHash>::Data>,
}

impl<P: FieldBasedMerkleTreeParameters> NaiveMerkleTree<P> {
    
    pub fn new(height: usize) -> Self {
        NaiveMerkleTree {
            height,
            tree: Vec::new(),
            padding_tree: Vec::new(),
            root: None,
        }
    } 

    pub fn append(
        &mut self,
        leaves: &[<P::H as FieldBasedHash>::Data],
    ) -> Result<(), Error>
    {
        // Deal with edge cases
        if self.height == 0 {
            // If height = 0, return tree with only the root
            // set to be empty_node if leaves.is_empty() has been passed
            // or to the the first (and only) leaf if !leaves.is_empty().
            let root = if leaves.is_empty() {
                hash_empty::<P::H>()
            } else {
                if leaves.len() > 1 {
                    Err(MerkleTreeError::TooManyLeaves(self.height))?
                }
                Ok(leaves[0])
            }?;

            (*self).tree = vec![root.clone()];
            (*self).root = Some(root);

        } else {
            let new_time = start_timer!(|| "MerkleTree::New");
            let num_leaves = leaves.len();
            let last_level_size = if num_leaves <= 1 {
                2usize
            } else {
                num_leaves.next_power_of_two()
            };
            let tree_size = 2 * last_level_size - 1;
            let tree_height = tree_height(tree_size);
            if tree_height > self.height {
                Err(MerkleTreeError::TooManyLeaves(self.height))?
            }

            // Initialize the merkle tree.
            let mut tree = Vec::with_capacity(tree_size);
            let empty_hash = hash_empty::<P::H>()?;
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
            let mut cur_height = tree_height + 1;
            let mut padding_tree = Vec::new();
            let mut cur_hash = tree[0].clone();
            while cur_height <= self.height as usize {
                cur_hash = hash_inner_node::<P::H>(cur_hash, empty_hash)?;
                padding_tree.push((cur_hash.clone(), empty_hash.clone()));
                cur_height += 1;
            }

            let root_hash = cur_hash;

            end_timer!(new_time);

            *self = NaiveMerkleTree {
                height: self.height,
                tree,
                padding_tree,
                root: Some(root_hash),
            };
        }

        Ok(())
    }

    #[inline]
    pub fn root(&self) -> Option<<P::H as FieldBasedHash>::Data> {
        self.root.clone()
    }

    #[inline]
    pub fn leaves(&self) -> &[<P::H as FieldBasedHash>::Data] {
        if self.height == 1 {
            &self.tree[..1]
        } else {
            let leaf_index = convert_index_to_last_level(0, tree_height(self.tree.len()));
            &self.tree[leaf_index..]
        }
    }

    #[inline]
    pub fn height(&self) -> usize { self.height }

    pub fn generate_proof(
        &self,
        index: usize,
        leaf: &<P::H as FieldBasedHash>::Data,
    ) -> Result<FieldBasedBinaryMHTPath<P>, Error>
    {
        // Check that height is bigger than one
        if self.height == 0 {
            Err(MerkleTreeError::Other("Unable to prove: no existence proof defined for Merkle Tree of trivial height".to_owned()))?
        }

        // Check that index is not bigger than num_leaves
        if index >= 1 << self.height {
            Err(MerkleTreeError::IncorrectLeafIndex(index))?
        }

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

        if path.len() > self.height {
            Err(MerkleTreeError::IncorrectPathLength(path.len(), self.height))?
        }

        //Push the other elements of the padding tree
        for &(_, ref sibling_hash) in &self.padding_tree {
            path.push((sibling_hash.clone(), false));
        }

        end_timer!(prove_time);
        if path.len() != self.height as usize {
            Err(MerkleTreeError::IncorrectPathLength(path.len(), self.height))?
        } else {
            Ok(FieldBasedBinaryMHTPath::<P>::new(path))
        }
    }
}

/// Returns the output hash, given a left and right hash value.
pub(crate) fn hash_inner_node<H: FieldBasedHash>(
    left: H::Data,
    right: H::Data,
) -> Result<H::Data, Error> {
    H::init_constant_length(2, None)
        .update(left)
        .update(right)
        .finalize()
}

pub(crate) fn hash_empty<H: FieldBasedHash>() -> Result<H::Data, Error> {
    Ok(<H::Data as Field>::zero())
}

#[cfg(test)]
mod test {
    use crate::{
        crh::parameters::{MNT4PoseidonHash, MNT4BatchPoseidonHash},
        merkle_tree::field_based_mht::*, FieldBasedHash
    };
    use algebra::{
        fields::mnt4753::Fr as MNT4753Fr, Field,
        UniformRand, ToBytes, to_bytes, FromBytes,
    };
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    const TEST_HEIGHT: usize = 5;

    #[derive(Debug, Clone)]
    struct MNT4753FieldBasedMerkleTreeParams;
    impl FieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
        type Data = MNT4753Fr;
        type H = MNT4PoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const ZERO_NODE_CST: Option<FieldBasedMerkleTreePrecomputedZeroConstants<'static, Self::H>> =
            Some(MNT4753_MHT_POSEIDON_PARAMETERS);
    }

    impl BatchFieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
        type BH = MNT4BatchPoseidonHash;
    }

    type MNT4753FieldBasedMerkleTree = NaiveMerkleTree<MNT4753FieldBasedMerkleTreeParams>;
    type MNT4PoseidonMHT = FieldBasedOptimizedMHT<MNT4753FieldBasedMerkleTreeParams>;

    fn generate_merkle_tree<P: FieldBasedMerkleTreeParameters>(leaves: &[P::Data], height: usize)
    {
        let mut tree = NaiveMerkleTree::<P>::new(height);
        tree.append(&leaves).unwrap();
        let root = tree.root().unwrap();
        if tree.height() > 0 {
            for (i, leaf) in leaves.iter().enumerate() {
                let proof = tree.generate_proof(i, leaf).unwrap();
                assert!(proof.verify(tree.height(), &leaf, &root).unwrap());

                // Negative test: pass non existing index
                assert!(tree.generate_proof(i + 1 << tree.height(), leaf).is_err());

                // Check leaf index is the correct one
                assert_eq!(i, proof.leaf_index());

                if i == 0 { assert!(proof.is_leftmost()); } // leftmost check
                else if i == 2usize.pow(height as u32) - 1 { assert!(proof.is_rightmost()) }  //rightmost check
                else { assert!(!proof.is_leftmost()); assert!(!proof.is_rightmost()); } // other cases check

                // Serialization/deserialization test
                let proof_serialized = to_bytes!(proof).unwrap();
                let proof_deserialized = FieldBasedBinaryMHTPath::<P>::read(proof_serialized.as_slice()).unwrap();
                assert_eq!(proof, proof_deserialized);
            }
        } else {
            // Assert error when trying to create merkle proofs for trees of trivial height
            for (i, leaf) in leaves.iter().enumerate() {
                assert!(tree.generate_proof(i, leaf).is_err());
            }

            // Check correct root value
            if leaves.len() == 0 {
                assert_eq!(root, P::Data::zero());
            } else {
                assert_eq!(root, leaves[0]);
            }
        }
    }

    #[test]
    fn good_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..4 {
            leaves.push(
                MNT4PoseidonHash::init_constant_length(
                    1,
                    None
                )
                .update(MNT4753Fr::rand(&mut rng))
                .finalize()
                .unwrap()
            );
        }
        generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, TEST_HEIGHT);

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for _ in 0..16 {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, TEST_HEIGHT);

        //Test #leaves == 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..32 {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, TEST_HEIGHT);
    }

    fn bad_merkle_tree_verify<P: FieldBasedMerkleTreeParameters>(leaves: &[P::Data], height: usize)
    {
        let mut tree = NaiveMerkleTree::<P>::new(height);
        tree.append(&leaves).unwrap();
        let root = <P::Data as Field>::zero();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, leaf).unwrap();
            assert!(!proof.verify(tree.height(), &leaf, &root).unwrap());

            // Check leaf index is the correct one
            assert_eq!(i, proof.leaf_index());

            if i == 0 { assert!(proof.is_leftmost()); } // leftmost check
            else if i == 2usize.pow(height as u32) - 1 { assert!(proof.is_rightmost()) }  //rightmost check
            else { assert!(!proof.is_leftmost()); assert!(!proof.is_rightmost()); } // other cases check

            // Serialization/deserialization test
            let proof_serialized = to_bytes!(proof).unwrap();
            let proof_deserialized = FieldBasedBinaryMHTPath::<P>::read(proof_serialized.as_slice()).unwrap();
            assert_eq!(proof, proof_deserialized);
        }
    }

    #[test]
    fn bad_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..4 {
            leaves.push(
                MNT4PoseidonHash::init_constant_length(
                    1,
                    None
                )
                    .update(MNT4753Fr::rand(&mut rng))
                    .finalize()
                    .unwrap()
            );
        }
        bad_merkle_tree_verify::<MNT4753FieldBasedMerkleTreeParams>(&leaves, TEST_HEIGHT);

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for _ in 0..16 {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        bad_merkle_tree_verify::<MNT4753FieldBasedMerkleTreeParams>(&leaves, TEST_HEIGHT);

        //Test #leaves == 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..32 {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        bad_merkle_tree_verify::<MNT4753FieldBasedMerkleTreeParams>(&leaves, TEST_HEIGHT);
    }

    #[test]
    fn compare_merkle_trees_mnt4() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        let num_leaves = 30;

        let mut leaves = Vec::new();
        for _ in 0..num_leaves {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        let mut tree = MNT4753FieldBasedMerkleTree::new(TEST_HEIGHT);
        tree.append(&leaves).unwrap();
        let root1 = tree.root().unwrap();

        let mut tree = MNT4PoseidonMHT::init(TEST_HEIGHT, num_leaves).unwrap();
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        for _ in 0..num_leaves {
            tree.append(MNT4753Fr::rand(&mut rng)).unwrap();
        }
        tree.finalize_in_place().unwrap();
        assert_eq!(tree.root().unwrap(), root1, "Outputs of the Merkle trees for MNT4 do not match.");
    }

    #[test]
    fn test_edge_cases() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        // HEIGHT > 1 with 0 or 1 leaves
        {
            // Generate empty Merkle Tree
            let mut leaves: Vec<MNT4753Fr> = vec![];
            generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, TEST_HEIGHT);

            // Generate Merkle Tree with only 1 leaf
            leaves.push(MNT4753Fr::rand(&mut rng));
            generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, TEST_HEIGHT);

            // Generate Merkle Tree with more leaves with respect to the specified height
            for _ in 0..1 << TEST_HEIGHT {
                leaves.push(MNT4753Fr::rand(&mut rng));
            }
            assert!(std::panic::catch_unwind(|| generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, TEST_HEIGHT)).is_err());
        }

        // HEIGHT == 1
        {
            // Generate empty Merkle Tree
            let mut leaves: Vec<MNT4753Fr> = vec![];
            generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, 1);

            // Generate Merkle Tree with only 1 leaf
            leaves.push(MNT4753Fr::rand(&mut rng));
            generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, 1);

            // Generate Merkle Tree with more leaves with respect to the specified height
            for _ in 0..1 << TEST_HEIGHT {
                leaves.push(MNT4753Fr::rand(&mut rng));
            }
            assert!(std::panic::catch_unwind(|| generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, 1)).is_err());
        }

        // HEIGHT == 0
        {
            let mut leaves: Vec<MNT4753Fr> = vec![];

            // Generate Merkle Tree with only the root, but without passing any leaf
            generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, 0);

            // Generate Merkle Tree with only the root, and passing one leaf
            leaves.push(MNT4753Fr::rand(&mut rng));
            generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, 0);

            // Generate Merkle Tree with only the root, passing more than one leaf. Assert error
            leaves.push(MNT4753Fr::rand(&mut rng));
            assert!(std::panic::catch_unwind(|| generate_merkle_tree::<MNT4753FieldBasedMerkleTreeParams>(&leaves, 0)).is_err());
        }
    }
}