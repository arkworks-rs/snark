use algebra::Field;
use crate::crh::FieldBasedHash;
use crate::merkle_tree::*;
use crate::Error;

/// Merkle Tree whose leaves are field elements, best with hash functions
/// that works with field elements, such as Poseidon. This implementation
/// works with leaves of size 1 field element.
/// Leaves passed when creating a MerkleTree/MerklePath proof won't be
/// hashed, it's responsibility of the caller to do it, if desired.
pub struct NaiveMerkleTree<P: FieldBasedMerkleTreeParameters> {
    tree:         Vec<<P::H as FieldBasedHash>::Data>,
    padding_tree: Vec<(
        <P::H as FieldBasedHash>::Data,
        <P::H as FieldBasedHash>::Data,
    )>,
    root:         Option<<P::H as FieldBasedHash>::Data>,
}

impl<P: FieldBasedMerkleTreeParameters> NaiveMerkleTree<P> {
    pub const HEIGHT: u8 = P::HEIGHT as u8;

    pub fn blank() -> Self {
        NaiveMerkleTree {
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

        Ok(NaiveMerkleTree {
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
    ) -> Result<FieldBasedBinaryMHTPath<P>, Error>
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
            Err(MerkleTreeError::IncorrectPathLength(path.len(), (Self::HEIGHT - 1) as usize))?
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
    Ok(H::init(None)
        .update(left)
        .update(right)
        .finalize())
}

pub(crate) fn hash_empty<H: FieldBasedHash>() -> Result<H::Data, Error> {
    Ok(<H::Data as Field>::zero())
}

#[cfg(test)]
mod test {
    use crate::{crh::{
        MNT4PoseidonHash, batched_crh::MNT4BatchPoseidonHash,
    }, merkle_tree::field_based_mht::*, FieldBasedHash};
    use algebra::{
        fields::mnt4753::Fr as MNT4753Fr, Field,
        UniformRand
    };
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    #[derive(Clone)]
    struct MNT4753FieldBasedMerkleTreeParams;
    impl FieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
        type Data = MNT4753Fr;
        type H = MNT4PoseidonHash;
        const HEIGHT: usize = 6;
        const MERKLE_ARITY: usize = 2;
        const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> = Some(MNT4753_MHT_POSEIDON_PARAMETERS);
    }

    impl BatchFieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
        type BH = MNT4BatchPoseidonHash;
    }

    type MNT4753FieldBasedMerkleTree = NaiveMerkleTree<MNT4753FieldBasedMerkleTreeParams>;
    type MNT4PoseidonMHT = FieldBasedOptimizedMHT<MNT4753FieldBasedMerkleTreeParams>;

    fn generate_merkle_tree(leaves: &[MNT4753Fr])
    {
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();
        let root = tree.root();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, leaf).unwrap();
            assert!(proof.verify(&leaf, &root).unwrap());
        }
    }

    #[test]
    fn good_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..4 {
            leaves.push(MNT4PoseidonHash::init(None).update(MNT4753Fr::rand(&mut rng)).finalize());
        }
        generate_merkle_tree(&leaves);

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for _ in 0..16 {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        generate_merkle_tree(&leaves);

        //Test #leaves == 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..32 {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        generate_merkle_tree(&leaves);
    }

    fn bad_merkle_tree_verify(leaves: &[MNT4753Fr])
    {
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();
        let root = MNT4753Fr::zero();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, leaf).unwrap();
            assert!(!proof.verify(&leaf, &root).unwrap());
        }
    }

    #[test]
    fn bad_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);

        //Test #leaves << 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..4 {
            leaves.push(MNT4PoseidonHash::init(None).update(MNT4753Fr::rand(&mut rng)).finalize());
        }
        bad_merkle_tree_verify(&leaves);

        //Test #leaves = 2^HEIGHT - 1
        let mut leaves = Vec::new();
        for _ in 0..16 {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        bad_merkle_tree_verify(&leaves);

        //Test #leaves == 2^HEIGHT
        let mut leaves = Vec::new();
        for _ in 0..32 {
            let f = MNT4753Fr::rand(&mut rng);
            leaves.push(f);
        }
        bad_merkle_tree_verify(&leaves);
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
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();
        let root1 = tree.root();

        let mut tree = MNT4PoseidonMHT::init();
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        for _ in 0..num_leaves {
            tree.append(MNT4753Fr::rand(&mut rng));
        }
        tree.finalize_in_place();
        assert_eq!(tree.root().unwrap(), root1, "Outputs of the Merkle trees for MNT4 do not match.");
    }
}