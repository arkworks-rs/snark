use crate::{
    crh::FieldBasedHash, Error
};
use super::*;
use algebra::ToConstraintField;

pub trait FieldBasedMerkleTreeConfig {
    const HEIGHT: usize;
    type H: FieldBasedHash;
}

/// Stores the hashes of a particular path (in order) from leaf to root.
/// Our path `is_left_child()` if the boolean in `path` is false.
#[derive(Derivative)]
#[derivative(
Clone(bound = "P: FieldBasedMerkleTreeConfig"),
Debug(bound = "P: FieldBasedMerkleTreeConfig, <P::H as FieldBasedHash>::Data: fmt::Debug")
)]
pub struct FieldBasedMerkleTreePath<P: FieldBasedMerkleTreeConfig> {
    pub(crate) path: Vec<(
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
    pub fn verify<L>(
        &self,
        root_hash: &<P::H as FieldBasedHash>::Data,
        leaf: &L,
    ) -> Result<bool, Error>
    where
        L: ToConstraintField<<P::H as FieldBasedHash>::Data>
    {
        if self.path.len() != (P::HEIGHT - 1) as usize {
            return Ok(false);
        }

        if !self.path.is_empty() {
            let mut prev = hash_leaf::<P::H, L>(leaf)?;

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
            Ok(false)
        }
    }
}

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

    pub fn new<L>(
        leaves: &[L],
    ) -> Result<Self, Error>
    where
        L: ToConstraintField<<P::H as FieldBasedHash>::Data>
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

        // Compute and store the hash values for each leaf.
        let last_level_index = level_indices.pop().unwrap();
        for (i, leaf) in leaves.iter().enumerate() {
            tree[last_level_index + i] = hash_leaf::<P::H, _>(leaf)?;
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
        while cur_height < (Self::HEIGHT - 1) as usize {
            cur_hash = hash_inner_node::<P::H>(cur_hash, empty_hash)?;
            padding_tree.push((cur_hash.clone(), empty_hash.clone()));
            cur_height += 1;
        }

        let root_hash = hash_inner_node::<P::H>(cur_hash, empty_hash)?;

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

    pub fn generate_proof<L>(
        &self,
        index: usize,
        leaf: &L,
    ) -> Result<FieldBasedMerkleTreePath<P>, Error>
        where
            L: ToConstraintField<<P::H as FieldBasedHash>::Data>
    {
        let prove_time = start_timer!(|| "MerkleTree::GenProof");
        let mut path = Vec::new();

        let leaf_hash = hash_leaf::<P::H, _>(leaf)?;
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

        assert!(path.len() < Self::HEIGHT as usize);
        if path.len() != (Self::HEIGHT - 1) as usize {

            //Push the sibling of the tree root which is empty hash.
            let empty_hash = hash_empty::<P::H>()?;
            path.push((empty_hash, false));

            //Then push the other elements of the padding tree
            for &(_, ref sibling_hash) in &self.padding_tree {
                path.push((sibling_hash.clone(), false));
            }
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

/// Returns the hash of a leaf.
pub(crate) fn hash_leaf<H, L>(leaf: &L) -> Result<H::Data, Error>
where
    H: FieldBasedHash,
    L: ToConstraintField<H::Data>

{
    let input = leaf.to_field_elements()?;
    H::evaluate(input.as_slice())
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
        Error,
    };
    use algebra::{
        fields::mnt4753::Fr, ToConstraintField, Field,
        curves::mnt6753::{G1Affine, G1Projective}, ProjectiveCurve,
        UniformRand
    };
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;

    #[derive(Derivative)]
    #[derivative(Clone, Eq, PartialEq)]
    struct MNT4753Utxo {
        public_key: G1Affine,
        amount:     Fr,
    }

    impl ToConstraintField<Fr> for MNT4753Utxo {
        fn to_field_elements(&self) -> Result<Vec<Fr>, Error> {
            let mut f = self.public_key.to_field_elements()?;
            f.push(self.amount);
            Ok(f)
        }
    }

    impl UniformRand for MNT4753Utxo {
        fn rand<R: Rng + ?Sized>(rng: &mut R) -> Self {
            let public_key = G1Projective::rand(rng).into_affine();
            let amount: Fr = rng.gen();
            MNT4753Utxo{public_key, amount}
        }
    }

    struct MNT4753FieldBasedMerkleTreeParams;

    impl FieldBasedMerkleTreeConfig for MNT4753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 32;
        type H = MNT4PoseidonHash;
    }

    type MNT4753FieldBasedMerkleTree = FieldBasedMerkleHashTree<MNT4753FieldBasedMerkleTreeParams>;

    fn generate_merkle_tree<L>(leaves: &[L])
    where
        L: ToConstraintField<Fr> + Clone + Eq
    {
        let tree = MNT4753FieldBasedMerkleTree::new(&leaves).unwrap();
        let root = tree.root();
        for (i, leaf) in leaves.iter().enumerate() {
            let proof = tree.generate_proof(i, leaf).unwrap();
            assert!(proof.verify(&root, leaf).unwrap());
        }
    }

    #[test]
    fn good_root_test() {
        let mut rng = XorShiftRng::seed_from_u64(9174123u64);
        let mut leaves = Vec::new();
        for _ in 0..4 {
            let f = MNT4753Utxo::rand(&mut rng);
            leaves.push(f);
        }
        generate_merkle_tree(&leaves);
        let mut leaves = Vec::new();
        for _ in 0..100 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        generate_merkle_tree(&leaves);
    }

    fn bad_merkle_tree_verify<L>(leaves: &[L])
        where
            L: ToConstraintField<Fr> + Clone + Eq
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
        let mut leaves = Vec::new();
        for _ in 0..4 {
            let f = MNT4753Utxo::rand(&mut rng);
            leaves.push(f);
        }
        generate_merkle_tree(&leaves);
        let mut leaves = Vec::new();
        for _ in 0..100 {
            let f: Fr = rng.gen();
            leaves.push(f);
        }
        bad_merkle_tree_verify(&leaves);
    }
}