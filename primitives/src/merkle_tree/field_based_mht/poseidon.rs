extern crate rand;

use algebra::{PrimeField, MulShort};
use crate::{crh::{
    BatchFieldBasedHash,
    poseidon::{PoseidonParameters, batched_crh::PoseidonBatchHash}
}, merkle_tree::field_based_mht::FieldBasedMerkleTree, PoseidonHash, FieldBasedHash};
use std::{marker::PhantomData, clone::Clone};

#[derive(Clone, Debug)]
pub struct PoseidonMerklePath<F: PrimeField + MulShort>(
    pub Vec<(F, bool)>
);

#[derive(Clone)]
pub struct PoseidonMerkleTree<
    F: PrimeField + MulShort,
    T: FieldBasedMerkleTreeParameters<Data = F>,
    P: PoseidonParameters<Fr = F>
>{
    root: F,
    // Stores all MT nodes
    array_nodes: Vec<F>,
    processing_step: usize,
    // Stores the initial index of each level of the MT
    initial_pos: Vec<usize>,
    // Stores the final index of each level of the MT
    final_pos: Vec<usize>,
    // Stores the index up until the nodes were already hashed, for each level
    processed_pos: Vec<usize>,
    // Stores the next available index for each level of the MT
    new_elem_pos: Vec<usize>,
    levels: usize,
    rate: usize,
    num_leaves: usize,
    finalized: bool,
    _tree_parameters: PhantomData<T>,
    _hash_parameters: PhantomData<P>,
}

impl<
    F: PrimeField + MulShort,
    T: FieldBasedMerkleTreeParameters<Data = F>,
    P: PoseidonParameters<Fr = F>
> PoseidonMerkleTree<F, T, P> {

    pub fn compute_subtree(&mut self, starting_index: usize) {
        let mut index = starting_index;
        for i in 0..self.levels  {
            if (self.new_elem_pos[i] - self.processed_pos[i]) >= self.rate {
                let num_groups_leaves = (self.new_elem_pos[i] - self.processed_pos[i]) / self.rate;
                let last_pos_to_process = self.processed_pos[i] + num_groups_leaves * self.rate;

                let (input_vec, output_vec) =
                    self.array_nodes[self.initial_pos[i]..self.final_pos[i + 1]].split_at_mut(self.final_pos[i] - self.initial_pos[i]);

                let new_pos_parent = self.new_elem_pos[i + 1] + num_groups_leaves;

                Self::batch_hash(
                    &mut input_vec[(self.processed_pos[i] - self.initial_pos[i])..(last_pos_to_process - self.initial_pos[i])],
                   &mut output_vec[(self.new_elem_pos[i + 1] - self.initial_pos[i + 1])..(new_pos_parent - self.initial_pos[i + 1])],
                    i + 1,
                    index,
                );

                self.new_elem_pos[i + 1] += num_groups_leaves;
                self.processed_pos[i] += num_groups_leaves * self.rate;
                index += self.num_leaves + (index / P::R);
                // Get the parent index, i.e. the node from which certainly all the followings are empty
            }
        }
    }

    fn batch_hash(input: &mut [F], output: &mut [F], parent_level: usize, starting_index: usize) {
        //let starting_index = Self::get_last_non_empty_index(input, &T::EMPTY_HASH_CST[parent_level - 1]);

        let starting_index = if starting_index != input.len() {
            let mut temp_index = starting_index;
            // If it's not the first children, we must go to the next chunk
            if temp_index % P::R != 0 {
                temp_index += P::R - (temp_index % P::R);
            }
            // If it was part of the last chunk, then all chunks are non-empty
            if temp_index >= input.len() { None } else { Some(temp_index) }
        } else {
            None // All chunks are non-empty
        };

        if starting_index.is_some() {
            let starting_index = starting_index.unwrap();
            let (input_to_process, _) = input.split_at_mut(starting_index);
            let (output_to_process, output_precomputed) = output.split_at_mut(starting_index / P::R);
            // Must compute the hash for the non-empty chunks
            if input_to_process.len() != 0 {
                PoseidonBatchHash::<F, P>::batch_evaluate_in_place(input_to_process, output_to_process);
            }
            // For the empty ones, instead, we can get the precomputed
            output_precomputed.iter_mut().for_each(|out| *out = T::EMPTY_HASH_CST[parent_level]);
        }
        else {
            PoseidonBatchHash::<F, P>::batch_evaluate_in_place(input, output);
        }
    }
}

impl<
    F: PrimeField + MulShort,
    T: FieldBasedMerkleTreeParameters<Data = F>,
    P: PoseidonParameters<Fr = F>
> FieldBasedMerkleTree for PoseidonMerkleTree<F, T, P> {

    type Data = F;
    type MerklePath = PoseidonMerklePath<F>;
    type Parameters = T;

    fn init(num_leaves: usize) -> Self {
        let rate = P::R;
        let processing_step = num_leaves;

        let mut initial_pos = Vec::new();
        let mut final_pos = Vec::new();
        let mut processed_pos = Vec::new();
        let mut new_elem_pos = Vec::new();

        let last_level_size = num_leaves.next_power_of_two();
        let mut size = last_level_size;

        let mut initial_idx = 0;
        let mut final_idx = last_level_size;

        let mut level_idx = 1;

        while size >= 1 {
            initial_pos.push(initial_idx);
            final_pos.push(final_idx);
            processed_pos.push(initial_idx);
            new_elem_pos.push(initial_idx);

            initial_idx += size;
            size /= rate;
            final_idx = initial_idx + size;
            level_idx += 1;
        }

        let tree_size = *final_pos.last().unwrap();

        let mut array_nodes = Vec::with_capacity(tree_size);
        for _i in 0..tree_size {
            array_nodes.push(F::zero());
        }

        let cpus = rayon::current_num_threads();
        let mut chunk_size = processing_step / (cpus * rate);
        let mut processing_block_size = chunk_size * cpus * rate;
        if processing_step < cpus * rate {
            chunk_size = processing_step / rate;
            if chunk_size == 0 {
                chunk_size = 1;
            }
            processing_block_size = chunk_size * rate;
        }

        if processing_block_size > last_level_size {
            processing_block_size = last_level_size;
        }

        Self {
            root: { F::zero() },
            array_nodes: { array_nodes },
            processing_step: { processing_block_size },
            initial_pos: { initial_pos },
            final_pos: { final_pos },
            processed_pos: { processed_pos },
            new_elem_pos: { new_elem_pos },
            levels: { level_idx - 2 },
            rate: { P::R },
            num_leaves: {last_level_size},
            finalized: false,
            _tree_parameters: PhantomData,
            _hash_parameters: PhantomData
        }
    }

    fn reset(&mut self) -> &mut Self {

        for i in 0..self.new_elem_pos.len() {
            self.new_elem_pos[i] = self.initial_pos[i];
            self.processed_pos[i] = self.initial_pos[i];
        }
        self.finalized = false;

        self
    }

    // Note: `Field` implements the `Copy` trait, therefore invoking this function won't
    // cause a moving of ownership for `leaf`, but just a copy. Another copy is
    // performed below in `self.array_nodes[self.new_elem_pos_subarray[0]] = leaf;`
    // We can reduce this to one copy by passing a reference to leaf, but from an
    // interface point of view this is not logically correct: someone calling this
    // functions will likely not use the `leaf` anymore in most of the cases
    // (in the other cases he can just clone it).
    fn append(&mut self, leaf: Self::Data) -> &mut Self {
        if self.new_elem_pos[0] < self.final_pos[0] {
            self.array_nodes[self.new_elem_pos[0]] = leaf;
            self.new_elem_pos[0] += 1;
        }

        if self.new_elem_pos[0] == self.final_pos[0] {
            self.compute_subtree(self.new_elem_pos[0]);
        }

        if (self.new_elem_pos[0] - self.processed_pos[0]) >= self.processing_step {
            self.compute_subtree(self.new_elem_pos[0]);
        }
        self
    }

    fn finalize(&self) -> Self {
        let mut copy = (*self).clone();
        let new_elem_pos_copy = copy.new_elem_pos[0];
        copy.new_elem_pos[0] = copy.final_pos[0];
        copy.compute_subtree(new_elem_pos_copy);
        copy.finalized = true;
        copy.root = *copy.array_nodes.last().unwrap();
        copy
    }

    fn finalize_in_place(&mut self) -> &mut Self {
        let new_elem_pos_copy = self.new_elem_pos[0];
        self.new_elem_pos[0] = self.final_pos[0];
        self.compute_subtree(new_elem_pos_copy);
        self.finalized = true;
        self.root = *self.array_nodes.last().unwrap();
        self
    }

    fn root(&self) -> Option<Self::Data> {
        match self.finalized {
            true => Some(self.root),
            false => None
        }
    }

    fn get_merkle_path(&self, leaf_index: usize) -> Option<Self::MerklePath> {
        match self.finalized {
            true => {
                let mut merkle_path = Vec::with_capacity(self.levels);

                let mut node_index = leaf_index;
                for _ in 0..self.levels {
                    if node_index % P::R == 0 { // Node is a left child
                        merkle_path.push((self.array_nodes[node_index + 1], false));
                    } else { // Node is a right child
                        merkle_path.push((self.array_nodes[node_index - 1], true));
                    }
                    node_index = self.num_leaves + (node_index/P::R); // Get parent index
                }
                assert_eq!(self.array_nodes[node_index], self.root); // Sanity check
                Some(PoseidonMerklePath(merkle_path))
            },
            false => None,
        }
    }

    fn verify_merkle_path(
        &self,
        leaf: &Self::Data,
        path: &Self::MerklePath
    ) -> Option<bool> {
        match self.finalized {
            true => {
                assert_eq!(path.0.len(), self.levels);
                let mut digest = PoseidonHash::<F, P>::init(None);
                let mut sibling = leaf.clone();
                let mut level = 0;
                for &(node, direction) in &path.0 {
                    if sibling == node && node == T::EMPTY_HASH_CST[level] {
                        sibling = T::EMPTY_HASH_CST[level + 1]
                    } else {
                        if direction {
                            digest.update(node).update(sibling);
                        } else {
                            digest.update(sibling).update(node);
                        }
                        sibling = digest.finalize();
                        digest.reset(None);
                    }
                    level += 1;
                }
                if sibling == self.root {
                    Some(true)
                } else {
                    Some(false)
                }
            }
            false => None,
        }
    }
}


#[cfg(test)]
mod test {
    use rand_xorshift::XorShiftRng;
    use super::rand::SeedableRng;
    use crate::crh::poseidon::parameters::{
        MNT4753PoseidonParameters, MNT6753PoseidonParameters,
    };
    use algebra::fields::mnt4753::Fr as MNT4753Fr;
    use algebra::fields::mnt6753::Fr as MNT6753Fr;
    use algebra::UniformRand;
    use algebra::biginteger::BigInteger768;
    use algebra::fields::Field;

    use crate::merkle_tree::field_based_mht::poseidon::PoseidonMerkleTree;
    use crate::merkle_tree::field_based_mht::FieldBasedMerkleTree;
    use super::{MNT4753MHTPoseidonParameters, MNT6753MHTPoseidonParameters};
    use crate::{NaiveMerkleTree, FieldBasedMerkleTreeConfig};
    use crate::crh::poseidon::{MNT4PoseidonHash, MNT6PoseidonHash};

    struct MNT4753FieldBasedMerkleTreeParams;
    impl FieldBasedMerkleTreeConfig for MNT4753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 7;
        type H = MNT4PoseidonHash;
    }
    type NaiveMNT4PoseidonMHT = NaiveMerkleTree<MNT4753FieldBasedMerkleTreeParams>;

    struct MNT6753FieldBasedMerkleTreeParams;
    impl FieldBasedMerkleTreeConfig for MNT6753FieldBasedMerkleTreeParams {
        const HEIGHT: usize = 7;
        type H = MNT6PoseidonHash;
    }
    type NaiveMNT6PoseidonMHT = NaiveMerkleTree<MNT6753FieldBasedMerkleTreeParams>;

    type MNT4PoseidonMHT = PoseidonMerkleTree<MNT4753Fr, MNT4753MHTPoseidonParameters, MNT4753PoseidonParameters>;
    type MNT6PoseidonMHT = PoseidonMerkleTree<MNT6753Fr, MNT6753MHTPoseidonParameters, MNT6753PoseidonParameters>;

    #[test]
    fn merkle_tree_test_mnt4() {
        // running time for 1048576 leaves
        // processing_step = 1024 => 90278 ms
        // processing_step = 1024 * 64 => 40753 ms
        // processing_step = 1024 * 1024 => 38858 ms

        let expected_output = MNT4753Fr::new(BigInteger768([8181981188982771303, 9834648934716236448, 6420360685258842467, 14258691490360951478, 10642011566662929522, 16918207755479993617, 3581400602871836321, 14012664850056020974, 16755211538924649257, 4039951447678776727, 12365175056998155257, 119677729692145]));
        let num_leaves = 1024 * 1024;
        let mut tree = MNT4PoseidonMHT::init(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for _ in 0..num_leaves {
            tree.append(MNT4753Fr::rand(&mut rng));
        }
        tree.finalize_in_place();
        assert_eq!(tree.root().unwrap(), expected_output, "Output of the Merkle tree computation for MNT4 does not match to the expected value.");
    }

    #[test]
    fn merkle_tree_test_mnt6() {
        // running time for 1048576 leaves
        // processing_step = 1024 => 87798 ms
        // processing_step = 1024 * 64 => 39199 ms
        // processing_step = 1024 * 1024 => 37450 ms

        let expected_output = MNT6753Fr::new(BigInteger768([18065863015580309240, 1059485854425188866, 1479096878827665107, 6899132209183155323, 1829690180552438097, 7395327616910893705, 16132683753083562833, 8528890579558218842, 9345795575555751752, 8161305655297462527, 6222078223269068637, 401142754883827]));
        let num_leaves = 1024 * 1024;
        let mut tree = MNT6PoseidonMHT::init(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for _ in 0..num_leaves {
            tree.append(MNT6753Fr::rand(&mut rng));
        }
        tree.finalize_in_place();
        assert_eq!(tree.root().unwrap(), expected_output, "Output of the Merkle tree computation for MNT6 does not match to the expected value.");
    }

    #[test]
    fn merkle_tree_test_mnt4_empty_leaves() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        let max_leaves = 64;

        for num_leaves in 1..=max_leaves {

            // Generate random leaves
            let mut leaves = Vec::with_capacity(num_leaves);
            for _ in 0..num_leaves {
                leaves.push(MNT4753Fr::rand(&mut rng))
            }

            // Push them in a Naive Poseidon Merkle Tree and get the root
            leaves.extend_from_slice(vec![MNT4753Fr::zero(); max_leaves - num_leaves].as_slice());
            let naive_mt = NaiveMNT4PoseidonMHT::new(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = MNT4PoseidonMHT::init(max_leaves);
            leaves[0..num_leaves].iter().for_each(|&leaf| { mt.append(leaf); });
            let root = mt.finalize_in_place().root().unwrap();

            assert_eq!(naive_root, root);
        }

        // Test the case in which there are empty leaves interleaved with non empty ones
        // (i.e. misbehaviour of the user or non append only merkle tree)
        for num_leaves in 1..=max_leaves {
            //println!("Num leaves: {}", num_leaves);
            // Make half of the added leaves empty
            let mut leaves = Vec::with_capacity(num_leaves);
            for _ in 0..num_leaves/2 {
                leaves.push(MNT4753Fr::zero())
            }
            for _ in num_leaves/2..num_leaves {
                leaves.push(MNT4753Fr::rand(&mut rng))
            }

            // Push them in a Naive Poseidon Merkle Tree and get the root
            leaves.extend_from_slice(vec![MNT4753Fr::zero(); max_leaves - num_leaves].as_slice());
            let naive_mt = NaiveMNT4PoseidonMHT::new(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = MNT4PoseidonMHT::init(max_leaves);
            leaves[..].iter().for_each(|&leaf| { mt.append(leaf); });
            let root = mt.finalize_in_place().root().unwrap();

            assert_eq!(naive_root, root);
        }
    }

    #[test]
    fn merkle_tree_test_mnt6_empty_leaves() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        let max_leaves = 64;

        for num_leaves in 1..=max_leaves {

            // Generate random leaves
            let mut leaves = Vec::with_capacity(num_leaves);
            for _ in 0..num_leaves {
                leaves.push(MNT6753Fr::rand(&mut rng))
            }

            // Push them in a Naive Poseidon Merkle Tree and get the root
            leaves.extend_from_slice(vec![MNT6753Fr::zero(); max_leaves - num_leaves].as_slice());
            let naive_mt = NaiveMNT6PoseidonMHT::new(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = MNT6PoseidonMHT::init(max_leaves);
            leaves[..].iter().for_each(|&leaf| { mt.append(leaf); });
            let root = mt.finalize_in_place().root().unwrap();

            assert_eq!(naive_root, root);
        }

        // Test the case in which there are empty leaves interleaved with non empty ones
        // (i.e. misbehaviour of the user or non append only merkle tree)
        for num_leaves in 1..=max_leaves {

            // Make half of the added leaves empty
            let mut leaves = Vec::with_capacity(num_leaves);
            for _ in 0..num_leaves/2 {
                leaves.push(MNT6753Fr::zero())
            }
            for _ in num_leaves/2..num_leaves {
                leaves.push(MNT6753Fr::rand(&mut rng))
            }

            // Push them in a Naive Poseidon Merkle Tree and get the root
            leaves.extend_from_slice(vec![MNT6753Fr::zero(); max_leaves - num_leaves].as_slice());
            let naive_mt = NaiveMNT6PoseidonMHT::new(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = MNT6PoseidonMHT::init(max_leaves);
            leaves[0..num_leaves].iter().for_each(|&leaf| { mt.append(leaf); });
            let root = mt.finalize_in_place().root().unwrap();

            assert_eq!(naive_root, root);
        }
    }

    #[test]
    fn merkle_tree_path_test_mnt4() {

        let num_leaves = 64;
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut tree = MNT4PoseidonMHT::init(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for _ in 0..num_leaves/2 {
            let leaf = MNT4753Fr::rand(&mut rng);
            tree.append(leaf);
            leaves.push(leaf);
        }
        for _ in num_leaves/2..num_leaves {
            let leaf = MNT4753Fr::zero();
            tree.append(leaf);
            leaves.push(leaf);
        }
        tree.finalize_in_place();

        let naive_tree = NaiveMNT4PoseidonMHT::new(leaves.as_slice()).unwrap();

        let root = tree.root().unwrap();
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {
            let path = tree.get_merkle_path(i).unwrap();
            assert!(tree.verify_merkle_path(&leaves[i], &path).unwrap());

            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(&naive_root, &leaves[i]).unwrap());

            assert_eq!(path.0, naive_path.path);
        }
    }

    #[test]
    fn merkle_tree_path_test_mnt6() {

        let num_leaves = 64;
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut tree = MNT6PoseidonMHT::init(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for _ in 0..num_leaves/2 {
            let leaf = MNT6753Fr::rand(&mut rng);
            tree.append(leaf);
            leaves.push(leaf);
        }
        for _ in num_leaves/2..num_leaves {
            let leaf = MNT6753Fr::zero();
            tree.append(leaf);
            leaves.push(leaf);
        }
        tree.finalize_in_place();

        let naive_tree = NaiveMNT6PoseidonMHT::new(leaves.as_slice()).unwrap();

        let root = tree.root().unwrap();
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {
            let path = tree.get_merkle_path(i).unwrap();
            assert!(tree.verify_merkle_path(&leaves[i], &path).unwrap());

            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(&naive_root, &leaves[i]).unwrap());

            assert_eq!(path.0, naive_path.path);
        }
    }
}

use crate::merkle_tree::field_based_mht::FieldBasedMerkleTreeParameters;

use algebra::{
    field_new,
    biginteger::BigInteger768,
    fields::{
        mnt6753::Fr as MNT6753Fr,
        mnt4753::Fr as MNT4753Fr
    },
};

#[derive(Clone)]
pub struct MNT4753MHTPoseidonParameters;

impl FieldBasedMerkleTreeParameters for MNT4753MHTPoseidonParameters {
    type Data = MNT4753Fr;

    const MERKLE_ARITY: usize = 2;

    const EMPTY_HASH_CST: &'static [Self::Data] = &[
        field_new!(MNT4753Fr,BigInteger768([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
        field_new!(MNT4753Fr,BigInteger768([16176523673378133334, 16956778265158276239, 4733072280094522678, 4192433122035539299, 15174036148465069565, 1064698154771993694, 11910330977019062149, 6319782654186851533, 16700914731313914268, 5445268834830289758, 1400505001002362360, 471965325761536])),
        field_new!(MNT4753Fr,BigInteger768([4277889709792061947, 16378415495394357880, 7044435121910608910, 14821584290848291872, 7395794966562039924, 16993255131573590866, 5619183296459254929, 13099791189349530546, 10480178751052362035, 2637008226322547176, 3781718695868806314, 12315153018770])),
        field_new!(MNT4753Fr,BigInteger768([12602160878026473095, 6703218031846246807, 3261593422652552071, 5744898705605418863, 568988017168646903, 11287816066753693968, 10074855438408190401, 15629691871152873479, 15780837983413394866, 9677498205745683353, 2710763567607654741, 62078701997150])),
        field_new!(MNT4753Fr,BigInteger768([3953781730258530358, 6110249182111184881, 5451270732317761177, 6943311044605255929, 3250905642119052198, 499062261769816654, 10259239560863736633, 14722428766718467461, 1595932533524212605, 12636118867227484435, 2535281826946782271, 28256121815614])),
        field_new!(MNT4753Fr,BigInteger768([12881414061557326842, 2307157185518412145, 16729611732227341092, 13661588739810946674, 18234379091090570831, 13728387961021650531, 1569717769337480380, 3151613579174483992, 1013749114920705196, 18081081464301638855, 17662018612700913184, 270340369174749])),
        field_new!(MNT4753Fr,BigInteger768([3605040078201732329, 15148117684514223001, 2205284832938577511, 14243812222834942509, 4721805097502261162, 12477409373644239295, 2133828053265029846, 8064415567340640846, 9139665707130095289, 2327875196086669640, 7490804028311302854, 39972867252844])),
        field_new!(MNT4753Fr,BigInteger768([13474454514586827491, 10425256355538130430, 16532298926701769119, 6952341694048969702, 4310318211674876121, 17985109358262381251, 7018608422116534547, 856082501862907948, 10201142017622815327, 4014460796003869169, 18010264820431769424, 292492626559474])),
        field_new!(MNT4753Fr,BigInteger768([13139671956879947102, 11795592347454983858, 7434654901760918990, 15991845460188159262, 8926398827090641285, 11990848706687304309, 926869728656438786, 4306382339328932384, 7269461104939176981, 1807966305486065051, 4264822993380428385, 293851919357767])),
        field_new!(MNT4753Fr,BigInteger768([7342491680395881561, 4489070398956337739, 1841666804269854648, 6979644289443721867, 6071690188598067424, 3997963076382743365, 3039024652836142348, 14814384789583809377, 1655550027026311714, 9713495234097377698, 2137662669331280893, 88937546849400])),
        field_new!(MNT4753Fr,BigInteger768([16373180430190352304, 4568986638504379799, 11510087407464449489, 16156967732575703950, 6050766378378947199, 859918639245209489, 6772545515530587766, 13110119935032867725, 7036298432403266305, 17344570648793268650, 821472278359650501, 413286128211952])),
        field_new!(MNT4753Fr,BigInteger768([10064204450670869041, 6485002129695176050, 13703715585507433387, 6358528581148388798, 8317188422852594958, 16249587742229080712, 15320018320869903910, 6519053266577151375, 13012784699871617737, 5320294235171818985, 454727173294227562, 5938002571396])),
        field_new!(MNT4753Fr,BigInteger768([5970737698044827110, 12553390334279345905, 17425479610490698473, 10934643973186563254, 15507168698998809488, 12161816561081931432, 14261361224889855408, 16994719105388441699, 17735618555000182797, 9609548097148711049, 17046188423488746023, 179196001269934])),
        field_new!(MNT4753Fr,BigInteger768([16709709138701649262, 12385797615346539441, 15793702613532409439, 11792904936680731911, 12448525568586080213, 18093815366893683166, 17458183917083722104, 15386763209861070357, 15840640717629650353, 1633965272590164610, 6406553317189723786, 453724791849806])),
        field_new!(MNT4753Fr,BigInteger768([2251450467303238906, 18229227672555802069, 2321869528003277491, 5043025073640215091, 8897028974870798639, 1992482375368886049, 1127757339681285238, 4970155236198853509, 5531972124480915021, 11062705402253819529, 15409860081207635910, 63367562695789])),
        field_new!(MNT4753Fr,BigInteger768([12665956519007526801, 1843696063471531093, 4421617427066245380, 16793646830533228386, 9608384167236543090, 16246507232290209049, 8388997809959181075, 10840864805493139596, 8450266449020314034, 18230781272203046274, 14551001191227447302, 39411984747300])),
        field_new!(MNT4753Fr,BigInteger768([17645993219199470363, 4456989415638893456, 4408935425422831546, 11479469842080797658, 12461503208035377620, 8803849468867694413, 7077438435087789253, 18001987624227409493, 14865983638481701230, 12976208785697116336, 6171183493110737796, 275982657298262])),
        field_new!(MNT4753Fr,BigInteger768([12899372984440970302, 4799739196163719468, 18281200174184323078, 14946209440607152291, 18038633657381931669, 13067278118534518262, 5933725396271864440, 17619708040456630714, 2298350945899000737, 8542906148384703647, 4228741343160794490, 179821403404139])),
        field_new!(MNT4753Fr,BigInteger768([4909989095880468007, 17111380839816086128, 13050894464535600155, 17809974587333109512, 3086375548608236017, 4319641169843690011, 8065063935915414369, 473711108534546470, 16419368658227689488, 12529592002690167987, 13604515417941051506, 474131523927671])),
        field_new!(MNT4753Fr,BigInteger768([2169765295327214927, 12853779434032519714, 4188854069571610948, 11176715314675409911, 6762724088888884422, 17971049055147305269, 4462649898332426762, 6320507659953355380, 10630149020511712781, 7236213148398521505, 12209435975744724638, 169609105942177])),
        field_new!(MNT4753Fr,BigInteger768([4039022148289766684, 18287223727808800157, 5028272483573475465, 15052031213537407210, 4871429795110173763, 13614938552192121828, 12468252496702382358, 357389198149731327, 9783015964887059718, 7249799829448192, 15969956552103665981, 81227891365913])),
        field_new!(MNT4753Fr,BigInteger768([8605248520836760614, 7698601903285982577, 7884851945210851008, 11137343029012187974, 14941286058184996238, 8942300516617698513, 12939320316637942536, 1083599900484092912, 6146753984613830917, 8826262634545573139, 8688587480911778755, 260284710393087])),
        field_new!(MNT4753Fr,BigInteger768([7819867473996257895, 1163119882964275677, 11350473626422120371, 8626456057335088473, 14345605516748953137, 16749457994730221606, 16544765006115835698, 8298164333055185490, 4406983934583628718, 14285702076873139552, 5994563637574463256, 116897006820475])),
        field_new!(MNT4753Fr,BigInteger768([10594331094316413862, 70589505093706233, 5947914170501868337, 9672947291446030947, 10471508875756623171, 5885944209881926849, 12238184998651009858, 2356768587459654214, 1976222780735616667, 7837959657680688622, 15005043837627313635, 325697264245299])),
        field_new!(MNT4753Fr,BigInteger768([11667320706060706659, 4497687445765901980, 2058017719541602176, 12722409913033257750, 10570174209277577513, 8050204849865648454, 961757207376139324, 7111416262996090474, 7042025878446682219, 9563801595559695088, 7959547555032992673, 299817378928680])),
        field_new!(MNT4753Fr,BigInteger768([17875883576651124852, 12019580853957703123, 14052735765400810986, 12134502857846062868, 3932260060058116996, 7412890163800836041, 7364941959471116085, 13719981274545202396, 3896197460809059188, 10937048028897553703, 14723811072279331346, 153865728544816])),
        field_new!(MNT4753Fr,BigInteger768([6377003792482571837, 9388472707381902393, 16257705780321388021, 3660250735951570492, 6398042051117445808, 1214492883228958801, 6488083732806553102, 16343892194218810768, 10336715715449488978, 12278416360730941598, 15120896212829175083, 394513338866242])),
        field_new!(MNT4753Fr,BigInteger768([9245962539596129756, 11308784860300487285, 9743504241464063619, 2074639663493309096, 15871180161301662049, 6557679360035650783, 15143489300651627352, 9812097890750485822, 1244537599719308281, 10545958110682494900, 17248790119900992474, 249218627631122])),
        field_new!(MNT4753Fr,BigInteger768([11835214222516122436, 15764240437657974407, 6420372599378240894, 16652019220914568700, 7589938208001402659, 13631631940696103757, 11230090654519768203, 13440030055936664621, 14280625899794233793, 15100178106414451409, 16623067409628724177, 212213341437214])),
        field_new!(MNT4753Fr,BigInteger768([6011277968708967059, 639496817632482309, 5929013035997986788, 16843575660916066646, 1149777341859285338, 14658042833551766188, 7634621850133386722, 7001742147011157176, 7359717955614460602, 6888539823450308368, 17116969251913704728, 204877760505258])),
        field_new!(MNT4753Fr,BigInteger768([675673922872315140, 512129263392578794, 9798690560912251207, 7054764921560592919, 686766133538325336, 5178995292196027033, 6583362491349948775, 16058868603513050375, 12947923073962473868, 2752063852071387331, 12806379407361838115, 15307678462607])),
        field_new!(MNT4753Fr,BigInteger768([11094226676306026356, 892488786370986741, 5437666239468371688, 16672274046028921911, 9239873541740716904, 12705472627429935031, 7206833212006685522, 8231183089056512015, 10801086054025107055, 1571241152085700920, 9342851935948733002, 359062727739383])),
        field_new!(MNT4753Fr,BigInteger768([12626192863498559769, 4968852371277424393, 2530250434805537910, 16162102870268636043, 547272346507800090, 12183453615274424913, 8499797127502674163, 8056429342665427148, 7496571365272431209, 181055189816194176, 14815307604617191394, 102402652534269])),
        field_new!(MNT4753Fr,BigInteger768([13377265285437446591, 3352135748620234639, 13499357599871397619, 1538567981079493584, 16835432028939005585, 1489545390044353676, 5802843886319864877, 3744352011146277769, 10138325315351804883, 9706994107202364805, 6080057850105416814, 21683659967889])),
        field_new!(MNT4753Fr,BigInteger768([14676049290848132453, 2439091224842193539, 11542050736983134836, 4100153475221219504, 2196785088215953222, 17010857509769574328, 14064058056664542171, 11748564706229565562, 7464979168051417429, 6595033253934784911, 18169911910223096937, 215103939901774])),
        field_new!(MNT4753Fr,BigInteger768([17635318468772890277, 6352239780588639586, 10723339619007833446, 4893451802181675528, 4644556819212868158, 10153277421006788286, 9681039424445901541, 4852861664528550119, 8340680302631027082, 1123576236662129079, 16034785489282305855, 186509610595629])),
        field_new!(MNT4753Fr,BigInteger768([12625048523149826889, 5522016944643377594, 7699716105964784333, 18011247809484228968, 4314570080726502960, 4534820968063172891, 15750612471035961368, 16546982327248344102, 13106434188126237333, 2436901965804455265, 5509854058561642711, 235667917842447])),
        field_new!(MNT4753Fr,BigInteger768([15473980367787724738, 14437267393882795992, 8110764953514076111, 433491927694759431, 3501549144853175578, 4001144872122038393, 6285601585004779307, 7480052501517952400, 17706668432898352939, 17525543761380897562, 18262056004935591255, 100296896212791])),
        field_new!(MNT4753Fr,BigInteger768([13004537459536567505, 11536021368999062598, 11476994106473535335, 7501290503318089205, 1503497677790037124, 8571448125114677630, 2093681935226603521, 1783977282134180201, 815416671438766434, 14033390703273309148, 13863148357185408708, 151116486379211])),
        field_new!(MNT4753Fr,BigInteger768([11520390923272391174, 3700965319616748151, 16924840217636033521, 18866038153491911, 2401619975857281608, 15177168781119161934, 12160998724327471520, 9013888042188307789, 3884831131364280062, 4885763940431702386, 2558247536923573977, 230028098468810])),
        field_new!(MNT4753Fr,BigInteger768([10952132645556577918, 14768047366047066018, 16559210792506962762, 7202327027750626448, 569029245718616885, 6466881656623105602, 15350553108611606365, 4894926269159096183, 1368096661189637446, 14450950852108076518, 13525020817388864926, 256218987628845])),
        field_new!(MNT4753Fr,BigInteger768([11220546759813133736, 5574904932062982166, 1822072885686664942, 10305847449656783355, 786989818258593805, 8496531785376290371, 17372104123217289082, 11874517411562081834, 692404362408103552, 3619223599209338443, 7313814121803648884, 414498108514135]))
    ];
}

#[derive(Clone)]
pub struct MNT6753MHTPoseidonParameters;

impl FieldBasedMerkleTreeParameters for MNT6753MHTPoseidonParameters {
    type Data = MNT6753Fr;

    const MERKLE_ARITY: usize = 2;

    const EMPTY_HASH_CST: &'static [Self::Data] = &[
        field_new!(MNT6753Fr,BigInteger768([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
        field_new!(MNT6753Fr,BigInteger768([17637649811893425297, 9356568928171551347, 15442933386895351042, 16037249952112432370, 13788858697879839278, 11025667786216957645, 7102709371094193238, 12422915737218740758, 10531053957796416782, 9696092897380549051, 1339040383975805668, 64128018321799])),
        field_new!(MNT6753Fr,BigInteger768([4472923475352101931, 14125459602252897885, 13088438548806363052, 16269053428991819111, 3116242482118239773, 8191706716067520569, 4230313181893216955, 5167780582185616020, 15631071562248952054, 11846663030841039828, 642062297174172342, 296312004800026])),
        field_new!(MNT6753Fr,BigInteger768([17560533737151908951, 7063542403418931348, 9946122940491405418, 7701715352861997786, 9566898456569095206, 5256256878378084730, 9233983441268881659, 6942880360212085679, 6283174021297026980, 8102155274089156953, 2439467320428321035, 332916539951005])),
        field_new!(MNT6753Fr,BigInteger768([5858795400903007273, 7438824890337000111, 11488323460733706521, 2430825815000466084, 1725829885235881856, 15876812446098322793, 17156214869442521015, 5137679012965207550, 9116445458721095545, 4712054789035514492, 14030095712657958516, 179133400783271])),
        field_new!(MNT6753Fr,BigInteger768([9196483772299324859, 2927602375838435979, 16763111848481482763, 17882384337132970051, 17207183906514250255, 9530001320382362332, 3830621480341297746, 4866359965257025230, 16727937285691906509, 7542198491188832218, 8110695628458743283, 465520569394182])),
        field_new!(MNT6753Fr,BigInteger768([10747845174812323570, 2352057383464167932, 1141492571882677871, 1470351414409109323, 10420683829849635095, 15838896170163265003, 15955496784306175007, 4512671111438446024, 14027406595267595230, 4821799888258843127, 8552400658778235926, 90196730269457])),
        field_new!(MNT6753Fr,BigInteger768([5256253866873843525, 272593912443857551, 28470486217631823, 5159435739716868363, 6008486422443066719, 13667749945326593098, 12159573967831293532, 13247815893026277720, 5156136969387447701, 18025485311606231696, 3021473936786719040, 369755762512790])),
        field_new!(MNT6753Fr,BigInteger768([12344211293249999590, 12957596404284010059, 13910609108784582371, 17948699145681135586, 15785513802737509917, 9521697374075267890, 6365430956477394167, 10489625328383467046, 3109920408419137721, 4776466127615241425, 13808572562066969960, 309948976116781])),
        field_new!(MNT6753Fr,BigInteger768([6119070406732032896, 9410219222515712831, 5619723902868157253, 1973644971311309067, 6350780697601099563, 9024750504788049395, 15043082738210571628, 5620029319140107743, 3847387393513625150, 14928494257596235235, 2953722154564751554, 351659788771135])),
        field_new!(MNT6753Fr,BigInteger768([6049475571987158620, 14682333365358852589, 18135530363728134391, 14238712786862549175, 18420038755059004258, 7100569247008813954, 15681379881268007011, 3148544300132614406, 7346490772185170040, 2440588656045554503, 16485398900108052412, 60014821790753])),
        field_new!(MNT6753Fr,BigInteger768([12856661607330728400, 7016632322912905163, 3480013918990612535, 4344326050178292196, 16734656594695639398, 619569984615794816, 508464437174924696, 17238871601243639011, 9357837317046245838, 8883981285570502927, 7265239135733430929, 39139143918093])),
        field_new!(MNT6753Fr,BigInteger768([14817103890877107626, 6987260744665995883, 14570226874221084604, 11556329245783721885, 7466968604079736709, 9776987887571675901, 16458942593670591818, 9816159862029094791, 3921792482916243118, 6716123312977174659, 1710272864040696433, 421734633970436])),
        field_new!(MNT6753Fr,BigInteger768([4282397222691614301, 9879141210547300880, 12086373471014813377, 12666137842608178202, 12136963581046123545, 13264356597674291141, 13279969693572861251, 10611706899587547649, 5822212096659034083, 14127644322758565237, 14315577465274012830, 39177977839243])),
        field_new!(MNT6753Fr,BigInteger768([5728462243094819006, 17630581552464220954, 11320732721357334166, 1900213285723218149, 16241842790685525577, 15426354003830150652, 13204149184110252705, 12782089451147648433, 1732724535640987276, 10481021796229132770, 10419395107867006910, 344012985558625])),
        field_new!(MNT6753Fr,BigInteger768([344683888409415049, 17104079267154355794, 16575788368412272898, 8623293262069993159, 7982910632910187999, 7467718100746598937, 13092614656658983286, 16652493132760118102, 721157856720696161, 13335533611612451260, 781487980213691662, 442787986655123])),
        field_new!(MNT6753Fr,BigInteger768([6074433373638373965, 12098141975132028903, 1772727838370092714, 17302232439660752033, 13973402180971253892, 10747553266078344056, 4629001220352976062, 3362177028682656276, 4093103855723114167, 11999448184397189566, 4836796932935328252, 386630740659090])),
        field_new!(MNT6753Fr,BigInteger768([13014586413437199011, 14604625218698575416, 5071797836894621793, 17933033024926760524, 17079269583255311986, 6087229422530194396, 13294680830864986388, 401049880672557278, 12491328108961811657, 262021736727899607, 4727246708058880707, 110415138534356])),
        field_new!(MNT6753Fr,BigInteger768([1810339956058039811, 5220319720510094406, 13881603660491110179, 13222612110129007479, 6312365029261375458, 4488597311969070097, 3786930340526411451, 3100227749896650684, 17888997887837665338, 11120050721883473292, 17333024085167386788, 122622432905739])),
        field_new!(MNT6753Fr,BigInteger768([12015855797949312766, 15089252312128931938, 4487075872730733012, 11064976590409484842, 7345635811734382058, 6654267483537341812, 7338965532001215980, 8910608120598620104, 11760124634612512865, 16634036057775533219, 3154848179763801, 252988057458120])),
        field_new!(MNT6753Fr,BigInteger768([12844051781387933694, 14867417428795205800, 9706792882858210784, 2837497735373956323, 8708230781562375881, 16216709488172319510, 16685678840904390572, 17518208216783484556, 7106074686008369725, 18119542548943855133, 14132259050522202957, 74750391842824])),
        field_new!(MNT6753Fr,BigInteger768([2095049432688090524, 14022991336226067587, 13064251048943312063, 4930626082426515917, 7398876182738915571, 12675065241304118357, 11599240147379022093, 12104736345308339081, 7326994241467324793, 9543773052250172798, 2775763972445393468, 458168979271920])),
        field_new!(MNT6753Fr,BigInteger768([1946334007581914066, 17651035588256778549, 16514069863291089961, 6837961831333698346, 14290465216201654988, 11695022680475175792, 5272328103994695478, 6615418303775738534, 7306040626339560814, 1068008787604355877, 4115635962454045834, 495898874262831])),
        field_new!(MNT6753Fr,BigInteger768([14402861386561446620, 14792084567357226272, 6866700736595022709, 4113741901118049627, 16974812291891495925, 2591187165974885152, 9652009272067146832, 17945392524410072521, 16355201864774011574, 17722464991226331207, 6711417801116586799, 64077825502541])),
        field_new!(MNT6753Fr,BigInteger768([14064439024053015058, 8311654208515298479, 9104801255536293439, 3435835109271397379, 2712106172220266833, 16109486991502526028, 5651721887369042452, 14572741680048077241, 17167805857692983973, 10790465472671670937, 17136899827689865975, 235713138184448])),
        field_new!(MNT6753Fr,BigInteger768([10921122272069803471, 12577796758815504414, 17894410592699986884, 15650444086007959556, 948214525341676896, 4413473555476543334, 400409551525995658, 9569271116300630353, 14300418054165299370, 2405861297007405656, 8593207143586674229, 317679646497041])),
        field_new!(MNT6753Fr,BigInteger768([5010217285040896027, 6118190716370278919, 5568761447497243608, 7232401845348398675, 11018189049174018358, 9063262008998899553, 15414255187942422510, 15366085234988945047, 749508538410508494, 4179612678808007669, 9379657808494207396, 332004587403462])),
        field_new!(MNT6753Fr,BigInteger768([15596855953553551339, 14558516921801695117, 1716340212069481454, 16931771333182645310, 664561412384369638, 2277135588311724679, 17019393950247019203, 7528389714662451145, 4710861527742606713, 9268887905114908216, 4277691410811219790, 358017827571860])),
        field_new!(MNT6753Fr,BigInteger768([1345110017514752937, 11172153568636292065, 11079027760107165009, 766116359110021478, 16056140398387994894, 3290268012245266281, 4800923857475440103, 2198169219803152265, 7420936784792270326, 772819559725172974, 235282743392278103, 208876952517211])),
        field_new!(MNT6753Fr,BigInteger768([13471797275036477227, 13061727920320491911, 6819174308936409884, 8716139006117964600, 1962085276963290769, 6096628119158067780, 8167753635222191055, 14542801129845756191, 12572804873626916137, 9323731354266459335, 2575323808597388829, 11574582105067])),
        field_new!(MNT6753Fr,BigInteger768([4099487087924561648, 17735253112812409865, 1472196735946662562, 14589900298475437019, 13077017353725781269, 4656951940440280130, 2006510094889487317, 16199931703350800082, 14440178232943249115, 1407142698233581642, 13367748001204569270, 368868378632133])),
        field_new!(MNT6753Fr,BigInteger768([8074651267210831167, 8119908304743691932, 10050269642681585298, 17235916655576462100, 13335057030995613925, 10953418346798654943, 3695463758673253539, 5941835118077151115, 3992757082309531504, 9920313873598345988, 13151267180037628287, 207067907840698])),
        field_new!(MNT6753Fr,BigInteger768([894664452991370310, 16640116934783334734, 14808118850566711269, 6854546902707263589, 4286409562415851882, 6376630963274680056, 1706592727108815514, 2827938002676295743, 4157439222499882505, 8938576560015078325, 17999318348661348728, 34731319506718])),
        field_new!(MNT6753Fr,BigInteger768([11652860094199701250, 3478742230776828196, 10882615542005694703, 3349819079712688914, 8867992918381147826, 18409048493474927536, 11246279017112615234, 9748413453493631722, 17953656205798507701, 15187396101076825886, 11399324590810188311, 260131766860778])),
        field_new!(MNT6753Fr,BigInteger768([1249906864722423399, 12380849303036860652, 9152918792602847508, 8633501350856661146, 9693350457758570990, 14247290248527674449, 17646894536337687694, 4172434561479646685, 14342150478683294750, 7577792009144541865, 4695577900183865210, 317076072551130])),
        field_new!(MNT6753Fr,BigInteger768([8243914980137266903, 12333230796586460241, 14444364149875431527, 10803094216152777629, 4245055689606596787, 11044932568629364169, 11143770400531777907, 9358590092664905825, 4464390783738523435, 10767199047452580923, 15229918729430798984, 462864070917908])),
        field_new!(MNT6753Fr,BigInteger768([1332194145515023222, 5989346566927373592, 2377432556512621922, 7197981111191741919, 8933561838464171794, 14784476823475522561, 8478394278159923626, 5550504370540145528, 9532923115841735693, 14268124095180929026, 13757307827788564501, 122834449169298])),
        field_new!(MNT6753Fr,BigInteger768([14315080839281141763, 12316729172954982736, 12240902985496785994, 15572988659996267515, 2853659555022734346, 6458457554023043088, 18242519349120488704, 9074731373769712100, 16887310820873581593, 10961259943832775631, 7703684579928629068, 172779365404209])),
        field_new!(MNT6753Fr,BigInteger768([404644793350804679, 10645061623026087383, 7037216068474108354, 6875015450266344877, 9691319874633432046, 8981392801744671263, 7935787978429651629, 15636241802405345506, 1798934145016865467, 15862544822185043035, 1853057481023778743, 48582979057739])),
        field_new!(MNT6753Fr,BigInteger768([13841377248024908588, 6443916402067134263, 4121795598949110874, 17604038008411622297, 13539797093572147657, 8129794625753716120, 1778292717222571517, 1274414933390228330, 12302491198307697454, 12915166626799987755, 10282112543429553668, 337540712877314])),
        field_new!(MNT6753Fr,BigInteger768([2623819668573158668, 18068100053054478991, 11510571701283286325, 5846159179552586599, 7116703367486880097, 14386197053801562981, 15156005615023383543, 9383873929378927496, 14930130200844661115, 17345013830483094608, 11202472323720105939, 67169217781307])),
        field_new!(MNT6753Fr,BigInteger768([4403869937029740070, 12780380221027360032, 8411002685782179801, 9673734652350638919, 1271006590149990751, 14713224810005625570, 11720231750293065950, 5216848339135916669, 2768835827765483017, 9531440940757609710, 6191833244911324787, 187122239473111])),
    ];
}
