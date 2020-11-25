use algebra::Field;
use crate::{BatchFieldBasedMerkleTreeParameters, BatchFieldBasedHash, FieldBasedMerkleTree, FieldBasedMerkleTreePath, FieldBasedMHTPath, FieldBasedHash, FieldBasedHashParameters, check_precomputed_parameters};
use std::marker::PhantomData;

/// An implementation of FieldBasedMerkleTree, optimized in time and memory,
/// and able to support any BatchFieldBasedHash and Merkle arity.
/// TODO: Test with arity > 2
#[derive(Clone)]
pub struct FieldBasedOptimizedMHT<T: BatchFieldBasedMerkleTreeParameters>{
    root: T::Data,
    // Stores all MT nodes
    array_nodes: Vec<T::Data>,
    processing_step: usize,
    // Stores the initial index of each level of the MT
    initial_pos: Vec<usize>,
    // Stores the final index of each level of the MT
    final_pos: Vec<usize>,
    // Stores the index up until the nodes were already hashed, for each level
    processed_pos: Vec<usize>,
    // Stores the next available index for each level of the MT
    new_elem_pos: Vec<usize>,
    rate: usize,
    height: usize,
    finalized: bool,
    _tree_parameters: PhantomData<T>,
}

impl<T: BatchFieldBasedMerkleTreeParameters> FieldBasedOptimizedMHT<T> {

    /// Creates a new tree given its `height` and `processing_step`, used to tune
    /// the memory usage of the tree. In particular, `processing_step` defines the
    /// number of leaves to store before triggering the computation of the hashes
    /// of the upper level.
    /// Decreasing `processing_step` leads to less memory consumption but
    /// significantly worsen performances, as the computation of upper levels' hashes
    /// is triggered more often, along with a batch hash that may perform poorly when
    /// there are few elements to batch; of course, conversely, increasing
    /// `processing_step` increases the memory usage too but improves performances.
    /// Meaningful values for `processing_step` are between 1 (i.e. update upper levels
    /// at each leaf), leading to best memory efficiency but worse performances, and
    /// the maximum number of leaves (or the mean number of leaves you plan to add),
    /// leading to worse memory efficiency but best performances (upper levels hashes
    /// are computed just once, but all the leaves must be kept in memory).
    /// You can edit the `mht_poseidon_tuning` benchmarks in `primitives/src/benches/poseidon_mht.rs`
    /// to properly tune the `processing_step` parameter according to your use case.
    pub fn init(height: usize, processing_step: usize) -> Self {
        //let height = height + 1;
        assert!(check_precomputed_parameters::<T>(height));

        let rate = <<T::H as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design.
        assert_eq!(rate, T::MERKLE_ARITY);

        let last_level_size = T::MERKLE_ARITY.pow(height as u32);
        assert!(processing_step > 0 && processing_step <= last_level_size);

        let mut initial_pos = Vec::new();
        let mut final_pos = Vec::new();
        let mut processed_pos = Vec::new();
        let mut new_elem_pos = Vec::new();

        let mut size = last_level_size;

        let mut initial_idx = 0;
        let mut final_idx = last_level_size;

        while size >= 1 {
            initial_pos.push(initial_idx);
            final_pos.push(final_idx);
            processed_pos.push(initial_idx);
            new_elem_pos.push(initial_idx);

            initial_idx += size;
            size /= rate;
            final_idx = initial_idx + size;
        }

        let tree_size = *final_pos.last().unwrap();

        let mut array_nodes = Vec::with_capacity(tree_size);
        for _i in 0..tree_size {
            array_nodes.push(<T::Data as Field>::zero());
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
            root: { <T::Data as Field>::zero() },
            array_nodes: { array_nodes },
            processing_step: { processing_block_size },
            initial_pos: { initial_pos },
            final_pos: { final_pos },
            processed_pos: { processed_pos },
            new_elem_pos: { new_elem_pos },
            rate,
            height,
            finalized: false,
            _tree_parameters: PhantomData,
        }
    }

    fn compute_subtree(&mut self) {
        for i in 0..=self.height  {
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
                );
                self.new_elem_pos[i + 1] += num_groups_leaves;
                self.processed_pos[i] += num_groups_leaves * self.rate;
            }
        }
    }

    fn batch_hash(input: &mut [T::Data], output: &mut [T::Data], parent_level: usize) {

        let mut i = 0;
        let empty = T::EMPTY_HASH_CST.unwrap().nodes[parent_level - 1];

        // Stores the chunk that must be hashed, i.e. the ones containing at least one non-empty node
        let mut to_hash = Vec::new();

        // Stores the positions in the output vec in which the hash of the above chunks must be
        // placed (assuming that it can happen to have non-consecutive all-empty chunks)
        let mut output_pos = Vec::new();

        // If the element of each chunk of "rate" size are all equals this means that the nodes
        // are all empty nodes (with overwhelming probability), therefore we already have the output,
        // otherwise it must be computed.
        input.chunks(T::MERKLE_ARITY).for_each(|input_chunk| {
            if input_chunk.iter().all(|&item| item == empty) {
                output[i] = T::EMPTY_HASH_CST.unwrap().nodes[parent_level];
            } else {
                to_hash.extend_from_slice(input_chunk);
                output_pos.push(i);
            }
            i += 1;
        });

        // Compute the hash of the non-all-empty chunks
        if to_hash.len() != 0 {
            let mut to_hash_out = vec![<T::Data as Field>::zero(); to_hash.len() / T::MERKLE_ARITY];
            <T::BH as BatchFieldBasedHash>::batch_evaluate_in_place(
                to_hash.as_mut_slice(),
                to_hash_out.as_mut_slice()
            );

            // Put the hashes in the correct positions in the output vec
            to_hash_out.iter().enumerate().for_each(|(i, &h)| output[output_pos[i]] = h);
        }
    }
}

impl<T: BatchFieldBasedMerkleTreeParameters> FieldBasedMerkleTree for FieldBasedOptimizedMHT<T> {
    type Parameters = T;
    type MerklePath = FieldBasedMHTPath<T>;

    /// Note: `Field` implements the `Copy` trait, therefore invoking this function won't
    /// cause a moving of ownership for `leaf`, but just a copy. Another copy is
    /// performed below in `self.array_nodes[self.new_elem_pos_subarray[0]] = leaf;`
    /// We can reduce this to one copy by passing a reference to leaf, but from an
    /// interface point of view this is not logically correct: someone calling this
    /// functions will likely not use the `leaf` anymore in most of the cases
    /// (in the other cases he can just clone it).
    fn append(&mut self, leaf: T::Data) -> &mut Self {
        if self.new_elem_pos[0] < self.final_pos[0] {
            self.array_nodes[self.new_elem_pos[0]] = leaf;
            self.new_elem_pos[0] += 1;
        }

        if self.new_elem_pos[0] == self.final_pos[0] {
            self.compute_subtree();
        }

        if (self.new_elem_pos[0] - self.processed_pos[0]) >= self.processing_step {
            self.compute_subtree();
        }
        self
    }

    fn finalize(&self) -> Self {
        let mut copy = (*self).clone();
        copy.new_elem_pos[0] = copy.final_pos[0];
        copy.compute_subtree();
        copy.finalized = true;
        copy.root = *copy.array_nodes.last().unwrap();
        copy
    }

    fn finalize_in_place(&mut self) -> &mut Self {
        self.new_elem_pos[0] = self.final_pos[0];
        self.compute_subtree();
        self.finalized = true;
        self.root = *self.array_nodes.last().unwrap();
        self
    }

    fn reset(&mut self) -> &mut Self {

        for i in 0..self.new_elem_pos.len() {
            self.new_elem_pos[i] = self.initial_pos[i];
            self.processed_pos[i] = self.initial_pos[i];
        }
        self.finalized = false;

        self
    }

    fn root(&self) -> Option<T::Data> {
        match self.finalized {
            true => Some(self.root.clone()),
            false => None
        }
    }

    fn get_merkle_path(&self, leaf_index: usize) -> Option<Self::MerklePath> {
        match self.finalized {
            true => {
                let num_leaves = T::MERKLE_ARITY.pow(self.height as u32);
                let mut merkle_path = Vec::with_capacity(self.height);

                let mut node_index = leaf_index;
                for _ in 0..self.height {
                    let mut siblings = Vec::with_capacity(T::MERKLE_ARITY - 1);

                    // Based on the index of the node, we must compute the index of the left-most children
                    let start_position = node_index - ( node_index % T::MERKLE_ARITY );

                    // Then, the right most children index is simply given by adding the arity
                    let end_position = start_position + T::MERKLE_ARITY;

                    // We must save the siblings of the actual node
                    for i in start_position..end_position {
                        if i != node_index {
                            siblings.push(self.array_nodes[i])
                        }
                    }

                    // Remember to normalize the node_index with respect to the length of siblings
                    merkle_path.push((siblings, node_index % T::MERKLE_ARITY));

                    // Get parent index for next iteration
                    node_index = num_leaves + (node_index/T::MERKLE_ARITY);
                }

                // Sanity check: the last node_index must be the one of the root
                assert_eq!(self.array_nodes[node_index], self.root);
                Some(
                    FieldBasedMHTPath::<T>::new(merkle_path)
                )
            },
            false => None,
        }
    }

    fn height(&self) -> usize { self.height }
}

#[cfg(test)]
mod test {
    use algebra::{
        biginteger::BigInteger768,
        fields::{
            Field,
            mnt4753::Fr as MNT4753Fr, mnt6753::Fr as MNT6753Fr
        },
        UniformRand,
        ToBytes, to_bytes, FromBytes,
    };
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    use crate::{
        crh::parameters::{MNT4PoseidonHash, MNT4BatchPoseidonHash, MNT6PoseidonHash, MNT6BatchPoseidonHash},
        merkle_tree::field_based_mht::{
        FieldBasedMerkleTree, NaiveMerkleTree,
        FieldBasedMerkleTreePath, FieldBasedMerkleTreeParameters,
        BatchFieldBasedMerkleTreeParameters, FieldBasedOptimizedMHT,
        parameters::{
            MNT4753_MHT_POSEIDON_PARAMETERS, MNT6753_MHT_POSEIDON_PARAMETERS
        }
    }, FieldBasedMerkleTreePrecomputedEmptyConstants, FieldBasedMHTPath};

    // OptimizedMHT definitions for tests below
    #[derive(Clone, Debug)]
    struct MNT4753FieldBasedOptimizedMerkleTreeParams;

    impl FieldBasedMerkleTreeParameters for MNT4753FieldBasedOptimizedMerkleTreeParams {
        type Data = MNT4753Fr;
        type H = MNT4PoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> =
            Some(MNT4753_MHT_POSEIDON_PARAMETERS);
    }

    impl BatchFieldBasedMerkleTreeParameters for MNT4753FieldBasedOptimizedMerkleTreeParams {
        type BH = MNT4BatchPoseidonHash;
    }

    type NaiveMNT4PoseidonMHT = NaiveMerkleTree<MNT4753FieldBasedOptimizedMerkleTreeParams>;
    type MNT4PoseidonMHT = FieldBasedOptimizedMHT<MNT4753FieldBasedOptimizedMerkleTreeParams>;
    type MNT4MerklePath = FieldBasedMHTPath<MNT4753FieldBasedOptimizedMerkleTreeParams>;

    #[derive(Clone, Debug)]
    struct MNT6753FieldBasedOptimizedMerkleTreeParams;

    impl FieldBasedMerkleTreeParameters for MNT6753FieldBasedOptimizedMerkleTreeParams {
        type Data = MNT6753Fr;
        type H = MNT6PoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> =
            Some(MNT6753_MHT_POSEIDON_PARAMETERS);
    }

    impl BatchFieldBasedMerkleTreeParameters for MNT6753FieldBasedOptimizedMerkleTreeParams {
        type BH = MNT6BatchPoseidonHash;
    }

    type NaiveMNT6PoseidonMHT = NaiveMerkleTree<MNT6753FieldBasedOptimizedMerkleTreeParams>;
    type MNT6PoseidonMHT = FieldBasedOptimizedMHT<MNT6753FieldBasedOptimizedMerkleTreeParams>;
    type MNT6MerklePath = FieldBasedMHTPath<MNT6753FieldBasedOptimizedMerkleTreeParams>;

    #[test]
    fn merkle_tree_test_mnt4() {
        let expected_output = MNT4753Fr::new(BigInteger768([8181981188982771303, 9834648934716236448, 6420360685258842467, 14258691490360951478, 10642011566662929522, 16918207755479993617, 3581400602871836321, 14012664850056020974, 16755211538924649257, 4039951447678776727, 12365175056998155257, 119677729692145]));
        let height = 20;
        let num_leaves = 2usize.pow(height as u32);
        let mut tree = MNT4PoseidonMHT::init(height, num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for _ in 0..num_leaves {
            tree.append(MNT4753Fr::rand(&mut rng));
        }
        tree.finalize_in_place();
        assert_eq!(tree.root().unwrap(), expected_output, "Output of the Merkle tree computation for MNT4 does not match to the expected value.");
    }

    #[test]
    fn merkle_tree_test_mnt6() {
        let expected_output = MNT6753Fr::new(BigInteger768([18065863015580309240, 1059485854425188866, 1479096878827665107, 6899132209183155323, 1829690180552438097, 7395327616910893705, 16132683753083562833, 8528890579558218842, 9345795575555751752, 8161305655297462527, 6222078223269068637, 401142754883827]));
        let height = 20;
        let num_leaves = 2usize.pow(height as u32);
        let mut tree = MNT6PoseidonMHT::init(height, num_leaves);
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
        let max_height = 6;
        let max_leaves = 2usize.pow(max_height as u32);

        for num_leaves in 1..=max_leaves {
            // Generate random leaves
            let mut leaves = Vec::with_capacity(num_leaves);
            for _ in 0..num_leaves {
                leaves.push(MNT4753Fr::rand(&mut rng))
            }

            // Push them in a Naive Poseidon Merkle Tree and get the root
            leaves.extend_from_slice(vec![MNT4753Fr::zero(); max_leaves - num_leaves].as_slice());
            let mut naive_mt = NaiveMNT4PoseidonMHT::new(max_height);
            naive_mt.append(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = MNT4PoseidonMHT::init(max_height, num_leaves);
            leaves[0..num_leaves].iter().for_each(|&leaf| { mt.append(leaf); });
            let root = mt.finalize_in_place().root().unwrap();

            assert_eq!(naive_root, root);
        }

        // Test the case in which there are empty leaves interleaved with non empty ones
        // (i.e. misbehaviour of the user or non append only merkle tree)
        for num_leaves in 1..=max_leaves {
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
            let mut naive_mt = NaiveMNT4PoseidonMHT::new(max_height);
            naive_mt.append(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = MNT4PoseidonMHT::init(max_height, num_leaves);
            leaves[..].iter().for_each(|&leaf| { mt.append(leaf); });
            let root = mt.finalize_in_place().root().unwrap();

            assert_eq!(naive_root, root);
        }
    }

    #[test]
    fn merkle_tree_test_mnt6_empty_leaves() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        let max_height = 6;
        let max_leaves = 2usize.pow(max_height as u32);

        for num_leaves in 1..=max_leaves {

            // Generate random leaves
            let mut leaves = Vec::with_capacity(num_leaves);
            for _ in 0..num_leaves {
                leaves.push(MNT6753Fr::rand(&mut rng))
            }

            // Push them in a Naive Poseidon Merkle Tree and get the root
            leaves.extend_from_slice(vec![MNT6753Fr::zero(); max_leaves - num_leaves].as_slice());
            let mut naive_mt = NaiveMNT6PoseidonMHT::new(max_height);
            naive_mt.append(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = MNT6PoseidonMHT::init(max_height, num_leaves);
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
            let mut naive_mt = NaiveMNT6PoseidonMHT::new(max_height);
            naive_mt.append(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = MNT6PoseidonMHT::init(max_height, num_leaves);
            leaves[0..num_leaves].iter().for_each(|&leaf| { mt.append(leaf); });
            let root = mt.finalize_in_place().root().unwrap();

            assert_eq!(naive_root, root);
        }
    }

    #[test]
    fn merkle_tree_path_test_mnt4() {

        let height = 6;
        let num_leaves = 2usize.pow(height as u32);
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut tree = MNT4PoseidonMHT::init(height, num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // Generate random leaves, half of which empty
        for _ in 0..num_leaves/2 {
            let leaf = MNT4753Fr::rand(&mut rng);
            tree.append(leaf);
            leaves.push(leaf);
        }
        for _ in num_leaves/2..num_leaves {
            let leaf = MNT4753Fr::zero();
            //tree.append(leaf);
            leaves.push(leaf);
        }

        // Compute the root of the tree, and do the same for a NaiveMHT, used here as reference
        tree.finalize_in_place();
        let mut naive_tree = NaiveMNT4PoseidonMHT::new(height);
        naive_tree.append(leaves.as_slice()).unwrap();
        let root = tree.root().unwrap();
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {

            // Create and verify a FieldBasedMHTPath
            let path = tree.get_merkle_path(i).unwrap();
            assert!(path.verify(tree.height(), &leaves[i], &root).unwrap());

            // Create and verify a Naive path
            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(naive_tree.height(), &leaves[i], &naive_root ).unwrap());

            // Assert the two paths are equal
            assert_eq!(naive_path, path);

            // Check leaf index is the correct one
            assert_eq!(i, path.leaf_index());

            if i == 0 { // leftmost check
                assert!(path.is_leftmost());
            }
            else if i == (num_leaves / 2) - 1 { // non-empty rightmost check
                assert!(path.is_non_empty_rightmost());
            }
            else if i == num_leaves - 1 { //rightmost check
                assert!(path.is_rightmost());
            }
            else { // Other cases check
                assert!(!path.is_leftmost());
                assert!(!path.is_rightmost());

                if i < (num_leaves / 2) - 1 {
                    assert!(!path.is_non_empty_rightmost());
                }
            }

            // Assert the two paths are equal
            assert_eq!(naive_path, path);

            // Serialization/deserialization test
            let path_serialized = to_bytes!(path).unwrap();
            let path_deserialized = MNT4MerklePath::read(path_serialized.as_slice()).unwrap();
            assert_eq!(path, path_deserialized);
        }
    }

    #[test]
    fn merkle_tree_path_is_non_empty_rightmost_test_mnt4() {

        let height = 6;
        let num_leaves = 2usize.pow(height as u32);
        let mut tree = MNT4PoseidonMHT::init(height, num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // Generate random leaves
        for i in 0..num_leaves {
            let leaf = MNT4753Fr::rand(&mut rng);
            tree.append(leaf);

            let tree_copy = tree.finalize();
            let path = tree_copy.get_merkle_path(i).unwrap();
            assert!(path.is_non_empty_rightmost());
        }
    }

    #[test]
    fn merkle_tree_path_test_mnt6() {

        let height = 6;
        let num_leaves = 2usize.pow(height as u32);
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut tree = MNT6PoseidonMHT::init(height, num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // Generate random leaves, half of which empty
        for _ in 0..num_leaves/2 {
            let leaf = MNT6753Fr::rand(&mut rng);
            tree.append(leaf);
            leaves.push(leaf);
        }
        for _ in num_leaves/2..num_leaves {
            let leaf = MNT6753Fr::zero();
            //tree.append(leaf);
            leaves.push(leaf);
        }

        // Compute the root of the tree, and do the same for a NaiveMHT, used here as reference
        tree.finalize_in_place();
        let mut naive_tree = NaiveMNT6PoseidonMHT::new(height);
        naive_tree.append(leaves.as_slice()).unwrap();
        let root = tree.root().unwrap();
        let naive_root = naive_tree.root();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {

            // Create and verify a FieldBasedMHTPath
            let path = tree.get_merkle_path(i).unwrap();
            assert!(path.verify(tree.height(), &leaves[i], &root).unwrap());

            // Create and verify a Naive path
            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.verify(naive_tree.height(), &leaves[i], &naive_root ).unwrap());

            // Assert the two paths are equal
            assert_eq!(naive_path, path);

            // Check leaf index is the correct one
            assert_eq!(i, path.leaf_index());

            if i == 0 { // leftmost check
                assert!(path.is_leftmost());
            }
            else if i == (num_leaves / 2) - 1 { // non-empty rightmost check
                assert!(path.is_non_empty_rightmost());
            }
            else if i == num_leaves - 1 { //rightmost check
                assert!(path.is_rightmost());
            }
            else { // Other cases check
                assert!(!path.is_leftmost());
                assert!(!path.is_rightmost());

                if i < (num_leaves / 2) - 1 {
                    assert!(!path.is_non_empty_rightmost());
                }
            }

            // Serialization/deserialization test
            let path_serialized = to_bytes!(path).unwrap();
            let path_deserialized = MNT6MerklePath::read(path_serialized.as_slice()).unwrap();
            assert_eq!(path, path_deserialized);
        }
    }
}