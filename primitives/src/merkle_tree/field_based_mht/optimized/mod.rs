use algebra::Field;
use crate::{Error, BatchFieldBasedMerkleTreeParameters, BatchFieldBasedHash, FieldBasedMerkleTree, FieldBasedMerkleTreePath, FieldBasedMHTPath, FieldBasedHash, FieldBasedHashParameters, check_precomputed_parameters, MerkleTreeError};
use std::marker::PhantomData;

/// An implementation of FieldBasedMerkleTree, optimized in time and memory,
/// and able to support any BatchFieldBasedHash and Merkle arity.
/// WARNING. This Merkle Tree implementation:
/// 1) Stores all the nodes in memory, so please retain from using it if
///    the available amount of memory is limited compared to the number
///    of leaves to be stored;
/// 2) Leaves and nodes are hashed without using any kind of domain separation:
///    while this is ok for use cases where the Merkle Trees have always the
///    same height, it's not for all the others.
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

    /// Creates a new tree given its `height` and `processing_step`, that defines the
    /// number of leaves to store before triggering the computation of the hashes
    /// of the upper levels. Changing this parameter will affect the performances of
    /// the batch hash: you can edit the `mht_poseidon_tuning` benchmarks in
    /// `primitives/src/benches/poseidon_mht.rs` to properly tune the `processing_step`
    /// parameter according to your use case.
    pub fn init(height: usize, processing_step: usize) -> Result<Self, Error> {
        if !check_precomputed_parameters::<T>(height) {
            Err(Box::new(MerkleTreeError::Other(format!(
                "Unsupported height. Max supported height is: {}",
                T::ZERO_NODE_CST.unwrap().nodes.len()
            ).to_owned())))?
        }

        let rate = <<T::H as FieldBasedHash>::Parameters as FieldBasedHashParameters>::R;
        // Rate may also be smaller than the arity actually, but this assertion
        // is reasonable and simplify the design.
        assert_eq!(rate, T::MERKLE_ARITY);

        let last_level_size = T::MERKLE_ARITY.pow(height as u32);
        if processing_step == 0 || processing_step > last_level_size {
            Err(Box::new(MerkleTreeError::Other(format!(
                "Invalid processing step. Must be between 1 and {}",
                last_level_size
            ).to_owned())))?
        }

        let mut initial_pos = Vec::new();
        let mut final_pos = Vec::new();
        let mut processed_pos = Vec::new();
        let mut new_elem_pos = Vec::new();

        let mut size = last_level_size;

        let mut initial_idx = 0;
        let mut final_idx = last_level_size;

        // Compute indexes
        while size >= 1 {
            initial_pos.push(initial_idx);
            final_pos.push(final_idx);
            processed_pos.push(initial_idx);
            new_elem_pos.push(initial_idx);

            initial_idx += size;
            size /= T::MERKLE_ARITY;
            final_idx = initial_idx + size;
        }

        let tree_size = *final_pos.last().unwrap();

        // Initialize to zero all tree nodes
        let mut array_nodes = Vec::with_capacity(tree_size);
        for _i in 0..tree_size {
            array_nodes.push(<T::Data as Field>::zero());
        }

        // Decide optimal processing block size based on the number of available
        // threads and the chosen processing step
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

        Ok(Self {
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
        })
    }

    /// Recursively update the nodes in the tree affected by some changes in the value of their children,
    /// starting from the leaves (as long as the number of affected nodes at a certain level is bigger
    /// equal than the rate of the hash function, therefore this function doesn't necessarily update
    /// all the nodes up until the root).
    fn compute_subtree(&mut self) -> Result<(), Error> {
        if self.height != 0 {
            for i in 0..=self.height  {

                // Enter only if the number of nodes to process at this level is bigger than the rate
                if (self.new_elem_pos[i] - self.processed_pos[i]) >= self.rate {

                    // The number of chunks of rate nodes to be processed
                    let num_groups_nodes = (self.new_elem_pos[i] - self.processed_pos[i]) / self.rate;

                    // Take as input vec all the nodes in the current level and all their parents
                    // (i.e. all the nodes at the next level)
                    let (input_vec, output_vec) =
                        self.array_nodes[self.initial_pos[i]..self.final_pos[i + 1]]
                            .split_at_mut(self.final_pos[i] - self.initial_pos[i]);

                    // The position of the last node in this level that will be affected by the changes.
                    // It's recomputed in this way as num_groups_nodes may have a remainder if
                    // the number of nodes to process it's not multiple of the rate.
                    let last_pos_to_process = self.processed_pos[i] + num_groups_nodes * self.rate;

                    // The position of the last node at parent level that will be affected by the changes.
                    let new_pos_parent = self.new_elem_pos[i + 1] + num_groups_nodes;

                    // From the input and the output vectors, use last_pos_to_process and new_pos_parent
                    // to isolate the nodes in this level and at parent level that are affected
                    // by changes, leaving the other ones out of the computation.
                    Self::batch_hash(
                        &mut input_vec[(self.processed_pos[i] - self.initial_pos[i])..(last_pos_to_process - self.initial_pos[i])],
                        &mut output_vec[(self.new_elem_pos[i + 1] - self.initial_pos[i + 1])..(new_pos_parent - self.initial_pos[i + 1])],
                        i + 1,
                    )?;

                    // Update new_elem_pos and processed_pos (in a consistent way as we did with
                    // new_pos_parent and last_pos_to_process.
                    self.new_elem_pos[i + 1] += num_groups_nodes;
                    self.processed_pos[i] += num_groups_nodes * self.rate;
                }
            }
        }

        Ok(())
    }

    pub fn get_leaves(&self) -> &[T::Data] {
        &self.array_nodes[self.initial_pos[0]..self.new_elem_pos[0]]
    }

    fn batch_hash(input: &mut [T::Data], output: &mut [T::Data], parent_level: usize) -> Result<(), Error> {

        let mut i = 0;
        let empty = T::ZERO_NODE_CST.unwrap().nodes[parent_level - 1];

        // Stores the chunk that must be hashed, i.e. the ones containing at least one non-empty node
        let mut to_hash = Vec::new();

        // Stores the positions in the output vec in which the hash of the above chunks must be
        // placed (assuming that it can happen to have non-consecutive all-empty chunks)
        let mut output_pos = Vec::new();

        // If the element of each chunk of "rate" size are all equals to the empty node at that level,
        // therefore we already have the output, otherwise it must be computed.
        input.chunks(T::MERKLE_ARITY).for_each(|input_chunk| {
            if input_chunk.iter().all(|&item| item == empty) {
                output[i] = T::ZERO_NODE_CST.unwrap().nodes[parent_level];
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
            )?;

            // Put the hashes in the correct positions in the output vec
            to_hash_out.iter().enumerate().for_each(|(i, &h)| output[output_pos[i]] = h);
        }

        Ok(())
    }
}

impl<T: BatchFieldBasedMerkleTreeParameters> FieldBasedMerkleTree for FieldBasedOptimizedMHT<T> {
    type Parameters = T;
    type MerklePath = FieldBasedMHTPath<T>;

    fn append(&mut self, leaf: T::Data) -> Result<&mut Self, Error> {

        // We can't take more leaves
        if self.processed_pos[0] == self.final_pos[0] {
            Err(MerkleTreeError::TooManyLeaves(self.height))?
        }

        // There is still room for more leaves
        if self.new_elem_pos[0] < self.final_pos[0] {
            self.array_nodes[self.new_elem_pos[0]] = leaf;
            self.new_elem_pos[0] += 1;

            // If the tree only has one leaf, we finished
            if self.height == 0 {
                self.processed_pos[0] = self.final_pos[0];
                return Ok(self);
            }
        }

        // With the previous update we reached the maximum capacity for the leaves
        if self.new_elem_pos[0] == self.final_pos[0] {
            self.compute_subtree()?;
        }

        // There is still room, but we can start updating the tree, according to
        // how we set processing_step parameter
        if (self.new_elem_pos[0] - self.processed_pos[0]) >= self.processing_step {
            self.compute_subtree()?;
        }

        Ok(self)
    }

    fn finalize(&self) -> Result<Self, Error> {
        let mut copy = (*self).clone();
        copy.new_elem_pos[0] = copy.final_pos[0];
        copy.compute_subtree()?;
        copy.finalized = true;
        if copy.array_nodes.len() == 0 {
            Err(MerkleTreeError::Other("Merkle tree is empty".to_owned()))?
        }
        copy.root = *copy.array_nodes.last().unwrap();
        Ok(copy)
    }

    fn finalize_in_place(&mut self) -> Result<&mut Self, Error> {
        self.new_elem_pos[0] = self.final_pos[0];
        self.compute_subtree()?;
        self.finalized = true;
        if self.array_nodes.len() == 0 {
            Err(MerkleTreeError::Other("Merkle tree is empty".to_owned()))?
        }
        self.root = *self.array_nodes.last().unwrap();
        Ok(self)
    }

    fn reset(&mut self) -> &mut Self {
        // Reset indices
        for i in 0..self.new_elem_pos.len() {
            self.new_elem_pos[i] = self.initial_pos[i];
            self.processed_pos[i] = self.initial_pos[i];
        }

        // Reset all nodes values
        self.array_nodes.iter_mut().for_each(|leaf| *leaf = <T::Data as Field>::zero());

        // Reset finalized value
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
        let num_leaves = T::MERKLE_ARITY.pow(self.height as u32);
        if leaf_index >= num_leaves {
            eprintln!("Invalid leaf index {} for num leaves {}", leaf_index, num_leaves);
            return None;
        }
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
                debug_assert_eq!(self.array_nodes[node_index], self.root);
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
        ToBytes, to_bytes, FromBytes, SemanticallyValid
    };
    use rand::{SeedableRng, RngCore, thread_rng};
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
        }, FieldBasedMerkleTreePrecomputedZeroConstants, FieldBasedMHTPath};

    // OptimizedMHT definitions for tests below
    #[derive(Clone, Debug)]
    struct MNT4753FieldBasedOptimizedMerkleTreeParams;

    impl FieldBasedMerkleTreeParameters for MNT4753FieldBasedOptimizedMerkleTreeParams {
        type Data = MNT4753Fr;
        type H = MNT4PoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const ZERO_NODE_CST: Option<FieldBasedMerkleTreePrecomputedZeroConstants<'static, Self::H>> =
            Some(MNT4753_MHT_POSEIDON_PARAMETERS);
    }

    impl BatchFieldBasedMerkleTreeParameters for MNT4753FieldBasedOptimizedMerkleTreeParams {
        type BH = MNT4BatchPoseidonHash;
    }

    #[derive(Clone, Debug)]
    struct MNT6753FieldBasedOptimizedMerkleTreeParams;

    impl FieldBasedMerkleTreeParameters for MNT6753FieldBasedOptimizedMerkleTreeParams {
        type Data = MNT6753Fr;
        type H = MNT6PoseidonHash;
        const MERKLE_ARITY: usize = 2;
        const ZERO_NODE_CST: Option<FieldBasedMerkleTreePrecomputedZeroConstants<'static, Self::H>> =
            Some(MNT6753_MHT_POSEIDON_PARAMETERS);
    }

    impl BatchFieldBasedMerkleTreeParameters for MNT6753FieldBasedOptimizedMerkleTreeParams {
        type BH = MNT6BatchPoseidonHash;
    }

    fn merkle_tree_root_test<T: BatchFieldBasedMerkleTreeParameters, R: RngCore>(
        height:        usize,
        num_leaves:    usize,
        expected_root: T::Data,
        mut rng:       &mut R
    )
    {
        // Init in memory optimized tree
        let mut tree = FieldBasedOptimizedMHT::<T>::init(height, num_leaves).unwrap();

        // Init naive merkle tree used as comparison
        let mut naive_mt = NaiveMerkleTree::<T>::new(height);

        // Create leaves at random
        let leaves = (0..num_leaves).map(|_| T::Data::rand(&mut rng)).collect::<Vec<_>>();

        // Append leaves to tree
        leaves.iter().for_each(|leaf| {
            tree.append(leaf.clone()).unwrap();
        });

        // Append leaves to naive_mt
        naive_mt.append(leaves.as_slice()).unwrap();

        // Exceeding maximum leaves will result in an error
        assert!(tree.append(leaves[0].clone()).is_err());

        // Finalize tree and get roots of both
        tree.finalize_in_place().unwrap();

        let optimized_root = tree.root().unwrap();
        let naive_root = naive_mt.root().unwrap();
        assert_eq!(naive_root, optimized_root);
        assert_eq!(
            tree.root().unwrap(),
            expected_root,
        );
    }

    /// Tests that effectively all the nodes of the tree are zeroed after a reset
    fn merkle_tree_reset_test<T: BatchFieldBasedMerkleTreeParameters, R: RngCore>(
        height:        usize,
        num_leaves:    usize,
        mut rng:       &mut R
    )
    {
        // Init in memory optimized tree
        let mut tree = FieldBasedOptimizedMHT::<T>::init(height, num_leaves).unwrap();

        // Create leaves at random
        let leaves = (0..num_leaves).map(|_| T::Data::rand(&mut rng)).collect::<Vec<_>>();

        // Add leaves to tree (don't fill the tree completely)
        leaves[..num_leaves/2].iter().for_each(|leaf| {
            tree.append(leaf.clone()).unwrap();
        });

        // This is the root we should get after having reset the tree if all the nodes
        // have been zeroed.
        let expected_root = tree.finalize().unwrap().root().unwrap();

        // Finish filling the tree
        leaves[num_leaves/2..].iter().for_each(|leaf| {
            tree.append(leaf.clone()).unwrap();
        });

        // Reset the tree
        tree.finalize_in_place().unwrap().reset();

        // Add the same leaves as we did initially
        leaves[..num_leaves/2].iter().for_each(|leaf| {
            tree.append(leaf.clone()).unwrap();
        });

        // Now, if all the nodes have been zeroed, than the assertion below will pass;
        // otherwise, this means that nodes still retained their old values, so the
        // computed root will be different.
        assert_eq!(expected_root, tree.finalize().unwrap().root().unwrap());
    }

    fn merkle_tree_test_edge_cases<T: BatchFieldBasedMerkleTreeParameters>() {
        let rng = &mut thread_rng();

        // HEIGHT > 1
        {
            // Generate empty tree and attempt to finalize
            let tree = FieldBasedOptimizedMHT::<T>::init(5, 1).unwrap();
            tree.finalize().unwrap();

            // Generate tree with only 1 leaf and attempt to finalize
            let mut tree = FieldBasedOptimizedMHT::<T>::init(5, 1).unwrap();
            assert!(tree.append(T::Data::rand(rng)).is_ok());
            tree.finalize().unwrap();
        }

        // HEIGHT == 1
        {
            // Generate empty tree and attempt to finalize
            let tree = FieldBasedOptimizedMHT::<T>::init(1, 2).unwrap();
            tree.finalize().unwrap();

            // Generate tree with only 1 leaf and attempt to finalize
            let mut tree = FieldBasedOptimizedMHT::<T>::init(1, 2).unwrap();
            assert!(tree.append(T::Data::rand(rng)).is_ok());
            tree.finalize().unwrap();

            // Generate tree with exactly 2 leaves and attempt to finalize.
            // Assert error if trying to add another leaf
            let mut tree = FieldBasedOptimizedMHT::<T>::init(1, 2).unwrap();
            assert!(tree.append(T::Data::rand(rng)).is_ok());
            assert!(tree.append(T::Data::rand(rng)).is_ok());
            assert!(tree.append(T::Data::rand(rng)).is_err());
            tree.finalize().unwrap();
        }

        // HEIGHT == 0
        {
            // Generate empty tree and attempt to finalize
            let mut tree = FieldBasedOptimizedMHT::<T>::init(0, 1).unwrap();
            tree.finalize_in_place().unwrap();
            assert_eq!(tree.root().unwrap(), T::Data::zero());

            // Check that reset() works properly also in this case
            tree.reset();
            tree.finalize_in_place().unwrap();
            assert_eq!(tree.root().unwrap(), T::Data::zero());

            // Generate tree and add one leaf. Assert adding one more leaf causes an error.
            let mut tree = FieldBasedOptimizedMHT::<T>::init(0, 1).unwrap();
            let fe = T::Data::rand(rng);
            assert!(tree.append(fe).is_ok());
            assert!(tree.append(T::Data::rand(rng)).is_err());
            tree.finalize_in_place().unwrap();
            assert_eq!(fe, tree.root().unwrap());

            // Check that reset() works properly also in this case
            tree.reset();
            let fe = T::Data::rand(rng);
            assert!(tree.append(fe).is_ok());
            tree.finalize_in_place().unwrap();
            assert_eq!(fe, tree.root().unwrap());
        }
    }

    #[test]
    fn merkle_tree_test_mnt4() {
        let expected_output = MNT4753Fr::new(BigInteger768([
            11737642701305799951,
            16779001331075430197,
            11819169129328038354,
            11423404101688341353,
            13644857877536036127,
            136974075146428157,
            13736146501659167139,
            15457726208981564885,
            16287955982068396368,
            2574770790166887043,
            15847921958357229891,
            431926751316706
        ]));
        let height = 10;
        let num_leaves = 2usize.pow(height as u32);
        let rng = &mut XorShiftRng::seed_from_u64(1231275789u64);

        merkle_tree_root_test::<MNT4753FieldBasedOptimizedMerkleTreeParams, _>(height, num_leaves, expected_output, rng);
        merkle_tree_reset_test::<MNT4753FieldBasedOptimizedMerkleTreeParams, _>(height, num_leaves,rng);
        merkle_tree_test_edge_cases::<MNT4753FieldBasedOptimizedMerkleTreeParams>();
    }

    #[test]
    fn merkle_tree_test_mnt6() {
        let expected_output = MNT6753Fr::new(BigInteger768([
            8485425859071260580,
            10496086997731513209,
            4252500720562453591,
            2141019788822111914,
            14051983083211686650,
            1024951982785915663,
            15435931545111578451,
            10317608288193115884,
            14391757241795953360,
            10971839229749467698,
            17614506209597433225,
            374251447408225
        ]));
        let height = 10;
        let num_leaves = 2usize.pow(height as u32);

        let rng = &mut XorShiftRng::seed_from_u64(1231275789u64);

        merkle_tree_root_test::<MNT6753FieldBasedOptimizedMerkleTreeParams, _>(height, num_leaves, expected_output, rng);
        merkle_tree_reset_test::<MNT6753FieldBasedOptimizedMerkleTreeParams, _>(height, num_leaves,rng);
        merkle_tree_test_edge_cases::<MNT6753FieldBasedOptimizedMerkleTreeParams>();
    }

    fn merkle_tree_test_empty_leaves<T: BatchFieldBasedMerkleTreeParameters, R: RngCore>(
        max_height: usize,
        max_leaves: usize,
        mut rng: &mut R,
    ) {
        for num_leaves in 1..=max_leaves {
            // Generate random leaves
            let mut leaves = Vec::with_capacity(num_leaves);
            for _ in 0..num_leaves {
                leaves.push(T::Data::rand(&mut rng))
            }

            // Push them in a Naive Poseidon Merkle Tree and get the root
            leaves.extend_from_slice(vec![<T::Data as Field>::zero(); max_leaves - num_leaves].as_slice());
            let mut naive_mt = NaiveMerkleTree::<T>::new(max_height);
            naive_mt.append(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root().unwrap();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = FieldBasedOptimizedMHT::<T>::init(max_height, num_leaves).unwrap();
            leaves[0..num_leaves].iter().for_each(|&leaf| { mt.append(leaf).unwrap(); });
            let root = mt.finalize_in_place().unwrap().root().unwrap();

            assert_eq!(naive_root, root);
        }

        // Test the case in which there are empty leaves interleaved with non empty ones
        // (i.e. misbehaviour of the user or non append only merkle tree)
        for num_leaves in 1..=max_leaves {
            // Make half of the added leaves empty
            let mut leaves = Vec::with_capacity(num_leaves);
            for _ in 0..num_leaves/2 {
                leaves.push(<T::Data as Field>::zero())
            }
            for _ in num_leaves/2..num_leaves {
                leaves.push(T::Data::rand(&mut rng))
            }

            // Push them in a Naive Poseidon Merkle Tree and get the root
            leaves.extend_from_slice(vec![<T::Data as Field>::zero(); max_leaves - num_leaves].as_slice());
            let mut naive_mt = NaiveMerkleTree::<T>::new(max_height);
            naive_mt.append(leaves.as_slice()).unwrap();
            let naive_root = naive_mt.root().unwrap();

            // Push them in a Poseidon Merkle Tree and get the root
            let mut mt = FieldBasedOptimizedMHT::<T>::init(max_height, num_leaves).unwrap();
            leaves[..].iter().for_each(|&leaf| { mt.append(leaf).unwrap(); });
            let root = mt.finalize_in_place().unwrap().root().unwrap();

            assert_eq!(naive_root, root);
        }
    }

    #[test]
    fn merkle_tree_test_mnt4_empty_leaves() {
        let rng = &mut XorShiftRng::seed_from_u64(1231275789u64);
        let max_height = 6;
        let max_leaves = 2usize.pow(max_height as u32);

        merkle_tree_test_empty_leaves::<MNT4753FieldBasedOptimizedMerkleTreeParams, _>(max_height, max_leaves, rng)
    }

    #[test]
    fn merkle_tree_test_mnt6_empty_leaves() {
        let rng = &mut XorShiftRng::seed_from_u64(1231275789u64);
        let max_height = 6;
        let max_leaves = 2usize.pow(max_height as u32);

        merkle_tree_test_empty_leaves::<MNT6753FieldBasedOptimizedMerkleTreeParams, _>(max_height, max_leaves, rng)
    }

    fn merkle_tree_path_test<T: BatchFieldBasedMerkleTreeParameters, R: RngCore>(
        height: usize,
        num_leaves: usize,
        mut rng: &mut R,
    ) {
        let mut leaves = Vec::with_capacity(num_leaves);
        let mut tree = FieldBasedOptimizedMHT::<T>::init(height, num_leaves).unwrap();

        // Generate random leaves, half of which empty
        for _ in 0..num_leaves/2 {
            let leaf = T::Data::rand(&mut rng);
            tree.append(leaf).unwrap();
            leaves.push(leaf);
        }
        for _ in num_leaves/2..num_leaves {
            let leaf = <T::Data as Field>::zero();
            leaves.push(leaf);
        }

        // Compute the root of the tree, and do the same for a NaiveMHT, used here as reference
        tree.finalize_in_place().unwrap();
        let mut naive_tree = NaiveMerkleTree::<T>::new(height);
        naive_tree.append(leaves.as_slice()).unwrap();
        let root = tree.root().unwrap();
        let naive_root = naive_tree.root().unwrap();
        assert_eq!(root, naive_root);

        for i in 0..num_leaves {

            // Create and verify a FieldBasedMHTPath
            let path = tree.get_merkle_path(i).unwrap();
            assert!(path.is_valid());
            assert!(path.verify(tree.height(), &leaves[i], &root).unwrap());

            // Create and verify a Naive path
            let naive_path = naive_tree.generate_proof(i, &leaves[i]).unwrap();
            assert!(naive_path.is_valid());
            assert!(naive_path.verify(naive_tree.height(), &leaves[i], &naive_root ).unwrap());

            // Assert the two paths are equal
            assert_eq!(naive_path, path);

            // Check leaf index is the correct one
            assert_eq!(i, path.leaf_index());

            if i == 0 { // leftmost check
                assert!(path.is_leftmost());
            }
            else if i == (num_leaves / 2) - 1 { // non-empty rightmost check
                assert!(path.are_right_leaves_empty());
            }
            else if i == num_leaves - 1 { //rightmost check
                assert!(path.is_rightmost());
            }
            else { // Other cases check
                assert!(!path.is_leftmost());
                assert!(!path.is_rightmost());

                if i < (num_leaves / 2) - 1 {
                    assert!(!path.are_right_leaves_empty());
                }
            }

            // Assert the two paths are equal
            assert_eq!(naive_path, path);

            // Serialization/deserialization test
            let path_serialized = to_bytes!(path).unwrap();
            let path_deserialized = FieldBasedMHTPath::<T>::read(path_serialized.as_slice()).unwrap();
            assert_eq!(path, path_deserialized);
        }
    }

    fn merkle_tree_path_are_right_leaves_empty_test<T: BatchFieldBasedMerkleTreeParameters, R: RngCore>(
        height: usize,
        num_leaves: usize,
        mut rng: &mut R,
    )
    {
        let mut tree = FieldBasedOptimizedMHT::<T>::init(height, num_leaves).unwrap();

        // Generate random leaves
        for i in 0..num_leaves {
            let leaf = T::Data::rand(&mut rng);
            tree.append(leaf).unwrap();

            let tree_copy = tree.finalize().unwrap();
            let path = tree_copy.get_merkle_path(i).unwrap();
            assert!(path.are_right_leaves_empty());
        }
    }

    #[test]
    fn merkle_tree_path_test_mnt4() {

        let height = 6;
        let num_leaves = 2usize.pow(height as u32);
        let rng = &mut XorShiftRng::seed_from_u64(1231275789u64);

        merkle_tree_path_test::<MNT4753FieldBasedOptimizedMerkleTreeParams, _>(height, num_leaves, rng);
        merkle_tree_path_are_right_leaves_empty_test::<MNT4753FieldBasedOptimizedMerkleTreeParams, _>(height, num_leaves, rng);
    }

    #[test]
    fn merkle_tree_path_test_mnt6() {

        let height = 6;
        let num_leaves = 2usize.pow(height as u32);
        let rng = &mut XorShiftRng::seed_from_u64(1231275789u64);

        merkle_tree_path_test::<MNT6753FieldBasedOptimizedMerkleTreeParams, _>(height, num_leaves, rng);
        merkle_tree_path_are_right_leaves_empty_test::<MNT6753FieldBasedOptimizedMerkleTreeParams, _>(height, num_leaves, rng);
    }
}