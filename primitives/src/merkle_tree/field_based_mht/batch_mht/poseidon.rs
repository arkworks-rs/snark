extern crate rand;

use crate::crh::BatchFieldBasedHash;
use algebra::{PrimeField, MulShort};
use crate::PoseidonParameters;
use crate::merkle_tree::field_based_mht::batch_mht::BatchMerkleTree;
use crate::crh::poseidon::batched_crh::PoseidonBatchHash;
use std::marker::PhantomData;
use std::clone::Clone;

#[derive(Clone)]
pub struct PoseidonBatchMerkleTree<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> {
    root: F,
    array_nodes: Vec<F>,
    processing_step: usize,
    initial_pos: Vec<usize>,
    final_pos: Vec<usize>,
    processed_pos: Vec<usize>,
    new_elem_pos: Vec<usize>,
    levels: usize,
    rate: usize,
    _parameters: PhantomData<P>
}

impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> PoseidonBatchMerkleTree<F, P> {
    pub fn root(&self) -> F {
        self.root.clone()
    }

    pub fn new(size_leaves: usize, processing_step: usize) -> Self {
        let rate = P::R;

        let mut initial_pos = Vec::new();
        let mut final_pos = Vec::new();
        let mut processed_pos = Vec::new();
        let mut new_elem_pos = Vec::new();

        let last_level_size = size_leaves.next_power_of_two();
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
            _parameters: PhantomData
        }
    }

    fn compute_subtree(&mut self) {
        for i in 0..self.levels  {

            if (self.new_elem_pos[i] - self.processed_pos[i]) >= self.rate {
                let num_groups_leaves = (self.new_elem_pos[i] - self.processed_pos[i]) / self.rate;
                let last_pos_to_process = self.processed_pos[i] + num_groups_leaves * self.rate;

                let (input_vec, output_vec) =
                    self.array_nodes[self.initial_pos[i]..self.final_pos[i + 1]].split_at_mut(self.final_pos[i] - self.initial_pos[i]);

                let new_pos_parent = self.new_elem_pos[i + 1] + num_groups_leaves;

                PoseidonBatchHash::<F,P>::batch_evaluate_in_place(&mut input_vec[(self.processed_pos[i] - self.initial_pos[i])..(last_pos_to_process - self.initial_pos[i])],
                                           &mut output_vec[(self.new_elem_pos[i + 1] - self.initial_pos[i + 1])..(new_pos_parent - self.initial_pos[i + 1])]);

                self.new_elem_pos[i + 1] += num_groups_leaves;
                self.processed_pos[i] += num_groups_leaves * self.rate;
            }
        }

    }
}

impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> BatchMerkleTree for PoseidonBatchMerkleTree<F,P> {

    type Hash = PoseidonBatchHash<F,P>;

    fn update(&mut self, leaf: F) {

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
    }

    fn finalize(&self) -> Self{
        let mut copy = (*self).clone();
        copy.new_elem_pos[0] = copy.final_pos[0];
        copy.compute_subtree();
        copy.root = *copy.array_nodes.last().unwrap();
        copy
    }

    fn finalize_in_place(&mut self) {
        self.new_elem_pos[0] = self.final_pos[0];
        self.compute_subtree();
        self.root = *self.array_nodes.last().unwrap();
    }

}


#[cfg(test)]
mod test {
    use rand_xorshift::XorShiftRng;
    use super::rand::SeedableRng;
    use crate::batched_crh::PoseidonBatchHash;
    use crate::crh::poseidon::parameters::MNT4753PoseidonParameters;
    use crate::crh::poseidon::parameters::MNT6753PoseidonParameters;
    use algebra::fields::mnt4753::Fr as MNT4753Fr;
    use algebra::fields::mnt6753::Fr as MNT6753Fr;
    use algebra::UniformRand;
    use std::time::Instant;

    use algebra::biginteger::BigInteger768;
    use crate::merkle_tree::field_based_mht::batch_mht::poseidon::PoseidonBatchMerkleTree;
    use crate::merkle_tree::field_based_mht::batch_mht::BatchMerkleTree;

    type MNT4BatchedMerkleTree = PoseidonBatchMerkleTree<MNT4753Fr, MNT4753PoseidonParameters>;
    type MNT6BatchedMerkleTree = PoseidonBatchMerkleTree<MNT6753Fr, MNT6753PoseidonParameters>;

    #[test]
    fn merkle_tree_test_mnt4() {
        // running time for 1048576 leaves
        // processing_step = 1024 => 90002 ms
        // processing_step = 1024 * 64 => 40612 ms
        // processing_step = 1024 * 1024 => 38617 ms

        let expected_output = MNT4753Fr::new(BigInteger768([8181981188982771303, 9834648934716236448, 6420360685258842467, 14258691490360951478, 10642011566662929522, 16918207755479993617, 3581400602871836321, 14012664850056020974, 16755211538924649257, 4039951447678776727, 12365175056998155257, 119677729692145]));

        let num_leaves = 1024 * 1024;
        let mut tree = MNT4BatchedMerkleTree::new(num_leaves, 1024 * 1024);

        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let now_4753 = Instant::now();
        for _ in 0..num_leaves {
            tree.update(MNT4753Fr::rand(&mut rng));
        }
        tree.finalize_in_place();
        let new_now_4753 = Instant::now();

        let duration_mnt4753 = new_now_4753.duration_since(now_4753);
        println!("Time for MT computation with {} leaves MNT4753 = {:?}", num_leaves, duration_mnt4753.as_millis());

        assert_eq!(tree.root(), expected_output, "Output of the Merkle tree computation for MNT4 does not match to the expected value.");
    }

    #[test]
    fn merkle_tree_test_mnt6() {
        // running time for 1048576 leaves
        // processing_step = 1024 => 91242 ms
        // processing_step = 1024 * 64 => 37168 ms
        // processing_step = 1024 * 1024 => 36838 ms

        let expected_output = MNT6753Fr::new(BigInteger768([18065863015580309240, 1059485854425188866, 1479096878827665107, 6899132209183155323, 1829690180552438097, 7395327616910893705, 16132683753083562833, 8528890579558218842, 9345795575555751752, 8161305655297462527, 6222078223269068637, 401142754883827]));

        let num_leaves = 1024 * 1024;
        let mut tree = MNT6BatchedMerkleTree::new(num_leaves, 1024 * 1024);

        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let now_6753 = Instant::now();
        for _ in 0..num_leaves {
            tree.update(MNT6753Fr::rand(&mut rng));
        }
        tree.finalize_in_place();
        let new_now_6753 = Instant::now();

        let duration_mnt6753 = new_now_6753.duration_since(now_6753);
        println!("Time for MT computation with {} leaves MNT6753 = {:?}", num_leaves, duration_mnt6753.as_millis());

        assert_eq!(tree.root(), expected_output, "Output of the Merkle tree computation for MNT6 does not match to the expected value.");
    }

}


#[derive(Clone)]
pub struct PoseidonBatchMerkleTreeMem<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>>{
    root: F,
    array_nodes: Vec<F>,
    processing_step: usize,
    final_pos: Vec<usize>,
    initial_pos_subarray: Vec<usize>,
    final_pos_subarray: Vec<usize>,
    processed_pos: Vec<usize>,
    processed_pos_subarray: Vec<usize>,
    new_elem_pos: Vec<usize>,
    new_elem_pos_subarray: Vec<usize>,
    levels: usize,
    rate: usize,
    _parameters: PhantomData<P>
}

impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> PoseidonBatchMerkleTreeMem<F, P> {
    pub fn root(&self) -> F {
        self.root.clone()
    }

    pub fn new(size_leaves: usize, processing_step: usize) -> Self {

        // change it to the rate of the hash function
        let rate = P::R;

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
        // the processing block is calculated to be a multiple of the number of cpus and the rate
        // if the processing step is not an exact multiple of the cpu and the rate then it is rounded down

        let last_level_size = size_leaves.next_power_of_two();
        // last_level_size is the minimum power of two that contains the leaves
        // if the number of leaves is smaller than a power of two then they will be filled by zero at the end.
        let mut size = last_level_size;

        assert!(last_level_size >= processing_block_size, "The number of leaves should be bigger than the processing chunk size.");

        if processing_block_size > last_level_size {
            // if the processing block is greater than the number of leaves, then
            // take this value as a maximum
            processing_block_size = last_level_size;
        }

        // final position of the full vector
        let mut final_pos = Vec::new();
        // initial position of the subarray
        let mut initial_pos_subarray = Vec::new();
        // final position of the subarray
        let mut final_pos_subarray = Vec::new();
        // postision until it has been processed relative to the subarray
        let mut processed_pos_subarray = Vec::new();
        // position until it has been processed relative to the full array
        let mut processed_pos = Vec::new();
        // position until new elements have been inserted in the subarray
        let mut new_elem_pos_subarray = Vec::new();
        // position until new elements have been inserted in the full array
        let mut new_elem_pos = Vec::new();


        // keeps track of the Merkle tree level
        let mut level_idx = 1;

        let mut initial_idx_subarray = 0;
        let mut initial_idx = 0;
        let mut final_idx = last_level_size;

        while size >= 1 {

            // the size of sub-block was calculate above
            let mut size_subblock = processing_block_size;
            // if the size the of the sub-array that is part of the tree is greater than
            // the corresponding real size, then use the real size to save memory space
            if size_subblock > size {
                size_subblock = size;
            }

            let final_idx_subarray = initial_idx_subarray + size_subblock;

            final_pos.push(final_idx);
            initial_pos_subarray.push(initial_idx_subarray);
            final_pos_subarray.push(final_idx_subarray);
            processed_pos.push(initial_idx);
            processed_pos_subarray.push(initial_idx_subarray);
            new_elem_pos.push(initial_idx);
            new_elem_pos_subarray.push(initial_idx_subarray);

            initial_idx += size;
            initial_idx_subarray += size_subblock;
            size /= rate; // the next level of the tree will have leaves / rate
            final_idx = initial_idx + size;
            level_idx += 1;
        }

        let tree_size = *final_pos_subarray.last().unwrap();

        let mut tree_nodes = Vec::with_capacity(tree_size);
        for _i in 0..tree_size {
            tree_nodes.push(F::zero());
        }

        Self {
            root: { F::zero() },
            array_nodes: { tree_nodes },
            processing_step: { processing_block_size },
            final_pos: { final_pos },
            initial_pos_subarray: { initial_pos_subarray },
            final_pos_subarray: { final_pos_subarray },
            processed_pos: { processed_pos },
            processed_pos_subarray: { processed_pos_subarray },
            new_elem_pos: { new_elem_pos },
            new_elem_pos_subarray: { new_elem_pos_subarray },
            levels: { level_idx - 2 },
            rate: { P::R },
            _parameters: PhantomData
        }
    }

    pub fn compute_subtree(&mut self) {
        for i in 0..self.levels {
            if (self.new_elem_pos_subarray[i] - self.processed_pos_subarray[i]) >= self.rate {
                let num_groups_leaves = (self.new_elem_pos_subarray[i] - self.processed_pos_subarray[i]) / self.rate;
                let last_pos_to_process = self.processed_pos_subarray[i] + num_groups_leaves * self.rate;


                let (input_vec, output_vec) =
                    self.array_nodes[self.initial_pos_subarray[i]..self.final_pos_subarray[i + 1]].split_at_mut(self.final_pos_subarray[i] - self.initial_pos_subarray[i]);

                let new_pos_parent = self.new_elem_pos_subarray[i + 1] + num_groups_leaves;

                PoseidonBatchHash::<F,P>::batch_evaluate_in_place(&mut input_vec[(self.processed_pos_subarray[i] - self.initial_pos_subarray[i])..(last_pos_to_process - self.initial_pos_subarray[i])],
                                           &mut output_vec[(self.new_elem_pos_subarray[i + 1] - self.initial_pos_subarray[i + 1])..(new_pos_parent - self.initial_pos_subarray[i + 1])]);

                self.new_elem_pos_subarray[i + 1] += num_groups_leaves;
                self.new_elem_pos[i + 1] += num_groups_leaves;
                self.processed_pos_subarray[i] += num_groups_leaves * self.rate;
                self.processed_pos[i] += num_groups_leaves * self.rate;

                if self.new_elem_pos_subarray[i] == self.final_pos_subarray[i] {
                    self.new_elem_pos_subarray[i] = self.initial_pos_subarray[i];
                }
                if self.processed_pos_subarray[i] == self.final_pos_subarray[i] {
                    self.processed_pos_subarray[i] = self.initial_pos_subarray[i];
                }
            }
        }
    }
}
impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> BatchMerkleTree for PoseidonBatchMerkleTreeMem<F,P> {

    type Hash = PoseidonBatchHash<F,P>;

    fn update(&mut self, leaf: F) {
        if self.new_elem_pos_subarray[0] < self.final_pos_subarray[0] {
            self.array_nodes[self.new_elem_pos_subarray[0]] = leaf;
            self.new_elem_pos[0] += 1;
            self.new_elem_pos_subarray[0] += 1;
        }

        if self.new_elem_pos_subarray[0] == self.final_pos_subarray[0] {
            self.compute_subtree();
        }

        if (self.new_elem_pos_subarray[0] - self.processed_pos_subarray[0]) >= self.processing_step {
            self.compute_subtree();
        }
    }

    fn finalize(&self) -> Self{
        let mut copy = (*self).clone();
        while copy.new_elem_pos[0] != copy.final_pos[0] {
            copy.update(F::zero());
        }
        copy.compute_subtree();
        copy.root = *copy.array_nodes.last().unwrap();
        copy
    }

    fn finalize_in_place(&mut self) {
        while self.new_elem_pos[0] != self.final_pos[0] {
            self.update(F::zero());
        }
        self.compute_subtree();
        self.root = (*self.array_nodes.last().unwrap()).clone();
    }
}


#[cfg(test)]
mod test_mem {
    use rand_xorshift::XorShiftRng;
    use super::rand::SeedableRng;
    use crate::batched_crh::PoseidonBatchHash;
    use crate::crh::poseidon::parameters::MNT4753PoseidonParameters;
    use crate::crh::poseidon::parameters::MNT6753PoseidonParameters;
    use algebra::fields::mnt4753::Fr as MNT4753Fr;
    use algebra::fields::mnt6753::Fr as MNT6753Fr;
    use algebra::UniformRand;

    use algebra::biginteger::BigInteger768;
    use crate::merkle_tree::field_based_mht::batch_mht::poseidon::PoseidonBatchMerkleTreeMem;
    use crate::merkle_tree::field_based_mht::batch_mht::BatchMerkleTree;
    use std::time::Instant;

    type MNT4BatchedMerkleTree = PoseidonBatchMerkleTreeMem<MNT4753Fr, MNT4753PoseidonParameters>;
    type MNT6BatchedMerkleTree = PoseidonBatchMerkleTreeMem<MNT6753Fr, MNT6753PoseidonParameters>;

    #[test]
    fn merkle_tree_test_mnt4() {
        // running time for 1048576 leaves
        // processing_step = 1024 => 90278 ms
        // processing_step = 1024 * 64 => 40753 ms
        // processing_step = 1024 * 1024 => 38858 ms

        let expected_output = MNT4753Fr::new(BigInteger768([8181981188982771303, 9834648934716236448, 6420360685258842467, 14258691490360951478, 10642011566662929522, 16918207755479993617, 3581400602871836321, 14012664850056020974, 16755211538924649257, 4039951447678776727, 12365175056998155257, 119677729692145]));

        let num_leaves = 1024 * 1024;
        let mut tree = MNT4BatchedMerkleTree::new(num_leaves, 1024 * 1024);

        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let now_4753 = Instant::now();
        for _ in 0..num_leaves {
            tree.update(MNT4753Fr::rand(&mut rng));
        }
        tree.finalize_in_place();
        let new_now_4753 = Instant::now();

        let duration_mnt4753 = new_now_4753.duration_since(now_4753);
        println!("Time for MT computation with {} leaves MNT4753 = {:?}", num_leaves, duration_mnt4753.as_millis());

        assert_eq!(tree.root(), expected_output, "Output of the Merkle tree computation for MNT4 does not match to the expected value.");
    }

    #[test]
    fn merkle_tree_test_mnt6() {
        // running time for 1048576 leaves
        // processing_step = 1024 => 87798 ms
        // processing_step = 1024 * 64 => 39199 ms
        // processing_step = 1024 * 1024 => 37450 ms

        let expected_output = MNT6753Fr::new(BigInteger768([18065863015580309240, 1059485854425188866, 1479096878827665107, 6899132209183155323, 1829690180552438097, 7395327616910893705, 16132683753083562833, 8528890579558218842, 9345795575555751752, 8161305655297462527, 6222078223269068637, 401142754883827]));

        let num_leaves = 1024 * 1024;
        let mut tree = MNT6BatchedMerkleTree::new(num_leaves, 1024 * 1024);

        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        let now_6753 = Instant::now();
        for _ in 0..num_leaves {
            tree.update(MNT6753Fr::rand(&mut rng));
        }
        tree.finalize_in_place();
        let new_now_6753 = Instant::now();

        let duration_mnt6753 = new_now_6753.duration_since(now_6753);
        println!("Time for MT computation with {} leaves MNT6753 = {:?}", num_leaves, duration_mnt6753.as_millis());

        assert_eq!(tree.root(), expected_output, "Output of the Merkle tree computation for MNT6 does not match to the expected value.");
    }

}
