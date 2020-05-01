extern crate rand;

use crate::crh::{
    BatchFieldBasedHash,
    FieldBasedHashParameters, poseidon::{
        parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters}
    }
};
use algebra::Field;

pub struct BatchMerkleTree<H: BatchFieldBasedHash>{
    root: H::Data,
    array_nodes: Vec<H::Data>,
    size_leaves: usize,
    processing_step: usize,
    initial_pos: Vec<usize>,
    final_pos: Vec<usize>,
    processed_pos: Vec<usize>,
    new_elem_pos: Vec<usize>,
    levels: usize,
    rate: usize
}

impl<H: BatchFieldBasedHash> BatchMerkleTree<H> {

    pub fn root(&self) -> H::Data {
        self.root.clone()
    }

    pub fn new(size_leaves: usize, processing_step: usize) -> Self {

        let rate = 2;

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
            array_nodes.push(H::Data::zero());
        }

        let cpus = rayon::current_num_threads();
        let chunk_size = processing_step / (cpus * rate);
        let mut processing_block_size = chunk_size * cpus * rate;
        if processing_block_size > last_level_size {
            processing_block_size = last_level_size;
        }

        Self {
            root: { H::Data::zero() },
            array_nodes: { array_nodes },
            size_leaves: { size_leaves },
            processing_step: { processing_block_size },
            initial_pos: { initial_pos },
            final_pos: { final_pos },
            processed_pos: { processed_pos },
            new_elem_pos: { new_elem_pos },
            levels: { level_idx - 2 },
            rate: { 2 }
        }
    }

    pub fn push(&mut self, elem: H::Data) {
        if self.new_elem_pos[0] < self.final_pos[0] {
            self.array_nodes[self.new_elem_pos[0]] = elem;
            self.new_elem_pos[0] += 1;
        }

        if self.new_elem_pos[0] == self.final_pos[0] {
            self.update();
        }

        if (self.new_elem_pos[0] - self.processed_pos[0]) >= self.processing_step {
            self.update();
        }
     }

    pub fn finalize(&mut self) {
        self.new_elem_pos[0] = self.final_pos[0];
        self.update();
        self.root = (*self.array_nodes.last().unwrap()).clone()
    }

    pub fn update(&mut self) {
        for i in 0..self.levels  {

            if (self.new_elem_pos[i] - self.processed_pos[i]) >= self.rate {
                let num_groups_leaves = (self.new_elem_pos[i] - self.processed_pos[i]) / self.rate;
                let last_pos_to_process = self.processed_pos[i] + num_groups_leaves * self.rate;

                let (input_vec, output_vec) =
                    self.array_nodes[self.initial_pos[i]..self.final_pos[i + 1]].split_at_mut(self.final_pos[i] - self.initial_pos[i]);

                let new_pos_parent = self.new_elem_pos[i + 1] + num_groups_leaves;

                H::batch_evaluate_in_place(&mut input_vec[(self.processed_pos[i] - self.initial_pos[i])..(last_pos_to_process - self.initial_pos[i])],
                                           &mut output_vec[(self.new_elem_pos[i + 1] - self.initial_pos[i + 1])..(new_pos_parent - self.initial_pos[i + 1])]);

                self.new_elem_pos[i + 1] += num_groups_leaves;
                self.processed_pos[i] += num_groups_leaves * self.rate;
            }
        }

    }

}

#[cfg(test)]
mod test {
    use rand_xorshift::XorShiftRng;
    use super::rand::SeedableRng;
    use crate::PoseidonBatchHash;
    use crate::crh::poseidon::parameters::MNT4753PoseidonParameters;
    use crate::crh::poseidon::parameters::MNT6753PoseidonParameters;
    use crate::merkle_tree::field_based_mht::batched_mht::BatchMerkleTree;
    use algebra::fields::mnt4753::Fr as MNT4753Fr;
    use algebra::fields::mnt6753::Fr as MNT6753Fr;
    use algebra::UniformRand;

    use algebra::biginteger::BigInteger768;

    type MNT4PoseidonBatchHash = PoseidonBatchHash<MNT4753Fr, MNT4753PoseidonParameters>;
    type MNT4BatchedMerkleTree = BatchMerkleTree<MNT4PoseidonBatchHash>;

    type MNT6PoseidonBatchHash = PoseidonBatchHash<MNT6753Fr, MNT6753PoseidonParameters>;
    type MNT6BatchedMerkleTree = BatchMerkleTree<MNT6PoseidonBatchHash>;

    #[test]
    fn merkle_tree_test_mnt4 () {

        let expected_output = MNT4753Fr::new(BigInteger768([8181981188982771303, 9834648934716236448, 6420360685258842467, 14258691490360951478, 10642011566662929522, 16918207755479993617, 3581400602871836321, 14012664850056020974, 16755211538924649257, 4039951447678776727, 12365175056998155257, 119677729692145]));

        let num_leaves = 1024*1024;
        let mut tree = MNT4BatchedMerkleTree::new(num_leaves, 1024 * 1024 / 4);

        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for _ in 0..num_leaves {
             tree.push(MNT4753Fr::rand(&mut rng));
        }
        tree.finalize();
        assert_eq!(tree.root(), expected_output, "Output of the Merkle tree computation for MNT4 does not match to the expected value.");
    }

    #[test]
    fn merkle_tree_test_mnt6 () {

        let expected_output = MNT6753Fr::new(BigInteger768([18065863015580309240, 1059485854425188866, 1479096878827665107, 6899132209183155323, 1829690180552438097, 7395327616910893705, 16132683753083562833, 8528890579558218842, 9345795575555751752, 8161305655297462527, 6222078223269068637, 401142754883827]));

        let num_leaves = 1024*1024;
        let mut tree = MNT6BatchedMerkleTree::new(num_leaves, 1024 * 1024);

        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for _ in 0..num_leaves {
            tree.push(MNT6753Fr::rand(&mut rng));
        }
        tree.finalize();
        assert_eq!(tree.root(), expected_output, "Output of the Merkle tree computation for MNT6 does not match to the expected value.");
    }

}