extern crate rand;

use crate::crh::BatchFieldBasedHash;
use algebra::{PrimeField, MulShort};
use crate::PoseidonParameters;
use crate::merkle_tree::field_based_mht::ramht::{RandomAccessMerkleTree, RandomAccessDeleteMerkleTree, NullableLeaf};
use crate::crh::poseidon::batched_crh::PoseidonBatchHash;
use std::marker::PhantomData;
use std::clone::Clone;

#[derive(Clone, Debug)]
pub struct PoseidonRAMTPath<F: PrimeField + MulShort>(
    pub Vec<(F, bool)>
);

#[derive(Clone)]
pub struct PoseidonRandomAccessMerkleTree<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>>{
    root: F,
    array_nodes: Vec<F>,
    processing_step: usize,
    final_pos: Vec<usize>,
    initial_pos: Vec<usize>,
    initial_pos_subarray: Vec<usize>,
    final_pos_subarray: Vec<usize>,
    processed_pos: Vec<usize>,
    processed_pos_subarray: Vec<usize>,
    new_elem_pos: Vec<usize>,
    new_elem_pos_subarray: Vec<usize>,
    levels: usize,
    rate: usize,
    num_leaves: usize,
    finalized: bool,
    _parameters: PhantomData<P>
}

impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> PoseidonRandomAccessMerkleTree<F, P> {

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

impl<F: PrimeField + MulShort, P: Clone + PoseidonParameters<Fr=F>> RandomAccessMerkleTree
for PoseidonRandomAccessMerkleTree<F, P> {

    type Data = F;
    type MerklePath = PoseidonRAMTPath<F>;

    fn init(num_leaves: usize) -> Self {
        // change it to the rate of the hash function
        let rate = P::R;
        let processing_step = num_leaves;

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

        let last_level_size = num_leaves.next_power_of_two();
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
        // initial position of the full array
        let mut initial_pos = Vec::new();


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
            initial_pos.push(initial_idx);

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
            initial_pos: {initial_pos},
            initial_pos_subarray: { initial_pos_subarray },
            final_pos_subarray: { final_pos_subarray },
            processed_pos: { processed_pos },
            processed_pos_subarray: { processed_pos_subarray },
            new_elem_pos: { new_elem_pos },
            new_elem_pos_subarray: { new_elem_pos_subarray },
            levels: { level_idx - 2 },
            rate: { P::R },
            num_leaves: { last_level_size },
            finalized: false,
            _parameters: PhantomData
        }
    }

    fn reset(&mut self) -> &mut Self {

        for i in 0..self.new_elem_pos.len() {
            self.new_elem_pos_subarray[i] = self.initial_pos_subarray[i];
            self.new_elem_pos[i] = self.initial_pos[i];
            self.processed_pos[i] = self.initial_pos[i];
            self.processed_pos_subarray[i] = self.initial_pos_subarray[i];
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

        self
    }

    fn set(&mut self, index: usize, new_leaf: Self::Data) -> &mut Self {
        unimplemented!()
    }

    fn finalize(&self) -> Self {
        let mut copy = (*self).clone();
        while copy.new_elem_pos[0] != copy.final_pos[0] {
            copy.append(F::zero());
        }
        copy.compute_subtree();
        copy.root = *copy.array_nodes.last().unwrap();
        copy.finalized = true;
        copy
    }

    fn finalize_in_place(&mut self) -> &mut Self {
        while self.new_elem_pos[0] != self.final_pos[0] {
            self.append(F::zero());
        }
        self.compute_subtree();
        self.root = (*self.array_nodes.last().unwrap()).clone();
        self.finalized = true;

        self
    }

    fn root(&self) -> Option<Self::Data> {
        match self.finalized {
            true => Some(self.root),
            false => None
        }
    }

    fn get_merkle_path(&self, leaf: &Self::Data) -> Option<Self::MerklePath> {
        unimplemented!()
    }

    fn verify_merkle_path(
        &self,
        leaf: &Self::Data,
        root: &Self::Data,
        path: &Self::MerklePath
    ) -> Option<bool> {
        unimplemented!()
    }
}

pub struct NullablePoseidonRADMTLeaf<F: PrimeField + MulShort>(PhantomData<F>);

impl<F: PrimeField + MulShort> NullableLeaf for NullablePoseidonRADMTLeaf<F>{
    type Data = F;

    fn get_null_value() -> Self::Data {
        F::zero()
    }
}

impl<F: PrimeField + MulShort, P: Clone + PoseidonParameters<Fr=F>> RandomAccessDeleteMerkleTree
for PoseidonRandomAccessMerkleTree<F, P> {
    type Leaf = NullablePoseidonRADMTLeaf<F>;
}


#[cfg(test)]
mod test {
    //TODO: Add test cases for new functions
    use rand_xorshift::XorShiftRng;
    use super::rand::SeedableRng;
    use crate::crh::poseidon::parameters::{
        MNT4753PoseidonParameters, MNT6753PoseidonParameters,
    };
    use algebra::fields::mnt4753::Fr as MNT4753Fr;
    use algebra::fields::mnt6753::Fr as MNT6753Fr;
    use algebra::UniformRand;

    use algebra::biginteger::BigInteger768;
    use crate::merkle_tree::field_based_mht::ramht::poseidon::PoseidonRandomAccessMerkleTree;
    use crate::merkle_tree::field_based_mht::ramht::RandomAccessMerkleTree;

    type MNT4PoseidonRAMT = PoseidonRandomAccessMerkleTree<MNT4753Fr, MNT4753PoseidonParameters>;
    type MNT6PoseidonRAMT = PoseidonRandomAccessMerkleTree<MNT6753Fr, MNT6753PoseidonParameters>;

    #[test]
    fn merkle_tree_test_mnt4() {
        // running time for 1048576 leaves
        // processing_step = 1024 => 90278 ms
        // processing_step = 1024 * 64 => 40753 ms
        // processing_step = 1024 * 1024 => 38858 ms

        let expected_output = MNT4753Fr::new(BigInteger768([8181981188982771303, 9834648934716236448, 6420360685258842467, 14258691490360951478, 10642011566662929522, 16918207755479993617, 3581400602871836321, 14012664850056020974, 16755211538924649257, 4039951447678776727, 12365175056998155257, 119677729692145]));
        let num_leaves = 1024 * 1024;
        let mut tree = MNT4PoseidonRAMT::init(num_leaves);
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
        let mut tree = MNT6PoseidonRAMT::init(num_leaves);
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
        for _ in 0..num_leaves {
            tree.append(MNT6753Fr::rand(&mut rng));
        }
        tree.finalize_in_place();
        assert_eq!(tree.root().unwrap(), expected_output, "Output of the Merkle tree computation for MNT6 does not match to the expected value.");
    }

}
