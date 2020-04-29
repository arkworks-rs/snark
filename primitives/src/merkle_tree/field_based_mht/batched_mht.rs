extern crate rand;

use crate::Error;

use crate::crh::{
    BatchFieldBasedHash,
    FieldBasedHashParameters, poseidon::{
        parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters}
    }
};
use algebra::Field;

pub struct BatchMerkleTree<H: BatchFieldBasedHash>{
    root: H::Data
}

impl<H: BatchFieldBasedHash> BatchMerkleTree<H> {

    pub fn root(&self) -> H::Data {
        self.root.clone()
    }

    fn merkle_tree_compute(input_vec: &mut[H::Data], output_vec: &mut[H::Data], input_size: usize){
        // Supporting function that processes the inputs and outputs in chunks

        use rayon::prelude::*;

        let num_cores = 16;

        if input_size < 2 * num_cores {
            input_vec.par_chunks_mut(2).zip(output_vec.par_chunks_mut(1)).for_each( |(p1,p2)| {
                H::batch_evaluate_in_place(p1, p2);
            });
            return;
        }

        input_vec.par_chunks_mut(input_size/num_cores).zip(output_vec.par_chunks_mut(input_size/num_cores/2)).for_each( |(p1, p2)| {
            H::batch_evaluate_in_place(p1, p2);
        });

    }

    pub fn new(mut input_vec: &mut [H::Data]) -> Self {

        // rate = 2
        // d01, d02, d11, d12, d21, d22

        let last_level_size = input_vec.len().next_power_of_two();
        let tree_size = 2 * last_level_size - 1;
        let size_output = tree_size - input_vec.len();

        let mut output = Vec::new();
        for i in 0..size_output {
            output.push(H::Data::zero());
        }

        let mut size_input = input_vec.len();

        let mut copy_output = &mut output[..];
        Self::merkle_tree_compute(input_vec, copy_output, size_input);
        size_input /= 2;

        while size_input > 1 {
            let (input_vec, output_vec) = copy_output.split_at_mut(size_input);
            Self::merkle_tree_compute(input_vec, output_vec, size_input);
            copy_output = output_vec;
            size_input /= 2;
        }

        Self{ root: {*copy_output.last().unwrap()}}
    }

}

#[cfg(test)]
mod test {
    use rand_xorshift::XorShiftRng;
    use super::rand::SeedableRng;
    use crate::PoseidonBatchHash;
    use crate::crh::poseidon::parameters::MNT4753PoseidonParameters;
    use crate::merkle_tree::field_based_mht::batched_mht::BatchMerkleTree;
    use algebra::fields::mnt4753::Fr as MNT4753Fr;
    use algebra::fields::mnt6753::Fr as MNT6753Fr;
    use algebra::UniformRand;

    type MNT4PoseidonBatchHash = PoseidonBatchHash<MNT4753Fr, MNT4753PoseidonParameters>;
    type MNT4BatchedMerkleTree = BatchMerkleTree<MNT4PoseidonBatchHash>;

    #[test]
    fn merkle_tree_test () {

        //  the number of leaves to test
        let num_leaves = 1024;

        // the vectors that store random input data
        let mut leaves = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_leaves {
            leaves.push(MNT4753Fr::rand(&mut rng));
        }

        let tree = MNT4BatchedMerkleTree::new(&mut leaves[..]);
        let root = tree.root();

        println!("Merkle root = {:?}", root);
    }

}