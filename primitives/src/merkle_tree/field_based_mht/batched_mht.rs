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
        //
        // let num_cores = 16;
        //
        // assert_eq!(input_vec.len() % 2, 0, "The length of the input to the hash is not even.");
        // assert_eq!(output_vec.len() >= input_vec.len() / 2, true,  "The length of the output is not greater or equal to half of the input length.");
        //
        // if input_size < 2 * num_cores {
        //     input_vec.par_chunks_mut(2).zip(output_vec.par_chunks_mut(1)).for_each( |(p1,p2)| {
        //         H::batch_evaluate(p1, p2);
        //     });
        //     return;
        // }
        //
        // use rayon::prelude::*;
        //
        // input_vec.par_chunks_mut(input_size/num_cores).zip(output_vec.par_chunks_mut(input_size/num_cores/2)).for_each( |(p1, p2)| {
        //     H::batch_evaluate(p1, p2);
        // });

    }

    pub fn new(array: &mut [H::Data]) -> Self {

        // rate = 2
        // d01, d02, d11, d12, d21, d22

        /*
        let mut size_input = array.len();
        size_input /= 2;
        let mut size_output = size_input;
        while (size_input > 1) {
            size_input /= 2;
            size_output += size_input;
        }
        */

        //let last_level_size = leaves.len().next_power_of_two();
        // let tree size = 2* last_level_size -1;
        // /* 2 ^ height -1 */
        //
        //
        // let mut output = Vec::with_capacity(tree_size);
        // let mut copy_vec = &mut array[..];
        // let mut size_input = input_size;
        //
        // while size_input > 1{
        //     Self::merkle_tree_compute(input_vec, output, size_input);
        //     copy_vec = output;
        //     let (input_vec, output) = copy_vec.split_at_mut(size_input);
        //     size_input = size_input / 2;
        // }
        let output = vec![H::Data::zero();1];

        Self{ root: {*output.last().unwrap()}}
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
            leaves.push(MNT4753Fr::rand(&mut rng));
        }

        let tree = MNT4BatchedMerkleTree::new(&mut leaves[..]);
        let root = tree.root();

        println!("Merkle root = {:?}", root);
    }

}