extern crate hex;
extern crate rand;
extern crate rayon;

use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::mnt4753::Fr as MNT4753Fr;

use crate::crh::{
    FieldBasedHashParameters, poseidon::{
        poseidon_original::PoseidonHash,
        parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters}
    }
};

pub mod poseidon_original;
pub mod parameters;

pub trait PoseidonParameters: 'static + FieldBasedHashParameters{

    const T: usize;  // Number of S-Boxesb
    const R_F:i32;   // Number of full rounds
    const R_P:i32;   // Number of partial rounds
    const R:usize;   // The rate of the hash function
    const ZERO:Self::Fr;   // The zero element in the field
    const C2:Self::Fr;     // The constant 3 to add in the position corresponding to the capacity
    const AFTER_ZERO_PERM: &'static[Self::Fr]; // State vector after a zero permutation
    const ROUND_CST: &'static[Self::Fr];  // Array of round constants
    const MDS_CST: &'static[Self::Fr];  // The MDS matrix
    const MDS_CST_SHORT: &'static[Self::Fr];  // The MDS matrix for fast matrix multiplication

}

pub type MNT4PoseidonHash = PoseidonHash<MNT4753Fr, MNT4753PoseidonParameters>;
pub type MNT6PoseidonHash = PoseidonHash<MNT6753Fr, MNT6753PoseidonParameters>;

#[cfg(test)]
mod test {
    use super::*;
    use rayon::prelude::*;
    use rand_xorshift::XorShiftRng;
    use std::str::FromStr;
    use crate::{FieldBasedHash, BatchFieldBasedHash};
    use crate::crh::poseidon::poseidon_original::PoseidonBatchHash;
    use super::rand::SeedableRng;
    use algebra::UniformRand;
    use std::time::Instant;

    #[test]
    fn test_poseidon_hash_mnt4() {
        let mut input = Vec::new();
        input.push(MNT4753Fr::from_str("1").unwrap());
        input.push(MNT4753Fr::from_str("2").unwrap());
        let output = MNT4PoseidonHash::evaluate(&input);

        println!("{:?}", output);
    }


    #[test]
    fn test_poseidon_hash_mnt6() {
        let mut input = Vec::new();
        input.push(MNT6753Fr::from_str("1").unwrap());
        input.push(MNT6753Fr::from_str("2").unwrap());
        let output = MNT6PoseidonHash::evaluate(&mut input);

        println!("{:?}", output);
    }

    #[test]
    fn test_hash_speed() {
        // =============================================================================
        // Computation for MNT4

        type Mnt4PoseidonHash = PoseidonHash<MNT4753Fr, MNT4753PoseidonParameters>;
        type Mnt4BatchPoseidonHash = PoseidonBatchHash<MNT4753Fr, MNT4753PoseidonParameters>;

        //  the number of rounds to test
        let num_rounds = 1000;

        // the vectors that store random input data
        let mut vec_vec_elem_4753 = Vec::new();

        let mut array_elem_4753 = Vec::new();

        // the random number generator to generate random input data
        // let mut rng = &mut thread_rng();
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_rounds {
            let mut vec_elem_4753 = Vec::new();
            let elem1 = MNT4753Fr::rand(&mut rng);
            let elem2 = MNT4753Fr::rand(&mut rng);
            vec_elem_4753.push(elem1.clone());
            vec_elem_4753.push(elem2.clone());
            vec_vec_elem_4753.push(vec_elem_4753);
            array_elem_4753.push(elem1.clone());
            array_elem_4753.push(elem2.clone());
        }

        // =============================================================================
        // Calculate Poseidon Hash for mnt4753
        let now_4753 = Instant::now();

        let mut output_4753 = Vec::new();

        for i in 0..num_rounds {

            // Call the poseidon hash
            let output = Mnt4PoseidonHash::evaluate(&vec_vec_elem_4753[i]);
            output_4753.push(output.unwrap());
        }

        let new_now_4753 = Instant::now();

        // =============================================================================
        // Calculate Poseidon Hash for mnt4753 batch evaluation

        let mut array1 = Vec::new();
        let mut array2 = Vec::new();
        let mut array3 = Vec::new();
        let mut array4 = Vec::new();

        for i in 0..(num_rounds / 4) {
            array1.push(vec_vec_elem_4753[i][0].clone());
            array1.push(vec_vec_elem_4753[i][1].clone());
        }
        for i in (num_rounds / 4)..(num_rounds / 2) {
            array2.push(vec_vec_elem_4753[i][0].clone());
            array2.push(vec_vec_elem_4753[i][1].clone());
        }
        for i in (num_rounds / 2)..(num_rounds * 3 / 4) {
            array3.push(vec_vec_elem_4753[i][0].clone());
            array3.push(vec_vec_elem_4753[i][1].clone());
        }
        for i in (num_rounds * 3 / 4)..(num_rounds) {
            array4.push(vec_vec_elem_4753[i][0].clone());
            array4.push(vec_vec_elem_4753[i][1].clone());
        }

        let mut array_array_input = Vec::new();
        array_array_input.push(array1);
        array_array_input.push(array2);
        array_array_input.push(array3);
        array_array_input.push(array4);


        let now_4753_batch = Instant::now();

        array_array_input.par_iter_mut().for_each(|mut p| Mnt4BatchPoseidonHash::batch_evaluate_2_1(&mut p));

        let new_now_4753_batch = Instant::now();

        // Call the poseidon batch hash
        let mut output_4753_batch = Vec::new();

        for i in 0..num_rounds / 4 {
            output_4753_batch.push(array_array_input[0][i]);
        }
        for i in 0..num_rounds / 4 {
            output_4753_batch.push(array_array_input[1][i]);
        }
        for i in 0..num_rounds / 4 {
            output_4753_batch.push(array_array_input[2][i]);
        }
        for i in 0..num_rounds / 4 {
            output_4753_batch.push(array_array_input[3][i]);
        }

        // =============================================================================
        // Compare results
        let output_batch = output_4753_batch;
        for i in 0..num_rounds {
            if output_4753[i] != output_batch[i] {
                println!("Hash outputs, position {}, for MNT4 are not equal.", i);
            }
        }
        println!("End comparison for MNT4.");

        // =============================================================================
        // Report the timing results

        let duration_4753_single = new_now_4753.duration_since(now_4753);
        println!("Time for {} rounds MNT4753 single = {:?}", num_rounds, duration_4753_single.as_millis());

        let duration_4753_batch = new_now_4753_batch.duration_since(now_4753_batch);
        println!("Time for {} rounds MNT4753 batch = {:?}", num_rounds, duration_4753_batch.as_millis());

        // =============================================================================
        // Computation for MNT6

        type Mnt6PoseidonHash = PoseidonHash<MNT6753Fr, MNT6753PoseidonParameters>;
        type Mnt6BatchPoseidonHash = PoseidonBatchHash<MNT6753Fr, MNT6753PoseidonParameters>;

        // the vectors that store random input data
        let mut vec_vec_elem_6753 = Vec::new();

        let mut array_elem_6753 = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_rounds {
            let mut vec_elem_6753 = Vec::new();
            let elem1 = MNT6753Fr::rand(&mut rng);
            let elem2 = MNT6753Fr::rand(&mut rng);
            vec_elem_6753.push(elem1.clone());
            vec_elem_6753.push(elem2.clone());
            vec_vec_elem_6753.push(vec_elem_6753);
            array_elem_6753.push(elem1.clone());
            array_elem_6753.push(elem2.clone());
        }

        // =============================================================================
        // Calculate Poseidon Hash for mnt6753
        let now_6753 = Instant::now();

        let mut output_6753 = Vec::new();

        for i in 0..num_rounds {

            // Call the poseidon hash
            let output = Mnt6PoseidonHash::evaluate(&vec_vec_elem_6753[i]);
            output_6753.push(output.unwrap());
        }

        let new_now_6753 = Instant::now();

        // =============================================================================
        // Calculate Poseidon Hash for mnt6753 batch evaluation

        let mut array1 = Vec::new();
        let mut array2 = Vec::new();
        let mut array3 = Vec::new();
        let mut array4 = Vec::new();

        for i in 0..(num_rounds / 4) {
            array1.push(vec_vec_elem_6753[i][0].clone());
            array1.push(vec_vec_elem_6753[i][1].clone());
        }
        for i in (num_rounds / 4)..(num_rounds / 2) {
            array2.push(vec_vec_elem_6753[i][0].clone());
            array2.push(vec_vec_elem_6753[i][1].clone());
        }
        for i in (num_rounds / 2)..(num_rounds * 3 / 4) {
            array3.push(vec_vec_elem_6753[i][0].clone());
            array3.push(vec_vec_elem_6753[i][1].clone());
        }
        for i in (num_rounds * 3 / 4)..(num_rounds) {
            array4.push(vec_vec_elem_6753[i][0].clone());
            array4.push(vec_vec_elem_6753[i][1].clone());
        }

        let mut array_array_input = Vec::new();
        array_array_input.push(array1);
        array_array_input.push(array2);
        array_array_input.push(array3);
        array_array_input.push(array4);


        let now_6753_batch = Instant::now();

        array_array_input.par_iter_mut().for_each(|mut p| Mnt6BatchPoseidonHash::batch_evaluate_2_1(&mut p));

        let new_now_6753_batch = Instant::now();

        // Call the poseidon batch hash
        let mut output_6753_batch = Vec::new();

        for i in 0..num_rounds / 4 {
            output_6753_batch.push(array_array_input[0][i]);
        }
        for i in 0..num_rounds / 4 {
            output_6753_batch.push(array_array_input[1][i]);
        }
        for i in 0..num_rounds / 4 {
            output_6753_batch.push(array_array_input[2][i]);
        }
        for i in 0..num_rounds / 4 {
            output_6753_batch.push(array_array_input[3][i]);
        }

        // =============================================================================
        // Compare results
        let output_batch = output_6753_batch;
        for i in 0..num_rounds {
            if output_6753[i] != output_batch[i] {
                println!("Hash outputs, position {}, for MNT6 are not equal.", i);
            }
        }
        println!("End comparison for MNT6.");

        // =============================================================================
        // Report the timing results

        let duration_6753_single = new_now_6753.duration_since(now_6753);
        println!("Time for {} rounds MNT6753 single = {:?}", num_rounds, duration_6753_single.as_millis());

        let duration_6753_batch = new_now_6753_batch.duration_since(now_6753_batch);
        println!("Time for {} rounds MNT6753 batch = {:?}", num_rounds, duration_6753_batch.as_millis());

    }

}