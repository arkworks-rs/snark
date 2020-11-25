extern crate rand;
extern crate rayon;

use algebra::PrimeField;
use std::marker::PhantomData;
use crate::crh::BatchFieldBasedHash;
use crate::{Error, PoseidonParameters, PoseidonHash, PoseidonBatchSBox};

pub struct PoseidonBatchHash<F: PrimeField, P: PoseidonParameters<Fr = F>, SB: PoseidonBatchSBox<P>>
{
    _field:      PhantomData<F>,
    _parameters: PhantomData<P>,
    _sbox:       PhantomData<SB>,
}

impl<F, P, SB> PoseidonBatchHash<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonBatchSBox<P>,
{

    fn poseidon_full_round(vec_state: &mut [Vec<P::Fr>], round_cst_idx: &mut usize, last: bool) {

        // For each of the element position of the state vector
        for j in 0..P::T {

            // get the constant associated to element position of the state vector
            let rc = P::ROUND_CST[*round_cst_idx];

            // go over each of the state vectors and add the constant
            for k in 0..vec_state.len() {
                vec_state[k][j] += &rc;
            }
            *round_cst_idx += 1;
        }

        // Apply the S-BOX to each of the elements of the state vector
        SB::apply_full_batch(vec_state, last);

    }

    fn poseidon_partial_round(vec_state: &mut [Vec<P::Fr>], round_cst_idx: &mut usize) {

        // For each of the state vector element position
        for j in 0..P::T {

            // get the constant associated to state vector element position
            let rc = P::ROUND_CST[*round_cst_idx];

            // go over each state vector
            for k in 0..vec_state.len() {
                vec_state[k][j] += &rc;
            }
            *round_cst_idx += 1;
        }

        // Apply the S-BOX to the first elements of each of the state vector
        SB::apply_partial_batch(vec_state);

    }

    pub fn poseidon_perm_gen(vec_state: &mut [Vec<P::Fr>]) {

        // index that goes over the round constants
        let mut round_cst_idx: usize = 0;

        // Full rounds
        for _i in 0..P::R_F {
            Self::poseidon_full_round(vec_state, &mut round_cst_idx, false);
        }

        // Partial rounds
        for _i in 0..P::R_P {
            Self::poseidon_partial_round(vec_state, &mut round_cst_idx);
        }

        // Full rounds
        // Last round does not contain the matrix mix
        for _i in 0..(P::R_F - 1) {
            Self::poseidon_full_round(vec_state, &mut round_cst_idx, false);
        }

        Self::poseidon_full_round(vec_state, &mut round_cst_idx, true);
    }
}


impl<F, P, SB> BatchFieldBasedHash for PoseidonBatchHash<F, P, SB>
    where
        F: PrimeField,
        P: PoseidonParameters<Fr = F>,
        SB: PoseidonBatchSBox<P>,
{
    type Data = F;
    type BaseHash = PoseidonHash<F, P, SB>;

    fn batch_evaluate(input_array: &[F]) -> Result<Vec<F>, Error> {

        // Input:
        // This function calculates the hashes of inputs by groups of the rate P::R.
        // The inputs are arranged in an array and arranged as consecutive chunks
        // Example:
        // For P::R = 2,
        // (d_00, d_01, d_10, d_11, d_20, d_21, ...
        // where d_00 and d_01 are the inputs to the first hash
        // d_10 and d_11 are the inputs to the second hash ... etc..
        // Output:
        // The output is returned as an array of size input length / P::R

        use rayon::prelude::*;

        // Checks that size of input/output vector
        let array_length = input_array.len() / P::R;
        assert_eq!(input_array.len() % P::R, 0, "The length of the input data array is not a multiple of the rate.");
        assert_ne!(input_array.len(), 0, "Input data array does not contain any data.");

        // Assign pre-computed values of the state vector equivalent to a permutation with zero element state vector
        let mut state_z = Vec::new();
        for i in 0..P::T {
            state_z.push(P::AFTER_ZERO_PERM[i]);
        }

        // Copy the result of the permutation to a vector of state vectors of the length equal to the length of the input
        // state is a vector of T-element state vector.
        let mut state = Vec::new();
        for _i in 0..array_length {
            state.push(state_z.clone());
        }

        // input_idx is to scan the input_array
        let mut input_idx = 0;

        // Add the input data in chunks of rate size to the state vectors
        for k in 0..array_length {
            for j in 0..P::R {
                state[k][j] += &input_array[input_idx];
                input_idx += 1;
            }
            // constant for m-ary Merkle tree
            state[k][P::R] += &P::C2;
        }

        // Calculate the chunk size to split the state vector
        let cpus = rayon::current_num_threads();
        let chunk_size = (array_length as f64/ cpus as f64).ceil() as usize;

        // apply permutation to different chunks in parallel
        state.par_chunks_mut(chunk_size)
            .for_each(| p1| {
                Self::poseidon_perm_gen(p1);
            });

        // write the result of the hash extracted from the state vector to the output vector
        let mut output_array=Vec::new();
        for k in 0..array_length {
            output_array.push(state[k][0]);
        }
        Ok(output_array)
    }

    fn batch_evaluate_in_place(input_array: &mut[F], output_array: &mut[F]) {

        // Input:
        // This function calculates the hashes of inputs by groups of the rate P::R.
        // The inputs are arranged in an array and arranged as consecutive chunks
        // Example:
        // For P::R = 2,
        // (d_00, d01, d_10, d_11, d_20, d_21, ...
        // Output:
        // The output will be placed in the output_array taking input length / P::R

        use rayon::prelude::*;

        // Checks that size of input/output vector
        let array_length = input_array.len() / P::R;
        assert_eq!(input_array.len() % P::R, 0, "The length of the input data array is not a multiple of the rate.");
        assert_ne!(input_array.len(), 0, "Input data array does not contain any data.");
        assert_eq!(output_array.len(), array_length,  "The size of the output vector is equal to the size of the input vector divided by the rate.");

        // Assign pre-computed values of the state vector equivalent to a permutation with zero element state vector
        let mut state_z = Vec::new();
        for i in 0..P::T {
            state_z.push(P::AFTER_ZERO_PERM[i]);
        }

        // Copy the result of the permutation to a vector of state vectors of the length equal to the length of the input
        // state is a vector of T-element state vector.
        let mut state = Vec::new();
        for _i in 0..array_length {
            state.push(state_z.clone());
        }

        // input_idx is to scan the input_array
        let mut input_idx = 0;

        for k in 0..array_length {
            for j in 0..P::R {
                state[k][j] += &input_array[input_idx];
                input_idx += 1;
            }
            // constant for m-ary Merkle tree
            state[k][P::R] += &P::C2;
        }

        // Calculate the chunk size to split the state vector
        let cpus = rayon::current_num_threads();
        let chunk_size = (array_length as f64/ cpus as f64).ceil() as usize;

        // apply permutation to different chunks in parallel
        state.par_chunks_mut(chunk_size)
            .for_each(| p1| {
                Self::poseidon_perm_gen(p1);
            });

        // write the result of the hash extracted from the state vector to the output vector
        for k in 0..array_length {
            output_array[k] = state[k][0];
        }
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use rand_xorshift::XorShiftRng;
    use std::str::FromStr;
    use crate::{FieldBasedHash, BatchFieldBasedHash};
    use super::rand::SeedableRng;
    use algebra::{
        fields::{
            mnt6753::Fr as MNT6753Fr,
            mnt4753::Fr as MNT4753Fr,
            bn_382::Fq as BN382Fq,
            bn_382::Fr as BN382Fr,
        },
        UniformRand
    };

    use crate::crh::poseidon::{
        MNT4PoseidonHash, MNT6PoseidonHash,
        MNT4BatchPoseidonHash, MNT6BatchPoseidonHash,
        BN382FqPoseidonHash, BN382FrPoseidonHash,
        BN382FqBatchPoseidonHash, BN382FrBatchPoseidonHash,
        parameters::{
            MNT4753PoseidonParameters, MNT6753PoseidonParameters,
            BN382FqPoseidonParameters, BN382FrPoseidonParameters
        }
    };

    #[test]
    fn test_batch_hash_mnt4() {

        //  the number of hashes to test
        let num_hashes = 1000;

        // the vectors that store random input data
        let mut input_serial = Vec::new();
        let mut input_batch = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_hashes {
            let mut pair_elem = Vec::new();
            let elem1 = MNT4753Fr::rand(&mut rng);
            let elem2 = MNT4753Fr::rand(&mut rng);
            pair_elem.push(elem1.clone());
            pair_elem.push(elem2.clone());
            input_serial.push(pair_elem);
            input_batch.push(elem1.clone());
            input_batch.push(elem2.clone());
        }

        // =============================================================================
        // Calculate Poseidon Hash for mnt4753
        let mut output_4753 = Vec::new();

        input_serial.iter().for_each(|p| {
            let mut digest = MNT4PoseidonHash::init(None);
            p.into_iter().for_each(|&f| { digest.update(f); });
            output_4753.push(digest.finalize());
        });

        // Calculate Poseidon Hash for mnt4753 batch evaluation
        let output_vec = (MNT4BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

        // =============================================================================
        // Compare results
        for i in 0..num_hashes {
            assert_eq!(output_4753[i], output_vec[i], "Hash outputs, position {}, for MNT4 are not equal.", i);
        }

        // Check with one single hash
        let single_output = MNT4PoseidonHash::init(None)
            .update(input_serial[0][0])
            .update(input_serial[0][1])
            .finalize();
        let single_batch_output = MNT4BatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

        assert_eq!(single_output, single_batch_output.unwrap()[0], "Single instance hash outputs are not equal for MNT4.");
    }


    #[test]
    fn test_batch_hash_bn382fq() {

        //  the number of hashes to test
        let num_hashes = 1000;

        // the vectors that store random input data
        let mut input_serial = Vec::new();
        let mut input_batch = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_hashes {
            let mut pair_elem = Vec::new();
            let elem1 = BN382Fq::rand(&mut rng);
            let elem2 = BN382Fq::rand(&mut rng);
            pair_elem.push(elem1.clone());
            pair_elem.push(elem2.clone());
            input_serial.push(pair_elem);
            input_batch.push(elem1.clone());
            input_batch.push(elem2.clone());
        }

        // =============================================================================
        let mut output = Vec::new();

        input_serial.iter().for_each(|p| {
            let mut digest = BN382FqPoseidonHash::init(None);
            p.into_iter().for_each(|&f| { digest.update(f); });
            output.push(digest.finalize());
        });

        let output_vec = (BN382FqBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

        // =============================================================================
        // Compare results
        for i in 0..num_hashes {
            assert_eq!(output[i], output_vec[i], "Hash outputs, position {}, for BN382Fr are not equal.", i);
        }

        // Check with one single hash
        let single_output = BN382FqPoseidonHash::init(None)
            .update(input_serial[0][0])
            .update(input_serial[0][1])
            .finalize();
        let single_batch_output = BN382FqBatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

        assert_eq!(single_output, single_batch_output.unwrap()[0], "Single instance hash outputs are not equal for BN382Fr.");
    }


    #[test]
    #[should_panic]
    fn test_batch_hash_mnt4_null_elem() {
        let input_batch = Vec::new();
        let output_vec = (MNT4BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();
        println!("{:?}", output_vec);
    }

    #[test]
    #[should_panic]
    fn test_batch_hash_mnt4_one_elem() {
        let mut input_batch = Vec::new();
        input_batch.push(MNT4753Fr::from_str("1").unwrap());
        let output_vec = (MNT4BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();
        println!("{:?}", output_vec);
    }

    #[test]
    #[should_panic]
    fn test_batch_hash_mnt4_three_elem() {
        let mut input_batch = Vec::new();
        input_batch.push(MNT4753Fr::from_str("1").unwrap());
        input_batch.push(MNT4753Fr::from_str("2").unwrap());
        input_batch.push(MNT4753Fr::from_str("3").unwrap());
        let output_vec = (MNT4BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();
        println!("{:?}", output_vec);
    }

    #[test]
    fn test_batch_hash_mnt6() {

        //  the number of hashes to test
        let num_hashes = 1000;

        // the vectors that store random input data
        let mut input_serial = Vec::new();
        let mut input_batch = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_hashes {
            let mut pair_elem = Vec::new();
            let elem1 = MNT6753Fr::rand(&mut rng);
            let elem2 = MNT6753Fr::rand(&mut rng);
            pair_elem.push(elem1.clone());
            pair_elem.push(elem2.clone());
            input_serial.push(pair_elem);
            input_batch.push(elem1.clone());
            input_batch.push(elem2.clone());
        }

        // =============================================================================
        // Calculate Poseidon Hash for mnt6753
        let mut output_6753 = Vec::new();

        input_serial.iter().for_each(|p| {
            let mut digest = MNT6PoseidonHash::init(None);
            p.into_iter().for_each(|&f| { digest.update(f); });
            output_6753.push(digest.finalize());
        });

        // Calculate Poseidon Hash for mnt4753 batch evaluation
        let output_vec = (MNT6BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

        // =============================================================================
        // Compare results
        for i in 0..num_hashes {
            assert_eq!(output_6753[i], output_vec[i], "Hash outputs, position {}, for MNT6 are not equal.", i);
        }

        // Check with one single hash
        let single_output = MNT6PoseidonHash::init(None)
            .update(input_serial[0][0])
            .update(input_serial[0][1])
            .finalize();
        let single_batch_output = MNT6BatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

        assert_eq!(single_output, single_batch_output.unwrap()[0], "Single instance hash outputs are not equal for MNT6.");
    }

    #[test]
    fn test_batch_hash_bn382fr() {

        //  the number of hashes to test
        let num_hashes = 1000;

        // the vectors that store random input data
        let mut input_serial = Vec::new();
        let mut input_batch = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_hashes {
            let mut pair_elem = Vec::new();
            let elem1 = BN382Fr::rand(&mut rng);
            let elem2 = BN382Fr::rand(&mut rng);
            pair_elem.push(elem1.clone());
            pair_elem.push(elem2.clone());
            input_serial.push(pair_elem);
            input_batch.push(elem1.clone());
            input_batch.push(elem2.clone());
        }

        // =============================================================================
        let mut output = Vec::new();

        input_serial.iter().for_each(|p| {
            let mut digest = BN382FrPoseidonHash::init(None);
            p.into_iter().for_each(|&f| { digest.update(f); });
            output.push(digest.finalize());
        });

        let output_vec = (BN382FrBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

        // =============================================================================
        // Compare results
        for i in 0..num_hashes {
            assert_eq!(output[i], output_vec[i], "Hash outputs, position {}, for BN382Fr are not equal.", i);
        }

        // Check with one single hash
        let single_output = BN382FrPoseidonHash::init(None)
            .update(input_serial[0][0])
            .update(input_serial[0][1])
            .finalize();
        let single_batch_output = BN382FrBatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

        assert_eq!(single_output, single_batch_output.unwrap()[0], "Single instance hash outputs are not equal for BN382Fr.");
    }


    #[test]
    #[should_panic]
    fn test_batch_hash_mnt6_null_elem() {
        let input_batch = Vec::new();
        let output_vec = (MNT6BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();
        println!("{:?}", output_vec);
    }

    #[test]
    #[should_panic]
    fn test_batch_hash_mnt6_one_elem() {
        let mut input_batch = Vec::new();
        input_batch.push(MNT6753Fr::from_str("1").unwrap());
        let output_vec = (MNT6BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();
        println!("{:?}", output_vec);
    }

    #[test]
    #[should_panic]
    fn test_batch_hash_mnt6_three_elem() {
        let mut input_batch = Vec::new();
        input_batch.push(MNT6753Fr::from_str("1").unwrap());
        input_batch.push(MNT6753Fr::from_str("2").unwrap());
        input_batch.push(MNT6753Fr::from_str("3").unwrap());
        let output_vec = (MNT6BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();
        println!("{:?}", output_vec);
    }

    #[test]
    fn test_batch_hash_mnt4_in_place() {
        //  the number of hashes to test
        let num_hashes = 1000;

        // the vectors that store random input data
        let mut input_batch = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_hashes {
            input_batch.push(MNT4753Fr::rand(&mut rng));
            input_batch.push(MNT4753Fr::rand(&mut rng));
        }

        // Calculate Poseidon Hash for mnt4753 batch evaluation
        let output_vec = (MNT4BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

        let mut output_vec_in_place = vec![MNT4753PoseidonParameters::ZERO; num_hashes];
        MNT4BatchPoseidonHash::batch_evaluate_in_place(&mut input_batch[..], &mut output_vec_in_place[..]);

        // =============================================================================
        // Compare results
        for i in 0..num_hashes {
            assert_eq!(output_vec_in_place[i], output_vec[i], "Hash outputs, position {}, for MNT6 are not equal.", i);
        }
    }

    #[test]
    fn test_batch_hash_mnt6_in_place() {

        //  the number of hashes to test
        let num_hashes = 1000;

        // the vectors that store random input data
        let mut input_batch = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_hashes {
            input_batch.push(MNT6753Fr::rand(&mut rng));
            input_batch.push(MNT6753Fr::rand(&mut rng));
        }

        // Calculate Poseidon Hash for mnt4753 batch evaluation
        let output_vec = (MNT6BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

        let mut output_vec_in_place = vec![MNT6753PoseidonParameters::ZERO; num_hashes];
        MNT6BatchPoseidonHash::batch_evaluate_in_place(&mut input_batch[..], &mut output_vec_in_place[..]);

        // =============================================================================
        // Compare results
        for i in 0..num_hashes {
            assert_eq!(output_vec_in_place[i], output_vec[i], "Hash outputs, position {}, for MNT6 are not equal.", i);
        }
    }

    #[test]
    fn test_batch_hash_bn382fq_in_place() {
        //  the number of hashes to test
        let num_hashes = 1000;

        // the vectors that store random input data
        let mut input_batch = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_hashes {
            input_batch.push(BN382Fq::rand(&mut rng));
            input_batch.push(BN382Fq::rand(&mut rng));
        }

        // Calculate Poseidon Hash for mnt4753 batch evaluation
        let output_vec = (BN382FqBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

        let mut output_vec_in_place = vec![BN382FqPoseidonParameters::ZERO; num_hashes];
        BN382FqBatchPoseidonHash::batch_evaluate_in_place(&mut input_batch[..], &mut output_vec_in_place[..]);

        // =============================================================================
        // Compare results
        for i in 0..num_hashes {
            assert_eq!(output_vec_in_place[i], output_vec[i], "Hash outputs, position {}, for BN382Fq are not equal.", i);
        }
    }

    #[test]
    fn test_batch_hash_bn382fr_in_place() {

        //  the number of hashes to test
        let num_hashes = 1000;

        // the vectors that store random input data
        let mut input_batch = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // we need the double of number of rounds because we have two inputs
        for _ in 0..num_hashes {
            input_batch.push(BN382Fr::rand(&mut rng));
            input_batch.push(BN382Fr::rand(&mut rng));
        }

        // Calculate Poseidon Hash for mnt4753 batch evaluation
        let output_vec = (BN382FrBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

        let mut output_vec_in_place = vec![BN382FrPoseidonParameters::ZERO; num_hashes];
        BN382FrBatchPoseidonHash::batch_evaluate_in_place(&mut input_batch[..], &mut output_vec_in_place[..]);

        // =============================================================================
        // Compare results
        for i in 0..num_hashes {
            assert_eq!(output_vec_in_place[i], output_vec[i], "Hash outputs, position {}, for BN382Fr are not equal.", i);
        }
    }
}