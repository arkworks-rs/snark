use crate::crh::BatchFieldBasedHash;
use crate::{BatchSBox, CryptoError, Error, PoseidonHash, PoseidonParameters};
use algebra::PrimeField;
use rayon::prelude::*;
use std::marker::PhantomData;

pub struct PoseidonBatchHash<
    F: PrimeField,
    P: PoseidonParameters<Fr = F>,
    SB: BatchSBox<Field = F, Parameters = P>,
> {
    _field: PhantomData<F>,
    _parameters: PhantomData<P>,
    _sbox: PhantomData<SB>,
}

impl<F, P, SB> PoseidonBatchHash<F, P, SB>
where
    F: PrimeField,
    P: PoseidonParameters<Fr = F>,
    SB: BatchSBox<Field = F, Parameters = P>,
{
    fn apply_permutation(input_array: &[F]) -> Vec<Vec<F>> {
        // Sanity checks
        let array_length = input_array.len() / P::R;
        assert_eq!(
            input_array.len() % P::R,
            0,
            "The length of the input data array is not a multiple of the rate."
        );
        assert_ne!(
            input_array.len(),
            0,
            "Input data array does not contain any data."
        );

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
        }

        // Calculate the chunk size to split the state vector
        let cpus = rayon::current_num_threads();
        let chunk_size = (array_length as f64 / cpus as f64).ceil() as usize;

        // apply permutation to different chunks in parallel
        state.par_chunks_mut(chunk_size).for_each(|p1| {
            Self::poseidon_perm_gen(p1);
        });

        state
    }

    fn poseidon_full_round(vec_state: &mut [Vec<P::Fr>], round_cst_idx: &mut usize) {
        // go over each of the state vectors and add the round constants
        for k in 0..vec_state.len() {
            let round_cst_idx_copy = &mut round_cst_idx.clone();
            P::add_round_constants(&mut vec_state[k], round_cst_idx_copy);
        }
        *round_cst_idx += P::T;

        // Apply the S-BOX to each of the elements of the state vector
        SB::apply_full_batch(vec_state);

        // Perform the matrix mix
        for i in 0..vec_state.len() {
            P::matrix_mix(&mut vec_state[i]);
        }
    }

    fn poseidon_partial_round(vec_state: &mut [Vec<P::Fr>], round_cst_idx: &mut usize) {
        // go over each of the state vectors and add the round constants
        for k in 0..vec_state.len() {
            let round_cst_idx_copy = &mut round_cst_idx.clone();
            P::add_round_constants(&mut vec_state[k], round_cst_idx_copy);
        }
        *round_cst_idx += P::T;

        // Apply the S-BOX to the first elements of each of the state vector
        SB::apply_partial_batch(vec_state);

        // Perform the matrix mix
        for i in 0..vec_state.len() {
            P::matrix_mix(&mut vec_state[i]);
        }
    }

    pub fn poseidon_perm_gen(vec_state: &mut [Vec<P::Fr>]) {
        // index that goes over the round constants
        let mut round_cst_idx: usize = 0;

        // Full rounds
        for _i in 0..P::R_F {
            Self::poseidon_full_round(vec_state, &mut round_cst_idx);
        }

        // Partial rounds
        for _i in 0..P::R_P {
            Self::poseidon_partial_round(vec_state, &mut round_cst_idx);
        }

        // Full rounds
        for _i in 0..P::R_F {
            Self::poseidon_full_round(vec_state, &mut round_cst_idx);
        }
    }
}

impl<F, P, SB> BatchFieldBasedHash for PoseidonBatchHash<F, P, SB>
where
    F: PrimeField,
    P: PoseidonParameters<Fr = F>,
    SB: BatchSBox<Field = F, Parameters = P>,
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

        if input_array.len() % P::R != 0 {
            Err(Box::new(CryptoError::Other(
                "The length of the input data array is not a multiple of the rate.".to_owned(),
            )))?
        }

        if input_array.len() == 0 {
            Err(Box::new(CryptoError::Other(
                "Input data array does not contain any data.".to_owned(),
            )))?
        }

        let state = Self::apply_permutation(input_array);

        // write the result of the hash extracted from the state vector to the output vector
        let mut output_array = Vec::new();
        for k in 0..input_array.len() / P::R {
            output_array.push(state[k][0]);
        }
        Ok(output_array)
    }

    fn batch_evaluate_in_place(input_array: &mut [F], output_array: &mut [F]) -> Result<(), Error> {
        // Input:
        // This function calculates the hashes of inputs by groups of the rate P::R.
        // The inputs are arranged in an array and arranged as consecutive chunks
        // Example:
        // For P::R = 2,
        // (d_00, d01, d_10, d_11, d_20, d_21, ...
        // Output:
        // The output will be placed in the output_array taking input length / P::R

        if input_array.len() % P::R != 0 {
            Err(Box::new(CryptoError::Other(
                "The length of the input data array is not a multiple of the rate.".to_owned(),
            )))?
        }

        if input_array.len() == 0 {
            Err(Box::new(CryptoError::Other(
                "Input data array does not contain any data.".to_owned(),
            )))?
        }

        if output_array.len() != input_array.len() / P::R {
            Err(Box::new(CryptoError::Other(format!(
                "Output array size must be equal to input_array_size/rate. Output array size: {}, Input array size: {}, Rate: {}",
                output_array.len(),
                input_array.len(),
                P::R
            ))))?
        }

        let state = Self::apply_permutation(input_array);

        // write the result of the hash extracted from the state vector to the output vector
        for k in 0..input_array.len() / P::R {
            output_array[k] = state[k][0];
        }

        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::{BatchFieldBasedHash, FieldBasedHash};
    use algebra::{Field, UniformRand};
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;
    use std::str::FromStr;

    #[cfg(feature = "mnt4_753")]
    mod mnt4_753 {
        use super::*;
        use crate::{MNT4BatchPoseidonHash, MNT4PoseidonHash};
        use algebra::fields::mnt4753::Fr as MNT4753Fr;

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
                let mut digest = MNT4PoseidonHash::init_constant_length(2, None);
                p.into_iter().for_each(|&f| {
                    digest.update(f);
                });
                output_4753.push(digest.finalize().unwrap());
            });

            // Calculate Poseidon Hash for mnt4753 batch evaluation
            let output_vec = (MNT4BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output_4753[i], output_vec[i],
                    "Hash outputs, position {}, for MNT4 are not equal.",
                    i
                );
            }

            // Check with one single hash
            let single_output = MNT4PoseidonHash::init_constant_length(2, None)
                .update(input_serial[0][0])
                .update(input_serial[0][1])
                .finalize()
                .unwrap();
            let single_batch_output = MNT4BatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

            assert_eq!(
                single_output,
                single_batch_output.unwrap()[0],
                "Single instance hash outputs are not equal for MNT4."
            );
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

            let mut output_vec_in_place = vec![MNT4753Fr::zero(); num_hashes];
            MNT4BatchPoseidonHash::batch_evaluate_in_place(
                &mut input_batch[..],
                &mut output_vec_in_place[..],
            )
            .unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output_vec_in_place[i], output_vec[i],
                    "Hash outputs, position {}, for MNT6 are not equal.",
                    i
                );
            }
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
    }

    #[cfg(feature = "mnt6_753")]
    mod mnt6_753 {
        use super::*;
        use crate::{MNT6BatchPoseidonHash, MNT6PoseidonHash};
        use algebra::fields::mnt6753::Fr as MNT6753Fr;

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
                let mut digest = MNT6PoseidonHash::init_constant_length(2, None);
                p.into_iter().for_each(|&f| {
                    digest.update(f);
                });
                output_6753.push(digest.finalize().unwrap());
            });

            // Calculate Poseidon Hash for mnt4753 batch evaluation
            let output_vec = (MNT6BatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output_6753[i], output_vec[i],
                    "Hash outputs, position {}, for MNT6 are not equal.",
                    i
                );
            }

            // Check with one single hash
            let single_output = MNT6PoseidonHash::init_constant_length(2, None)
                .update(input_serial[0][0])
                .update(input_serial[0][1])
                .finalize()
                .unwrap();
            let single_batch_output = MNT6BatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

            assert_eq!(
                single_output,
                single_batch_output.unwrap()[0],
                "Single instance hash outputs are not equal for MNT6."
            );
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

            let mut output_vec_in_place = vec![MNT6753Fr::zero(); num_hashes];
            MNT6BatchPoseidonHash::batch_evaluate_in_place(
                &mut input_batch[..],
                &mut output_vec_in_place[..],
            )
            .unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output_vec_in_place[i], output_vec[i],
                    "Hash outputs, position {}, for MNT6 are not equal.",
                    i
                );
            }
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
    }

    #[cfg(feature = "bn_382")]
    mod bn_382 {
        use super::*;
        use crate::{
            BN382FqBatchPoseidonHash, BN382FqPoseidonHash, BN382FrBatchPoseidonHash,
            BN382FrPoseidonHash,
        };
        use algebra::fields::bn_382::{Fq as BN382Fq, Fr as BN382Fr};

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
                let mut digest = BN382FqPoseidonHash::init_constant_length(2, None);
                p.into_iter().for_each(|&f| {
                    digest.update(f);
                });
                output.push(digest.finalize().unwrap());
            });

            let output_vec = (BN382FqBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output[i], output_vec[i],
                    "Hash outputs, position {}, for BN382Fq are not equal.",
                    i
                );
            }

            // Check with one single hash
            let single_output = BN382FqPoseidonHash::init_constant_length(2, None)
                .update(input_serial[0][0])
                .update(input_serial[0][1])
                .finalize()
                .unwrap();
            let single_batch_output = BN382FqBatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

            assert_eq!(
                single_output,
                single_batch_output.unwrap()[0],
                "Single instance hash outputs are not equal for BN382Fq."
            );
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
                let mut digest = BN382FrPoseidonHash::init_constant_length(2, None);
                p.into_iter().for_each(|&f| {
                    digest.update(f);
                });
                output.push(digest.finalize().unwrap());
            });

            let output_vec = (BN382FrBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output[i], output_vec[i],
                    "Hash outputs, position {}, for BN382Fr are not equal.",
                    i
                );
            }

            // Check with one single hash
            let single_output = BN382FrPoseidonHash::init_constant_length(2, None)
                .update(input_serial[0][0])
                .update(input_serial[0][1])
                .finalize()
                .unwrap();
            let single_batch_output = BN382FrBatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

            assert_eq!(
                single_output,
                single_batch_output.unwrap()[0],
                "Single instance hash outputs are not equal for BN382Fr."
            );
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

            let mut output_vec_in_place = vec![BN382Fq::zero(); num_hashes];
            BN382FqBatchPoseidonHash::batch_evaluate_in_place(
                &mut input_batch[..],
                &mut output_vec_in_place[..],
            )
            .unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output_vec_in_place[i], output_vec[i],
                    "Hash outputs, position {}, for BN382Fq are not equal.",
                    i
                );
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

            let mut output_vec_in_place = vec![BN382Fr::zero(); num_hashes];
            BN382FrBatchPoseidonHash::batch_evaluate_in_place(
                &mut input_batch[..],
                &mut output_vec_in_place[..],
            )
            .unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output_vec_in_place[i], output_vec[i],
                    "Hash outputs, position {}, for BN382Fr are not equal.",
                    i
                );
            }
        }
    }

    #[cfg(feature = "tweedle")]
    mod tweedle {
        use super::*;
        use crate::{
            TweedleFqBatchPoseidonHash, TweedleFqPoseidonHash, TweedleFrBatchPoseidonHash,
            TweedleFrPoseidonHash,
        };
        use algebra::fields::tweedle::{Fq as TweedleFq, Fr as TweedleFr};

        #[test]
        fn test_batch_hash_tweedlefq() {
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
                let elem1 = TweedleFq::rand(&mut rng);
                let elem2 = TweedleFq::rand(&mut rng);
                pair_elem.push(elem1.clone());
                pair_elem.push(elem2.clone());
                input_serial.push(pair_elem);
                input_batch.push(elem1.clone());
                input_batch.push(elem2.clone());
            }

            // =============================================================================
            let mut output = Vec::new();

            input_serial.iter().for_each(|p| {
                let mut digest = TweedleFqPoseidonHash::init_constant_length(2, None);
                p.into_iter().for_each(|&f| {
                    digest.update(f);
                });
                output.push(digest.finalize().unwrap());
            });

            let output_vec = (TweedleFqBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output[i], output_vec[i],
                    "Hash outputs, position {}, for TweedleFr are not equal.",
                    i
                );
            }

            // Check with one single hash
            let single_output = TweedleFqPoseidonHash::init_constant_length(2, None)
                .update(input_serial[0][0])
                .update(input_serial[0][1])
                .finalize()
                .unwrap();
            let single_batch_output =
                TweedleFqBatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

            assert_eq!(
                single_output,
                single_batch_output.unwrap()[0],
                "Single instance hash outputs are not equal for TweedleFq."
            );
        }

        #[test]
        fn test_batch_hash_tweedlefr() {
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
                let elem1 = TweedleFr::rand(&mut rng);
                let elem2 = TweedleFr::rand(&mut rng);
                pair_elem.push(elem1.clone());
                pair_elem.push(elem2.clone());
                input_serial.push(pair_elem);
                input_batch.push(elem1.clone());
                input_batch.push(elem2.clone());
            }

            // =============================================================================
            let mut output = Vec::new();

            input_serial.iter().for_each(|p| {
                let mut digest = TweedleFrPoseidonHash::init_constant_length(2, None);
                p.into_iter().for_each(|&f| {
                    digest.update(f);
                });
                output.push(digest.finalize().unwrap());
            });

            let output_vec = (TweedleFrBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output[i], output_vec[i],
                    "Hash outputs, position {}, for TweedleFr are not equal.",
                    i
                );
            }

            // Check with one single hash
            let single_output = TweedleFrPoseidonHash::init_constant_length(2, None)
                .update(input_serial[0][0])
                .update(input_serial[0][1])
                .finalize()
                .unwrap();
            let single_batch_output =
                TweedleFrBatchPoseidonHash::batch_evaluate(&input_batch[0..2]);

            assert_eq!(
                single_output,
                single_batch_output.unwrap()[0],
                "Single instance hash outputs are not equal for TweedleFr."
            );
        }

        #[test]
        fn test_batch_hash_tweedlefq_in_place() {
            //  the number of hashes to test
            let num_hashes = 1000;

            // the vectors that store random input data
            let mut input_batch = Vec::new();

            // the random number generator to generate random input data
            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            // we need the double of number of rounds because we have two inputs
            for _ in 0..num_hashes {
                input_batch.push(TweedleFq::rand(&mut rng));
                input_batch.push(TweedleFq::rand(&mut rng));
            }

            // Calculate Poseidon Hash for mnt4753 batch evaluation
            let output_vec = (TweedleFqBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

            let mut output_vec_in_place = vec![TweedleFq::zero(); num_hashes];
            TweedleFqBatchPoseidonHash::batch_evaluate_in_place(
                &mut input_batch[..],
                &mut output_vec_in_place[..],
            )
            .unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output_vec_in_place[i], output_vec[i],
                    "Hash outputs, position {}, for TweedleFq are not equal.",
                    i
                );
            }
        }

        #[test]
        fn test_batch_hash_tweedlefr_in_place() {
            //  the number of hashes to test
            let num_hashes = 1000;

            // the vectors that store random input data
            let mut input_batch = Vec::new();

            // the random number generator to generate random input data
            let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

            // we need the double of number of rounds because we have two inputs
            for _ in 0..num_hashes {
                input_batch.push(TweedleFr::rand(&mut rng));
                input_batch.push(TweedleFr::rand(&mut rng));
            }

            // Calculate Poseidon Hash for mnt4753 batch evaluation
            let output_vec = (TweedleFrBatchPoseidonHash::batch_evaluate(&input_batch)).unwrap();

            let mut output_vec_in_place = vec![TweedleFr::zero(); num_hashes];
            TweedleFrBatchPoseidonHash::batch_evaluate_in_place(
                &mut input_batch[..],
                &mut output_vec_in_place[..],
            )
            .unwrap();

            // =============================================================================
            // Compare results
            for i in 0..num_hashes {
                assert_eq!(
                    output_vec_in_place[i], output_vec[i],
                    "Hash outputs, position {}, for TweedleFr are not equal.",
                    i
                );
            }
        }
    }
}
