extern crate hex;
extern crate rand;
extern crate rayon;
extern crate num_cpus;

use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::{PrimeField, MulShort};

use std::marker::PhantomData;

use crate::crh::{
    FieldBasedHashParameters, poseidon::{
        parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters}
    }
};

use crate::crh::{FieldBasedHash, BatchFieldBasedHash};
use crate::Error;

pub mod parameters;

pub struct PoseidonHash<F: PrimeField, P: PoseidonParameters<Fr = F>>{
    _field:      PhantomData<F>,
    _parameters: PhantomData<P>,
}

pub struct PoseidonBatchHash<F: PrimeField, P: PoseidonParameters<Fr = F>>{
    _field:      PhantomData<F>,
    _parameters: PhantomData<P>,
}

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

impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> PoseidonBatchHash<F, P> {

    // Function that does the scalar multiplication
    // It uses Montgomery multiplication
    // Constants are defined such that the result is x * t * 2^768,
    // that is the Montgomery representation of the operand x * t, and t is the 64-bit constant
    #[inline]
    fn scalar_mul (res: &mut F, state: &mut[F], mut start_idx_cst: usize) {

        state.iter_mut().for_each(|x| {
            let elem = x.mul(&P::MDS_CST[start_idx_cst]);
            start_idx_cst += 1;
            *res += &elem;
        });
    }

    // Function that does the mix matrix
    #[allow(dead_code)]
    #[inline]
    fn matrix_mix (state: &mut Vec<F>) {

        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![F::zero(); P::T];

        let mut idx_cst = 0;
        for i in 0..P::T {
            Self::scalar_mul(&mut new_state[i], state, idx_cst);
            idx_cst += P::T;
        }
        *state = new_state;

    }

    // Function that does the scalar multiplication
    // It uses a partial Montgomery multiplication defined as PM(x, t) = x * t * 2^-64 mod M
    // t is a 64-bit matrix constant. In the algorithm, the constants are represented in
    // partial Montgomery representation, i.e. t * 2^64 mod M
    #[inline]
    fn scalar_mul_fast (res: &mut F, state: &mut[F], mut start_idx_cst: usize) {
        state.iter_mut().for_each(|x| {
            let elem = P::MDS_CST_SHORT[start_idx_cst].mul_short(&x);
            start_idx_cst += 1;
            *res += &elem;
        });
    }

    // Function that does the mix matrix with fast algorithm
    #[inline]
    fn matrix_mix_short (state: &mut Vec<F>) {

        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![F::zero(); P::T];

        let mut idx_cst = 0;
        for i in 0..P::T {
            Self::scalar_mul_fast(&mut new_state[i], state, idx_cst);
            idx_cst += P::T;
        }
        *state = new_state;
    }

    fn poseidon_full_round(vec_state: &mut Vec<Vec<P::Fr>>, round_cst_idx: &mut usize) {

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
        // Use Montgomery simultaneous inversion
        let mut w: Vec<P::Fr> = Vec::new();
        let mut accum_prod = P::Fr::one();

        w.push(accum_prod);
        // Calculate the intermediate partial products
        for i in 0..vec_state.len() {
            for j in 0..P::T {
                accum_prod = accum_prod * &vec_state[i][j];
                w.push(accum_prod);
            }
        }

        // if the accum_prod is zero, it means that one of the S-Boxes is zero
        // in that case compute the inverses individually
        if accum_prod == P::Fr::zero() {
            for i in 0..vec_state.len() {
                for j in 0..P::T {
                    if vec_state[i][j] != P::Fr::zero() {
                        vec_state[i][j] = vec_state[i][j].inverse().unwrap();
                    }
                }
            }
        } else {

            // Calculate the inversion of the products
            // The inverse always exists in this case
            let mut w_bar = accum_prod.inverse().unwrap();

            // Extract the individual inversions
            let mut idx: i64 = w.len() as i64 - 2;
            for i in (0..vec_state.len()).rev() {
                for j in (0..P::T).rev() {
                    let vec_1 = vec_state[i][j].clone();
                    vec_state[i][j] = w_bar * &w[idx as usize];
                    w_bar = w_bar * &vec_1;
                    idx -= 1;
                }
            }

        }
    }

    fn poseidon_partial_round(vec_state: &mut Vec<Vec<P::Fr>>, round_cst_idx: &mut usize) {

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
        let mut w: Vec<P::Fr> = Vec::new();
        let mut accum_prod = P::Fr::one();

        w.push(accum_prod);
        // Calculate the intermediate partial products
        for i in 0..vec_state.len() {
            accum_prod = accum_prod * &vec_state[i][0];
            w.push(accum_prod);
        }

        // if the accum_prod is zero, it means that one of the S-Boxes is zero
        // in that case compute the inverses individually
        if accum_prod == P::Fr::zero() {
            for i in 0..(vec_state.len() - 1) {
                if vec_state[i][0] != P::Fr::zero() {
                    vec_state[i][0] = vec_state[i][0].inverse().unwrap();
                }
            }
        } else {

            // Calculate the inversion of the products
            // Use Montgomery simultaneous inversion
            let mut w_bar = accum_prod.inverse().unwrap();

            // Extract the individual inversions
            let mut idx: i64 = w.len() as i64 - 2;
            for i in (0..vec_state.len()).rev() {
                let vec_1 = vec_state[i][0].clone();
                vec_state[i][0] = w_bar * &w[idx as usize];
                w_bar = w_bar * &vec_1;
                idx -= 1;
            }
        }
    }

    fn poseidon_perm_gen(vec_state: &mut Vec<Vec<P::Fr>>) {

        // index that goes over the round constants
        let mut round_cst_idx: usize = 0;

        // Full rounds
        for _i in 0..P::R_F {
            Self::poseidon_full_round(vec_state, &mut round_cst_idx);

            // Perform the matrix mix
            for i in 0..vec_state.len() {
                //Self::matrix_mix_short(&mut vec_state[i]);
                Self::matrix_mix_short(&mut vec_state[i]);
            }

        }

        // Partial rounds
        for _i in 0..P::R_P {
            Self::poseidon_partial_round(vec_state, &mut round_cst_idx);

            // Perform the matrix mix
            for i in 0..vec_state.len() {
                //Self::matrix_mix_short(&mut vec_state[i]);
                Self::matrix_mix_short(&mut vec_state[i]);
            }
        }

        // Full rounds
        for _i in 0..(P::R_F - 1) {
            Self::poseidon_full_round(vec_state, &mut round_cst_idx);

            // Perform the matrix mix
            for i in 0..vec_state.len() {
                //Self::matrix_mix_short(&mut vec_state[i]);
                Self::matrix_mix_short(&mut vec_state[i]);
            }
        }

        Self::poseidon_full_round(vec_state, &mut round_cst_idx);
    }
}


impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> PoseidonHash<F, P> {

    // Function that does the scalar multiplication
    // It uses Montgomery multiplication
    // Constants are defined such that the result is x * t * 2^768,
    // that is the Montgomery representation of the operand x * t, and t is the 64-bit constant
    #[inline]
    fn scalar_mul (res: &mut F, state: &mut[F], start_idx_cst: usize) {
        state.iter_mut().enumerate().for_each(|(i,x)| {
            let elem = x.mul(&P::MDS_CST[start_idx_cst + i]);
            *res += &elem;
        });
    }

    // Function that does the mix matrix
    #[inline]
    fn matrix_mix (state: &mut Vec<F>) {

        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![F::zero(); P::T];

        let mut idx_cst = 0;
        for i in 0..P::T {
            Self::scalar_mul(&mut new_state[i], state, idx_cst);
            idx_cst += P::T;
        }
        *state = new_state;

    }

    // Function that does the scalar multiplication with fast algorithm
    // It uses a partial Montgomery multiplication defined as PM(x, t) = x * t * 2^-64 mod M
    // t is a 64-bit matrix constant. In the algorithm, the constants are represented in
    // partial Montgomery representation, i.e. t * 2^64 mod M
    #[inline]
    fn scalar_mul_fast (res: &mut F, state: &mut[F], mut start_idx_cst: usize) {
        state.iter_mut().for_each(|x| {
            let elem = P::MDS_CST_SHORT[start_idx_cst].mul_short(&x);
            start_idx_cst += 1;
            *res += &elem;
        });
    }

    // Function that does the mix matrix with fast algorithm
    #[inline]
    fn matrix_mix_short (state: &mut Vec<F>) {

        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![F::zero(); P::T];

        let mut idx_cst = 0;
        for i in 0..P::T {
            Self::scalar_mul_fast(&mut new_state[i], state, idx_cst);
            idx_cst += P::T;
        }
        *state = new_state;
    }

    fn poseidon_perm (state: &mut Vec<F>) {

        let use_fast = true;

        // index that goes over the round constants
        let mut round_cst_idx = 0;

        // First full rounds
        for _i in 0..P::R_F {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            // Optimization for the inversion S-Box
            let w2 = state[0] * &state[1];
            let w = state[2] * &w2;
            if w == P::Fr::zero() {
                // At least one of the S-Boxes is zero
                for d in state.iter_mut() {
                    // The S-BOX is an inversion function
                    if *d != P::Fr::zero() {
                        *d = (*d).inverse().unwrap();
                    }
                }
            } else {
                let mut w_bar = w.inverse().unwrap();

                let z_2 = w_bar * &w2;
                w_bar = w_bar * &state[2];
                state[2] = z_2;
                let z_1 = w_bar * &state[0];
                state[0] = w_bar * &state[1];
                state[1] = z_1;
            }

            // Perform the matrix mix
            if use_fast {
                Self::matrix_mix_short(state);
            } else {
                Self::matrix_mix(state);
            }

        }

        // Partial rounds
        for _i in 0..P::R_P {

            // Add the round constants to the state vector
            for d in state.iter_mut() {
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply S-BOX only to the first element of the state vector
            if state[0]!=P::Fr::zero() {
                state[0] = state[0].inverse().unwrap();
            }

            // Apply the matrix mix
            if use_fast {
                Self::matrix_mix_short(state);
            } else {
                Self::matrix_mix(state);
            }
        }

        // Second full rounds
        // Process only to R_F -1 iterations. The last iteration does not contain a matrix mix
        for _i in 0..(P::R_F-1) {

            // Add the round constants
            for d in state.iter_mut() {
                //let rc = MNT4753Fr::from_str(ROUND_CST[round_cst_idx]).map_err(|_| ()).unwrap();
                let rc = P::ROUND_CST[round_cst_idx];
                *d += &rc;
                round_cst_idx += 1;
            }

            // Apply the S-BOX to each of the elements of the state vector
            // Optimization for the inversion S-Box
            let w2 = state[0] * &state[1];
            let w = state[2] * &w2;
            if w == P::Fr::zero() {
                // At least one of the S-Boxes is zero
                for d in state.iter_mut() {
                    // The S-BOX is an inversion function
                    if *d != P::Fr::zero() {
                        *d = (*d).inverse().unwrap();
                    }
                }
            } else {
                let mut w_bar = w.inverse().unwrap();

                let z_2 = w_bar * &w2;
                w_bar = w_bar * &state[2];
                state[2] = z_2;
                let z_1 = w_bar * &state[0];
                state[0] = w_bar * &state[1];
                state[1] = z_1;
            }

            // Apply matrix mix
            if use_fast {
                Self::matrix_mix_short(state);
            } else {
                Self::matrix_mix(state);
            }
        }

        // Last full round does not perform the matrix_mix
        // Add the round constants
        for d in state.iter_mut() {
            let rc = P::ROUND_CST[round_cst_idx];
            *d += &rc;
            round_cst_idx += 1;
        }

        // Apply the S-BOX to each of the elements of the state vector
        // Optimization for the inversion S-Box
        let w2 = state[0] * &state[1];
        let w = state[2] * &w2;
        if w == P::Fr::zero() {
            for d in state.iter_mut() {
                // The S-BOX is an inversion function
                if *d != P::Fr::zero() {
                    *d = (*d).inverse().unwrap();
                }
            }
        } else {
            let mut w_bar = w.inverse().unwrap();

            let z_2 = w_bar * &w2;
            w_bar = w_bar * &state[2];
            state[2] = z_2;
            let z_1 = w_bar * &state[0];
            state[0] = w_bar * &state[1];
            state[1] = z_1;
        }
    }
}


impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr = F>> FieldBasedHash for PoseidonHash<F, P> {
    type Data = F;
    type Parameters = P;

    fn evaluate(input: &[F]) -> Result<F, Error> {

        // state is a vector of 3 elements. They are initialized to constants that are obtained after applying a permutation to a zero elements vector
        let mut state = vec![P::AFTER_ZERO_PERM[0], P::AFTER_ZERO_PERM[1], P::AFTER_ZERO_PERM[2]];

        // calculate the number of cycles to process the input dividing in portions of rate elements
        let num_cycles = input.len() / P::R;
        // check if the input is a multiple of the rate by calculating the remainder of the division
        let rem = input.len() % P::R;

        // index to process the input
        let mut input_idx = 0;
        // iterate of the portions of rate elements
        for _i in 0..num_cycles {
            // add the elements to the state vector. Add rate elements
            for j in 0..P::R {
                state[j] += &input[input_idx];
                input_idx += 1;
            }
            // for application to a 2-1 Merkle tree, add the constant 3 to the third state vector
            state[P::R] += &P::C2;

            // apply permutation after adding the input vector
            Self::poseidon_perm(&mut state);
        }

        // in case the input is not a multiple of the rate process the remainder part padding a zero
        if rem != 0 {
            state[0] += &input[input_idx];
            state[P::R] += &P::C2;
            // apply permutation after adding the input vector
            Self::poseidon_perm(&mut state);
        }

        // return the first element of the state vector as the hash digest
        Ok(state[0])
    }
}

impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr = F>> BatchFieldBasedHash for PoseidonBatchHash<F, P> {
    type Data = F;
    type Parameters = P;

    fn batch_evaluate_2_1(input_array: &mut[F], output_array: &mut[F]) {

        // Input:
        // This function calculates the hashes of pairs of inputs.
        // The inputs are arranged in an array and arranged as pairs
        // Example:
        // (d_00, d01, d_10, d_11, d_20, d_21, ...
        // Output:
        // The output will be placed in the same array taking half of the positions
        // as the rate of the hash function is 2 field elements

        // Checks that input contains data
        assert_ne!(input_array.len(), 0, "Input to the hash has length 0.");
        assert_eq!(input_array.len() % 2, 0, "The length of the input to the hash is not even.");

        let input_length = input_array.len() / 2;

        // Assign pre-computed values of the state vector equivalent to a permutation with zero element state vector
        let state_z = vec![P::AFTER_ZERO_PERM[0], P::AFTER_ZERO_PERM[1], P::AFTER_ZERO_PERM[2]];

        // Copy the result of the permutation to a vector of state vectors of the length equal to the length of the input
        // state is a vector of 3-element state vector.
        let mut state = Vec::new();
        for _i in 0..input_length {
            state.push(state_z.clone());
        }

        // input_idx is to scan the input_array
        let mut input_idx = 0;

        for k in 0..input_length {
            state[k][0] += &input_array[input_idx];
            input_idx += 1;
            state[k][1] += &input_array[input_idx];
            input_idx += 1;
            // constant to add for a 2-1 Merkle tree
            state[k][2] += &P::C2;
        }

        // apply permutation after adding the input vector
        Self::poseidon_perm_gen(&mut state);

        // overwrite the input with the result of the hash
        for k in 0..input_array.len()/2 {
            output_array[k] = state[k][0];
        }
    }

    fn merkle_tree(input_vec: &mut[Self::Data], output_vec: &mut[Self::Data], input_size: usize){
        // Supporting function that processes the inputs and outputs in chunks

        let num_cores = num_cpus::get();

        assert_eq!(input_vec.len() % 2, 0, "The length of the input to the hash is not even.");
        assert_eq!(output_vec.len() >= (input_vec.len() / 2), true, "The length of the output is not long enough.");

        if input_size < 2 * num_cores {
            input_vec.par_chunks_mut(2).zip(output_vec.par_chunks_mut(1)).for_each( |(p1,p2)| {
                Self::batch_evaluate_2_1(p1,p2);
            });
            return;
        }

        use rayon::prelude::*;

        input_vec.par_chunks_mut(input_size/num_cores).zip(output_vec.par_chunks_mut(input_size/num_cores/2)).for_each( |(p1, p2)| {
            Self::batch_evaluate_2_1(p1,p2);
        });


    }

    fn merkle_tree_2_1(array: &mut [Self::Data], input_size: usize) {
        // Main function for the Merkle-tree computation

        let mut copy_vec = &mut array[..];
        let mut size_input = input_size;

        while size_input > 1{
            let (input_vec, output_vec) = copy_vec.split_at_mut(size_input);
            Self::merkle_tree(input_vec, output_vec, size_input);
            copy_vec = output_vec;
            size_input = size_input / 2;
        }
    }

}

pub type MNT4PoseidonHash = PoseidonHash<MNT4753Fr, MNT4753PoseidonParameters>;
pub type MNT6PoseidonHash = PoseidonHash<MNT6753Fr, MNT6753PoseidonParameters>;

#[cfg(test)]
mod test {
    use super::*;
    use rayon::prelude::*;
    use rand_xorshift::XorShiftRng;
    use std::str::FromStr;
    use crate::{FieldBasedHash, BatchFieldBasedHash, PoseidonBatchHash};
    use super::rand::SeedableRng;
    use algebra::{UniformRand, Field};
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

        // Fill the upper part of the vec with zeros. This will be the output
        array_elem_4753.resize(2 * num_rounds + num_rounds, MNT4753Fr::zero());

        // =============================================================================
        // Calculate Poseidon Hash for mnt4753
        let now_4753 = Instant::now();

        let mut output_4753 = Vec::new();

        vec_vec_elem_4753.iter_mut().for_each(|p| {
            let output = Mnt4PoseidonHash::evaluate(p);
            output_4753.push(output.unwrap());
        });

        let new_now_4753 = Instant::now();

        // Calculate Poseidon Hash for mnt4753 batch evaluation

        let (input_vec, output_vec) = array_elem_4753.split_at_mut(2 * num_rounds);


        let now_4753_batch = Instant::now();
        input_vec
            .par_chunks_mut(num_rounds/2)
            .zip(output_vec
                .par_chunks_mut(num_rounds/2/2))
            .for_each(|(p1,p2)| {
                Mnt4BatchPoseidonHash::batch_evaluate_2_1(p1, p2);
        });
        let new_now_4753_batch = Instant::now();

        // =============================================================================
        // Compare results
        for i in 0..num_rounds {
            if output_4753[i] != output_vec[i] {
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
        //  the number of rounds to test
        let num_rounds = 1000;

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

        // Fill the upper part of the vec with zeros. This will be the output
        array_elem_6753.resize(2 * num_rounds + num_rounds, MNT6753Fr::zero());

        // =============================================================================
        // Calculate Poseidon Hash for mnt4753
        let now_6753 = Instant::now();

        let mut output_6753 = Vec::new();

        vec_vec_elem_6753.iter_mut().for_each(|p| {
            let output = Mnt6PoseidonHash::evaluate(p);
            output_6753.push(output.unwrap());
        });

        let new_now_6753 = Instant::now();

        // Calculate Poseidon Hash for mnt4753 batch evaluation

        let (input_vec, output_vec) = array_elem_6753.split_at_mut(2 * num_rounds);


        let now_6753_batch = Instant::now();
        input_vec
            .par_chunks_mut(num_rounds/2)
            .zip(output_vec
                .par_chunks_mut(num_rounds/2/2))
            .for_each(|(p1,p2)| {
                Mnt6BatchPoseidonHash::batch_evaluate_2_1(p1, p2);
            });
        let new_now_6753_batch = Instant::now();

        // =============================================================================
        // Compare results
        for i in 0..num_rounds {
            if output_6753[i] != output_vec[i] {
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

    #[test]
    fn test_merkle_tree() {

        // =============================================================================
        // Computation for MNT4
        type Mnt4BatchPoseidonHash = PoseidonBatchHash<MNT4753Fr, MNT4753PoseidonParameters>;
        //  the size of the input
        let input_size = 1024*1024;

        // the vectors that store random input data
        let mut vec_elem_4753 = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // fill input vector with random elements
        vec_elem_4753.resize_with(input_size, || {MNT4753Fr::rand(&mut rng)});

        let mut sz1 = input_size/2;
        let mut sz2 = input_size + sz1;
        while sz1!=0 {
            vec_elem_4753.resize(sz2, MNT4753Fr::zero());
            sz1 /= 2;
            sz2 += sz1;
        }

        // =============================================================================
        // Calculate Merkle tree for mnt4753 with Batch Poseidon Hash
        let now_4753_batch = Instant::now();

        Mnt4BatchPoseidonHash::merkle_tree_2_1(&mut vec_elem_4753, input_size);

        let new_now_4753_batch = Instant::now();

        println!("Merkle root MNT4753 = {:?}", vec_elem_4753.last());

        // =============================================================================
        // Report the timing results

        let duration_4753_batch = new_now_4753_batch.duration_since(now_4753_batch);
        println!("Time for input_size={}, MNT4753 batch Merkle tree= {:?}", input_size, duration_4753_batch.as_millis());

        // =============================================================================
        // Computation for MNT6
        type Mnt6BatchPoseidonHash = PoseidonBatchHash<MNT6753Fr, MNT6753PoseidonParameters>;
        //  the size of the input
        let input_size = 1024*1024;

        // the vectors that store random input data
        let mut vec_elem_6753 = Vec::new();

        // the random number generator to generate random input data
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        // fill input vector with random elements
        vec_elem_6753.resize_with(input_size, || {MNT6753Fr::rand(&mut rng)});

        let mut sz1 = input_size/2;
        let mut sz2 = input_size + sz1;
        while sz1!=0 {
            vec_elem_6753.resize(sz2, MNT6753Fr::zero());
            sz1 /= 2;
            sz2 += sz1;
        }

        // =============================================================================
        // Calculate Merkle tree for mnt6753 with Batch Poseidon Hash
        let now_6753_batch = Instant::now();

        Mnt6BatchPoseidonHash::merkle_tree_2_1(&mut vec_elem_6753, input_size);

        let new_now_6753_batch = Instant::now();

        println!("Merkle root MNT6753 = {:?}", vec_elem_6753.last());

        // =============================================================================
        // Report the timing results

        let duration_6753_batch = new_now_6753_batch.duration_since(now_6753_batch);
        println!("Time for input_size={}, MNT6753 batch Merkle tree= {:?}", input_size, duration_6753_batch.as_millis());

    }

}