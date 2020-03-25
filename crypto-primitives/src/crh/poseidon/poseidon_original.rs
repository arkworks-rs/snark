use crate::crh::poseidon::PoseidonParameters;

use algebra::{PrimeField, MulShort};
//use algebra::arithmetic::{mac_with_carry, adc};
//use algebra::fields::FpParameters;
//use algebra::biginteger::BigInteger768;
//use algebra::BigInteger;

use std::marker::PhantomData;
use crate::crh::{FieldBasedHash, BatchFieldBasedHash};
use crate::Error;

pub struct PoseidonHash<F: PrimeField, P: PoseidonParameters<Fr = F>>{
    _field:      PhantomData<F>,
    _parameters: PhantomData<P>,
}

pub struct PoseidonBatchHash<F: PrimeField, P: PoseidonParameters<Fr = F>>{
    _field:      PhantomData<F>,
    _parameters: PhantomData<P>,
}

impl<F: PrimeField + MulShort, P: PoseidonParameters<Fr=F>> PoseidonBatchHash<F, P> {

    // Function that does the mix matrix
    // It uses Montgomery multiplication
    // Constants are defined such that the result is x * t * 2^768,
    // that is the Montgomery representation of the operand x * t, and t is the 64-bit constant
    fn matrix_mix (state: &mut Vec<F>) {

        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![F::zero(); P::T];

        let m_11 = P::MDS_CST[0];
        let m_12 = P::MDS_CST[1];
        let m_13 = P::MDS_CST[2];

        // scalar multiplication for position 0 of the state vector
        let elem_0 = state[0].mul(&m_11);
        let elem_1 = state[1].mul(&m_12);
        let elem_2 = state[2].mul(&m_13);

        new_state[0] = elem_0;
        new_state[0] += &elem_1;
        new_state[0] += &elem_2;

        // scalar multiplication for position 1 of the state vector
        let m_21 = P::MDS_CST[3];
        let m_22 = P::MDS_CST[4];
        let m_23 = P::MDS_CST[5];

        let elem_3 = state[0].mul(&m_21);
        let elem_4 = state[1].mul(&m_22);
        let elem_5 = state[2].mul(&m_23);

        new_state[1] = elem_3;
        new_state[1] += &elem_4;
        new_state[1] += &elem_5;

        // scalar multiplication for the position 2 of the state vector
        let m_31 = P::MDS_CST[6];
        let m_32 = P::MDS_CST[7];
        let m_33 = P::MDS_CST[8];

        let elem_6 = state[0].mul(&m_31);
        let elem_7 = state[1].mul(&m_32);
        let elem_8 = state[2].mul(&m_33);

        new_state[2] = elem_6;
        new_state[2] += &elem_7;
        new_state[2] += &elem_8;

        // copy the result to the state vector
        state[0] = new_state[0];
        state[1] = new_state[1];
        state[2] = new_state[2];

    }

    // Function that does the mix matrix with fast algorithm
    // It uses a partial Montgomery multiplication defined as PM(x, t) = x * t * 2^-64 mod M
    // t is a 64-bit matrix constant. In the algorithm, the constants are represented in
    // partial Montgomery representation, i.e. t * 2^64 mod M
    fn matrix_mix_short (state: &mut Vec<F>) {

        //use algebra::MulShort;

        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![F::zero(); P::T];

        let m_11 = P::MDS_CST_SHORT[0];
        let m_12 = P::MDS_CST_SHORT[1];
        let m_13 = P::MDS_CST_SHORT[2];

        let elem_0 = m_11.mul_short(&state[0]);
        let elem_1 = m_12.mul_short(&state[1]);
        let elem_2 = m_13.mul_short(&state[2]);

        new_state[0] = elem_0;
        new_state[0] += &elem_1;
        new_state[0] += &elem_2;

        // scalar multiplication for position 1 of the state vector
        let m_21 = P::MDS_CST_SHORT[3];
        let m_22 = P::MDS_CST_SHORT[4];
        let m_23 = P::MDS_CST_SHORT[5];
        // Because the matrix is symmetric, they can be replaced by
        // the following instructions.
        // let m_21 = m_12;
        // let m_22 = m_11;
        // let m_23 = P::MDS_CST_SHORT[5];


        let elem_3 = m_21.mul_short(&state[0]);
        let elem_4 = m_22.mul_short(&state[1]);
        let elem_5 = m_23.mul_short(&state[2]);

        new_state[1] = elem_3;
        new_state[1] += &elem_4;
        new_state[1] += &elem_5;

        // scalar multiplication for the position 2 of the state vector
        let m_31 = P::MDS_CST_SHORT[6];
        let m_32 = P::MDS_CST_SHORT[7];
        let m_33 = P::MDS_CST_SHORT[8];
        // Because the matrix is symmetric, they can be replaced by
        // the following instructions
        // let m_31 = m_13;
        // let m_32 = m_23;
        // let m_33 = m_11;

        let elem_6 = m_31.mul_short(&state[0]);
        let elem_7 = m_32.mul_short(&state[1]);
        let elem_8 = m_33.mul_short(&state[2]);

        new_state[2] = elem_6;
        new_state[2] += &elem_7;
        new_state[2] += &elem_8;

        // copy the result to the state vector
        state[0] = new_state[0];
        state[1] = new_state[1];
        state[2] = new_state[2];
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
        // vec_state = [ [z_1, z_2, z_3] [z_4, z_5, z_6] ]

        // w_0 = 1
        // w_1 = 1 * z_1
        // w_2 = 1 * z_1 * z_2
        // w_3 = 1 * z_1 * z_2 * z_3
        // w_4 = 1 * z_1 * z_2 * z_3 * z_4
        // w_5 = 1 * z_1 * z_2 * z_3 * z_4 * z_5
        // w_6 = 1 * z_1 * z_2 * z_3 * z_4 * z_5 * z_6

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

            // w_bar = (1 * z_1 * z_2 * z_3 * z_4 * z_5 * z_6) ^ -1
            // w_bar = z_1^-1 * z_2^-1 * z_3^-1 * z_4^-1 * z_5^-1 * z_6^-1

            // z_6 = w_bar * w_5 = z_6^-1
            // w_bar <= w_bar * z_6 = z_1^-1 * z_2^-1 * z_3^-1 * z_4^-1 * z_5^-1

            // z_5 = w_bar * w_4 = z_5^-1
            // w_bar <= w_bar * z_5 = z_1^-1 * z_2^-1 * z_3^-1 * z_4^-1

            // z_4 = w_bar * w_3 = z_4^-1
            // w_bar <= w_bar * z_4 = z_1^-1 * z_2^-1 * z_3^-1

            // z_3 = w_bar * w_2 = z_3^-1
            // w_bar <= w_bar * z_3 = z_1^-1 * z_2^-1

            // z_2 = w_bar * w_1 = z_2^-1
            // w_bar <= w_bar * z_2 = z_1^-1

            // z_1 = w_bar * w_0 = z_1^-1
            // w_bar <= w_bar * z_1 = 1

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

    // Function that does the mix matrix
    // It uses Montgomery multiplication
    // Constants are defined such that the result is x * t * 2^768,
    // that is the Montgomery representation of the operand x * t, and t is the 64-bit constant
    fn matrix_mix (state: &mut Vec<F>) {

        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![F::zero(); P::T];

        let m_11 = P::MDS_CST[0];
        let m_12 = P::MDS_CST[1];
        let m_13 = P::MDS_CST[2];

        // scalar multiplication for position 0 of the state vector
        let elem_0 = state[0].mul(&m_11);
        let elem_1 = state[1].mul(&m_12);
        let elem_2 = state[2].mul(&m_13);

        new_state[0] = elem_0;
        new_state[0] += &elem_1;
        new_state[0] += &elem_2;

        // scalar multiplication for position 1 of the state vector
        let m_21 = P::MDS_CST[3];
        let m_22 = P::MDS_CST[4];
        let m_23 = P::MDS_CST[5];

        let elem_3 = state[0].mul(&m_21);
        let elem_4 = state[1].mul(&m_22);
        let elem_5 = state[2].mul(&m_23);

        new_state[1] = elem_3;
        new_state[1] += &elem_4;
        new_state[1] += &elem_5;

        // scalar multiplication for the position 2 of the state vector
        let m_31 = P::MDS_CST[6];
        let m_32 = P::MDS_CST[7];
        let m_33 = P::MDS_CST[8];

        let elem_6 = state[0].mul(&m_31);
        let elem_7 = state[1].mul(&m_32);
        let elem_8 = state[2].mul(&m_33);

        new_state[2] = elem_6;
        new_state[2] += &elem_7;
        new_state[2] += &elem_8;

        // copy the result to the state vector
        state[0] = new_state[0];
        state[1] = new_state[1];
        state[2] = new_state[2];

    }

    // Function that does the mix matrix with fast algorithm
    // It uses a partial Montgomery multiplication defined as PM(x, t) = x * t * 2^-64 mod M
    // t is a 64-bit matrix constant. In the algorithm, the constants are represented in
    // partial Montgomery representation, i.e. t * 2^64 mod M
    fn matrix_mix_short (state: &mut Vec<F>) {

        //use algebra::MulShort;

        // the new state where the result will be stored initialized to zero elements
        let mut new_state = vec![F::zero(); P::T];

        let m_11 = P::MDS_CST_SHORT[0];
        let m_12 = P::MDS_CST_SHORT[1];
        let m_13 = P::MDS_CST_SHORT[2];

        let elem_0 = m_11.mul_short(&state[0]);
        let elem_1 = m_12.mul_short(&state[1]);
        let elem_2 = m_13.mul_short(&state[2]);

        new_state[0] = elem_0;
        new_state[0] += &elem_1;
        new_state[0] += &elem_2;

        // scalar multiplication for position 1 of the state vector
        let m_21 = P::MDS_CST_SHORT[3];
        let m_22 = P::MDS_CST_SHORT[4];
        let m_23 = P::MDS_CST_SHORT[5];
        // Because the matrix is symmetric, they can be replaced by
        // the following instructions:
        // let m_21 = m_12;
        // let m_22 = m_11;
        // let m_23 = P::MDS_CST_SHORT[5];


        let elem_3 = m_21.mul_short(&state[0]);
        let elem_4 = m_22.mul_short(&state[1]);
        let elem_5 = m_23.mul_short(&state[2]);

        new_state[1] = elem_3;
        new_state[1] += &elem_4;
        new_state[1] += &elem_5;

        // scalar multiplication for the position 2 of the state vector
        let m_31 = P::MDS_CST_SHORT[6];
        let m_32 = P::MDS_CST_SHORT[7];
        let m_33 = P::MDS_CST_SHORT[8];
        // Because the matrix is symmetric, they can be replaced by
        // the following instructions:
        // let m_31 = m_13;
        // let m_32 = m_23;
        // let m_33 = m_11;

        let elem_6 = m_31.mul_short(&state[0]);
        let elem_7 = m_32.mul_short(&state[1]);
        let elem_8 = m_33.mul_short(&state[2]);

        new_state[2] = elem_6;
        new_state[2] += &elem_7;
        new_state[2] += &elem_8;

        // copy the result to the state vector
        state[0] = new_state[0];
        state[1] = new_state[1];
        state[2] = new_state[2];
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
            // Assuming state vector of 3 elements
            // let bc = state[1]*&state[2];
            // let abc = state[0]*&bc;
            // if abc != P::Fr::zero() {
            //     let abc_inv = abc.inverse().unwrap();
            //
            //     let a_inv = abc_inv*&bc;
            //     let b_inv = abc_inv*&state[0]*&state[2];
            //     let c_inv = abc_inv*&state[0]*&state[1];
            //
            //     state[0] = a_inv;
            //     state[1] = b_inv;
            //     state[2] = c_inv;
            //
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
            // Assuming state vector of 3 elements

            // let bc = state[1]*&state[2];
            // let abc = state[0]*&bc;
            // if abc != P::Fr::zero() {
            //     let abc_inv = abc.inverse().unwrap();
            //
            //     let a_inv = abc_inv*&bc;
            //     let b_inv = abc_inv*&state[0]*&state[2];
            //     let c_inv = abc_inv*&state[0]*&state[1];
            //
            //     state[0] = a_inv;
            //     state[1] = b_inv;
            //     state[2] = c_inv;
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
        // Assuming state vector of 3 elements
        // let bc = state[1]*&state[2];
        // let abc = state[0]*&bc;
        // if abc != P::Fr::zero() {
        //     let abc_inv = abc.inverse().unwrap();
        //
        //     let a_inv = abc_inv*&bc;
        //     let b_inv = abc_inv*&state[0]*&state[2];
        //     let c_inv = abc_inv*&state[0]*&state[1];
        //
        //     state[0] = a_inv;
        //     state[1] = b_inv;
        //     state[2] = c_inv;

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

    fn batch_evaluate_2_1(input_array: &mut[F]) {

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
            input_array[k] = state[k][0];
        }
    }
}
