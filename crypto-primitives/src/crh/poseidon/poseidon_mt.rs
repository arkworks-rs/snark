use crate::crh::poseidon::PoseidonParameters;

use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::{Field, PrimeField, SquareRootField, UniformRand};
use std::ops::Mul;

use std::time::Instant;

use algebra::biginteger::BigInteger768;
use algebra::{to_bytes, ToBytes};
use algebra::field_new;

// Function that does the mix matrix
fn matrix_mix<P:PoseidonParameters> (state: &mut Vec<P::Fr>) {

    // the new state where the result will be stored initialized to zero elements
    let mut new_state = vec![P::Fr::zero(); P::T];

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

fn poseidon_full_round<T:PoseidonParameters> (vec_state: &mut Vec<Vec<T::Fr>>, round_cst_idx: &mut usize) {

    // Full rounds
    for _i in 0..T::R_F {

        // For each of the state vector element position
        for j in 0..T::T {

            // get the constant associated to state vector element position
            let rc = T::ROUND_CST[*round_cst_idx];

            // go over each state vector
            for k in 0..vec_state.len()-1 {
                vec_state[k][j] += &rc;
            }
            *round_cst_idx += 1;
        }
    }

    // Apply the S-BOX to each of the elements of the state vector
    let mut partial_prod: Vec<T::Fr> = Vec::new();
    //let mut accum_prod = T::Fr::from_repr(vec_state[vec_state.len()-1][T::T-1]);
    let mut accum_prod = T::Fr::one();

    partial_prod.push(accum_prod);
    // Calculate the intermediate partial products
    for i in (0..(vec_state.len() - 1)).rev() {
        for j in (0..(T::T - 1)).rev() {
            accum_prod = accum_prod * &vec_state[i][j];
            partial_prod.push(accum_prod);
        }
    }

    // partial_prod =
    // 0 => 1
    // 2 => x_6
    // 3 => x_6 * x_5
    // 4 => x_6 * x_5 * x_4
    // 5 => x_6 * x_5 * x_4 * x_3
    // 6 => x_6 * x_5 * x_4 * x_3 * x_2
    // 7 => x_6 * x_5 * x_4 * x_3 * x_2 * x_1

    // Calculate the inversion of the products
    let inv_prod = accum_prod.inverse();
    let mut v = T::Fr::one();
    match inv_prod {
        None => println!("Field inversion error of the accumulated products"),
        Some(inv) => {
            v = inv;
        }
    }

    //println!("partial_prod = {:?}", partial_prod);

    // Extract the individual inversions
    let mut idx = partial_prod.len()-2;
    let mut w = T::Fr::one();
    for i in 0..(vec_state.len() - 1) {
        for j in 0..(T::T - 1) {
            let vec_1 = vec_state[i][j].clone();
            vec_state[i][j] = v * &w * &partial_prod[idx];
            w = w * &vec_1;
            idx -= 1;
        }
    }
    // // Perform the matrix mix
    // for i in 0..(vec_state.len() - 1) {
    //     matrix_mix::<T>(&mut vec_state[i]);
    // }
}

fn poseidon_partial_round<T:PoseidonParameters> (vec_state: &mut Vec<Vec<T::Fr>>, round_cst_idx: &mut usize) {

    // Partial rounds
    for _i in 0..T::R_P {

        // For each of the state vector element position
        for j in 0..T::T {

            // get the constant associated to state vector element position
            let rc = T::ROUND_CST[*round_cst_idx];

            // go over each state vector
            for k in 0..vec_state.len()-1 {
                vec_state[k][j] += &rc;
            }
            *round_cst_idx += 1;
        }
    }

    // Apply the S-BOX to the first elements of each of the state vector
    let mut partial_prod: Vec<T::Fr> = Vec::new();
    let mut accum_prod = T::Fr::one();

    partial_prod.push(accum_prod);
    // Calculate the intermediate partial products
    for i in (0..(vec_state.len() - 1)).rev() {
        accum_prod = accum_prod * &vec_state[i][0];
        partial_prod.push(accum_prod);
    }

    // partial_prod =
    // 0 => 1
    // 2 => x_6
    // 3 => x_6 * x_5
    // 4 => x_6 * x_5 * x_4
    // 5 => x_6 * x_5 * x_4 * x_3
    // 6 => x_6 * x_5 * x_4 * x_3 * x_2
    // 7 => x_6 * x_5 * x_4 * x_3 * x_2 * x_1

    // Calculate the inversion of the products
    let inv_prod = accum_prod.inverse();
    let mut v = T::Fr::one();
    match inv_prod {
        None => println!("Field inversion error of the accumulated products"),
        Some(inv) => {
            v = inv;
        }
    }

    //println!("partial_prod = {:?}", partial_prod);

    // Extract the individual inversions
    let mut idx = partial_prod.len()-2;
    let mut w = T::Fr::one();
    for i in 0..(vec_state.len() - 1) {
        let vec_1 = vec_state[i][0].clone();
        vec_state[i][0] = v * &w * &partial_prod[idx];
        w = w * &vec_1;
        idx -= 1;
    }

    // // Perform the matrix mix
    // for i in 0..(vec_state.len() - 1) {
    //     matrix_mix::<T>(&mut vec_state[i]);
    // }

}

fn poseidon_perm_gen<T:PoseidonParameters> (vec_state: &mut Vec<Vec<T::Fr>>) {

    // index that goes over the round constants
    let mut round_cst_idx: usize = 0;

    poseidon_full_round::<T> (vec_state, &mut round_cst_idx);

    // Perform the matrix mix
    for i in 0..(vec_state.len() - 1) {
        matrix_mix::<T>(&mut vec_state[i]);
    }

    poseidon_partial_round::<T> (vec_state, &mut round_cst_idx);

    // Perform the matrix mix
    for i in 0..(vec_state.len() - 1) {
        matrix_mix::<T>(&mut vec_state[i]);
    }

    poseidon_full_round::<T> (vec_state, &mut round_cst_idx);

}

pub fn poseidon_engine_gen<T: PoseidonParameters>(input: &mut Vec<Vec<T::Fr>>) -> Vec<T::Fr> {

    // input:
    // [0] : (d_00, d_01, d_02, ...)
    // [1] : (d_10, d_11, d_12, ...)
    // [2] : (d_20, d_21, d_22, ...)
    // ...
    // It is assumed that the inputs have the same length

    // Checks that input contains data
    assert_ne!(input.len(),0, "Input to the hash has length 0.");

    // Checks that the inputs contain vectors of the same length
    let length_0 = input[0].len();
    for x in input.iter_mut() {
        assert_eq!((*x).len(),length_0, "Input vectors to hash do not have the same length.");
    }

    // calculate the number of cycles to process the input dividing in portions of rate elements
    // we assume the inputs lengths are the same for input1 and input2
    let num_cycles = length_0 / T::R;
    // check if the input is a multiple of the rate by calculating the remainder of the division
    let rem = length_0 % T::R;

    // First I initialized with a single state vector of zero and call the Poseidon hash and then
    // copy the result of the permutation to a vector of state vectors with the same length as the input
    let mut state_zero_vec = Vec::new();
    let state_z = vec![T::ZERO, T::ZERO, T::ZERO];
    state_zero_vec.push(state_z);

    // apply permutation to all zeros state vector
    poseidon_perm_gen::<T>(&mut state_zero_vec);

    // Copy the result of the permutation to a vector of state vectors of the length equal to the length of the input
    // state is a vector of 3-element state vector.
    let mut state = Vec::new();
    for i in 0..input.len() {
        let mut d = state_zero_vec[0].clone();
        state.push(d);
    }

    // input_idx is an index to process the inputs
    let mut input_idx = 0;

    // iterate of the portions of rate elements
    for _i in 0..num_cycles {
        // add the elements to the state vector. Add rate elements
        for j in 0..T::R {
            for k in 0..input.len() {
                state[k][j] += &input[k][input_idx];
            }
            input_idx += 1;
        }
        // for application to a 2-1 Merkle tree, add the constant 3 to the third element of the state vector
        for k in 0..input.len() {
            state[k][T::R] += &T::C2;
        }
        // apply permutation after adding the input vector
        poseidon_perm_gen::<T>(&mut state);
    }

    // in case the input is not a multiple of the rate process the remainder part padding a zero
    if rem != 0 {
        for k in 0..input.len() {
            state[k][0] += &input[k][input_idx];
            state[k][T::R] += &T::C2;
        }
        // apply permutation after adding the input vector
        poseidon_perm_gen::<T>(&mut state);
    }

    // the hashes of the inputs are the first elements of the state vectors
    // output is a vector of the hashes
    let mut output:Vec<T::Fr> = Vec::new();
    for k in 0..input.len() {
        output.push(state[k][0]);
    }

    // return output
    output
}