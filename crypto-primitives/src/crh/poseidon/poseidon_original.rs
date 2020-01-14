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

fn poseidon_perm<T:PoseidonParameters> (state: &mut Vec<T::Fr>) {

    // index that goes over the round constants
    let mut round_cst_idx = 0;

    // First full rounds
    for _i in 0..T::R_F {

        // Add the round constants to the state vector
        for d in state.iter_mut() {
            let rc = T::ROUND_CST[round_cst_idx];
            *d += &rc;
            round_cst_idx += 1;
        }

        // Apply the S-BOX to each of the elements of the state vector
        for d in state.iter_mut() {

            // The S-BOX is an inversion function
            let elem_state = (*d).inverse();
            match elem_state {
                None => println!("Field inversion error"),
                Some(inv) => {
                    *d = inv;
                },
            }
        }

        // Perform the matrix mix
        matrix_mix::<T> (state);
    }

    // Partial rounds
    for _i in 0..T::R_P {

        // Add the round constants to the state vector
        for d in state.iter_mut() {
            let rc = T::ROUND_CST[round_cst_idx];
            *d += &rc;
            round_cst_idx += 1;
        }

        // Apply S-BOX only to the first element of the state vector
        let t1 = (state[0]).inverse();
        match t1 {
            None => println!("Field inversion error"),
            Some(inv) => {
                state[0] = inv;
            },
        }

        // Apply the matrix mix
        matrix_mix::<T> (state);
    }

    // Second full rounds
    // Process only to R_F -1 iterations. The last iteration does not contain a matrix mix
    for _i in 0..(T::R_F-1) {

        // Add the round constants
        for d in state.iter_mut() {
            //let rc = MNT4753Fr::from_str(ROUND_CST[round_cst_idx]).map_err(|_| ()).unwrap();
            let rc = T::ROUND_CST[round_cst_idx];
            *d += &rc;
            round_cst_idx += 1;
        }

        // Apply the S-BOX to each element of the state vector
        for d in state.iter_mut() {
            let elem_state = (*d).inverse();
            match elem_state {
                None => println!("Field inversion error"),
                Some(inv) => {
                    *d = inv;
                },
            }
        }

        // Apply matrix mix
        matrix_mix::<T> (state);
    }

    // Last full round does not perform the matrix_mix
    // Add the round constants
    for d in state.iter_mut() {
        let rc = T::ROUND_CST[round_cst_idx];
        *d += &rc;
        round_cst_idx += 1;
    }

    // Apply the S-BOX to each element of the state vector
    for d in state.iter_mut() {
        let elem_state = (*d).inverse();
        match elem_state {
            None => println!("Field inversion error"),
            Some(inv) => {
                *d = inv;
            },
        }
    }

}

pub fn poseidon_engine<T: PoseidonParameters>(input: &mut Vec<T::Fr>) -> T::Fr {

    // state is a vector of 3 elements. They are initialized to zero elements
    let mut state = vec![T::ZERO, T::ZERO, T::ZERO];

    // calculate the number of cycles to process the input dividing in portions of rate elements
    let num_cycles = input.len() / T::R;
    // check if the input is a multiple of the rate by calculating the remainder of the division
    let rem = input.len() % T::R;

    // apply permutation to all zeros state vector
    poseidon_perm::<T>(&mut state);

    // index to process the input
    let mut input_idx = 0;
    // iterate of the portions of rate elements
    for _i in 0..num_cycles {
        // add the elements to the state vector. Add rate elements
        for j in 0..T::R {
            state[j] += &input[input_idx];
            input_idx += 1;
        }
        // for application to a 2-1 Merkle tree, add the constant 3 to the third state vector
        state[T::R] += &T::C2;
        // apply permutation after adding the input vector
        poseidon_perm::<T>(&mut state);

    }

    // in case the input is not a multiple of the rate process the remainder part padding a zero
    if rem != 0 {
        state[0] += &input[input_idx];
        state[T::R] += &T::C2;
        // apply permutation after adding the input vector
        poseidon_perm::<T>(&mut state);
    }

    // return the first element of the state vector as the hash digest
    state[0]
}