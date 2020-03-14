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

fn poseidon_perm_dual<T:PoseidonParameters> (state1: &mut Vec<T::Fr>, state2: &mut Vec<T::Fr>) {

    // index that goes over the round constants
    let mut round_cst_idx = 0;

    // First full rounds
    for _i in 0..T::R_F {

        // Add the round constants to the state vector
        for j in 0..state1.len() {
            //for d in state1.iter_mut() {
            let rc = T::ROUND_CST[round_cst_idx];
            //*d += &rc;
            state1[j] += &rc;
            state2[j] += &rc;
            round_cst_idx += 1;
        }

//        // Apply the S-BOX to each of the elements of the state vector
//        for d in state.iter_mut() {
//
//            // The S-BOX is an inversion function
//            let elem_state = (*d).inverse();
//            match elem_state {
//                None => println!("Field inversion error"),
//                Some(inv) => {
//                    *d = inv;
//                },
//            }
//        }

        // Apply the S-BOX to each of the elements of the state vector
        // Optimization for the inversion S-Box
        // Assuming state vector of 3 elements

//        let abc = state[0]*&state[1]*&state[2];
//        let elem_inv = abc.inverse();
//        let mut abc_inv = T::ZERO;
//        match elem_inv {
//            None => println!("Field inversion error"),
//            Some(inv) => {
//                abc_inv = inv;
//            }
//        }
//        let a_inv = abc_inv*&state[1]*&state[2];
//        let b_inv = abc_inv*&state[0]*&state[2];
//        let c_inv = abc_inv*&state[0]*&state[1];
//
//        state[0] = a_inv;
//        state[1] = b_inv;
//        state[2] = c_inv;
//
//        // Perform the matrix mix
//        matrix_mix::<T> (state);

        let abcdef = state1[0]*&state1[1]*&state1[2]*&state2[0]*&state2[1]*&state2[2];
        let elem_inv = abcdef.inverse();
        let mut abcdef_inv = T::ZERO;
        match elem_inv {
            None => println!("Field inversion error"),
            Some(inv) => {
                abcdef_inv = inv;
            }
        }
        let abc = state1[0]*&state1[1]*&state1[2];
        let cde = state2[0]*&state2[1]*&state2[2];
        let a_inv = abcdef_inv*&state1[1]*&state1[2]*&cde;
        let b_inv = abcdef_inv*&state1[0]*&state1[2]*&cde;
        let c_inv = abcdef_inv*&state1[0]*&state1[1]*&cde;
        let d_inv = abcdef_inv*&abc*&state2[1]*&state2[2];
        let e_inv = abcdef_inv*&abc*&state2[0]*&state2[2];
        let f_inv = abcdef_inv*&abc*&state2[0]*&state2[1];

        state1[0] = a_inv;
        state1[1] = b_inv;
        state1[2] = c_inv;

        state2[0] = d_inv;
        state2[1] = e_inv;
        state2[2] = f_inv;

        // Perform the matrix mix
        matrix_mix::<T> (state1);
        matrix_mix::<T> (state2);
    }

    // Partial rounds
    for _i in 0..T::R_P {

        // Add the round constants to the state vector
        //for d in state.iter_mut() {
        for j in 0..state1.len() {
            let rc = T::ROUND_CST[round_cst_idx];
            state1[j] += &rc;
            state2[j] += &rc;
            //*d += &rc;
            round_cst_idx += 1;
        }

        // Apply S-BOX only to the first element of the state vector
        let ab = state1[0]*&state2[0];
        let elem_inv = ab.inverse();
        let mut ab_inv = T::ZERO;
        match elem_inv {
            None => println!("Field inversion error"),
            Some(inv) => {
                ab_inv = inv;
            }
        }
        let a_inv = ab_inv*&state2[0];
        let b_inv = ab_inv*&state1[0];

        state1[0] = a_inv;
        state2[0] = b_inv;

//        let t1 = (state[0]).inverse();
//        match t1 {
//            None => println!("Field inversion error"),
//            Some(inv) => {
//                state[0] = inv;
//            },
//        }

        // Apply the matrix mix
        matrix_mix::<T> (state1);
        matrix_mix::<T> (state2);
    }

    // Second full rounds
    // Process only to R_F -1 iterations. The last iteration does not contain a matrix mix
    for _i in 0..(T::R_F-1) {

        // Add the round constants
        for j in 0..state1.len() {
            //for d in state1.iter_mut() {
            let rc = T::ROUND_CST[round_cst_idx];
            //*d += &rc;
            state1[j] += &rc;
            state2[j] += &rc;
            round_cst_idx += 1;
        }

//        // Apply the S-BOX to each element of the state vector
//        for d in state.iter_mut() {
//            let elem_state = (*d).inverse();
//            match elem_state {
//                None => println!("Field inversion error"),
//                Some(inv) => {
//                    *d = inv;
//                },
//            }
//        }

        // Optimization for the inversion S-Box
        // Assuming state vector of 3 elements

        let abcdef = state1[0]*&state1[1]*&state1[2]*&state2[0]*&state2[1]*&state2[2];
        let elem_inv = abcdef.inverse();
        let mut abcdef_inv = T::ZERO;
        match elem_inv {
            None => println!("Field inversion error"),
            Some(inv) => {
                abcdef_inv = inv;
            }
        }
        let abc = state1[0]*&state1[1]*&state1[2];
        let cde = state2[0]*&state2[1]*&state2[2];
        let a_inv = abcdef_inv*&state1[1]*&state1[2]*&cde;
        let b_inv = abcdef_inv*&state1[0]*&state1[2]*&cde;
        let c_inv = abcdef_inv*&state1[0]*&state1[1]*&cde;
        let d_inv = abcdef_inv*&abc*&state2[1]*&state2[2];
        let e_inv = abcdef_inv*&abc*&state2[0]*&state2[2];
        let f_inv = abcdef_inv*&abc*&state2[0]*&state2[1];

        state1[0] = a_inv;
        state1[1] = b_inv;
        state1[2] = c_inv;

        state2[0] = d_inv;
        state2[1] = e_inv;
        state2[2] = f_inv;

        // Apply matrix mix
        matrix_mix::<T> (state1);
        matrix_mix::<T> (state2);
    }

    // Last full round does not perform the matrix_mix
    // Add the round constants
    // Add the round constants
    for j in 0..state1.len() {
        //for d in state1.iter_mut() {
        let rc = T::ROUND_CST[round_cst_idx];
        //*d += &rc;
        state1[j] += &rc;
        state2[j] += &rc;
        round_cst_idx += 1;
    }

//    // Apply the S-BOX to each element of the state vector
//    for d in state.iter_mut() {
//        let elem_state = (*d).inverse();
//        match elem_state {
//            None => println!("Field inversion error"),
//            Some(inv) => {
//                *d = inv;
//            },
//        }
//    }

    // Apply the S-BOX to each element of the state vector
    // Optimization for the inversion S-Box
    // Assuming state vector of 3 elements
    let abcdef = state1[0]*&state1[1]*&state1[2]*&state2[0]*&state2[1]*&state2[2];
    let elem_inv = abcdef.inverse();
    let mut abcdef_inv = T::ZERO;
    match elem_inv {
        None => println!("Field inversion error"),
        Some(inv) => {
            abcdef_inv = inv;
        }
    }
    let abc = state1[0]*&state1[1]*&state1[2];
    let cde = state2[0]*&state2[1]*&state2[2];
    let a_inv = abcdef_inv*&state1[1]*&state1[2]*&cde;
    let b_inv = abcdef_inv*&state1[0]*&state1[2]*&cde;
    let c_inv = abcdef_inv*&state1[0]*&state1[1]*&cde;
    let d_inv = abcdef_inv*&abc*&state2[1]*&state2[2];
    let e_inv = abcdef_inv*&abc*&state2[0]*&state2[2];
    let f_inv = abcdef_inv*&abc*&state2[0]*&state2[1];

    state1[0] = a_inv;
    state1[1] = b_inv;
    state1[2] = c_inv;

    state2[0] = d_inv;
    state2[1] = e_inv;
    state2[2] = f_inv;

}

pub fn poseidon_engine_dual<T: PoseidonParameters>(input1: &mut Vec<T::Fr>, input2: &mut Vec<T::Fr>) -> Vec<T::Fr> {

    // state is a vector of 3 elements. They are initialized to zero elements
    let mut state1 = vec![T::ZERO, T::ZERO, T::ZERO];
    let mut state2 = vec![T::ZERO, T::ZERO, T::ZERO];

    // calculate the number of cycles to process the input dividing in portions of rate elements
    // we assume the inputs lengths are the same for input1 and input2
    let num_cycles = input1.len() / T::R;
    // check if the input is a multiple of the rate by calculating the remainder of the division
    let rem = input1.len() % T::R;

    // apply permutation to all zeros state vector
    poseidon_perm_dual::<T>(&mut state1, &mut state2);
    state2 = state1.clone();

    // index to process the input
    let mut input_idx = 0;
    // iterate of the portions of rate elements
    for _i in 0..num_cycles {
        // add the elements to the state vector. Add rate elements
        for j in 0..T::R {
            state1[j] += &input1[input_idx];
            state2[j] += &input2[input_idx];
            input_idx += 1;
        }
        // for application to a 2-1 Merkle tree, add the constant 3 to the third state vector
        state1[T::R] += &T::C2;
        state2[T::R] += &T::C2;
        // apply permutation after adding the input vector
        poseidon_perm_dual::<T>(&mut state1, &mut state2);
    }

    // in case the input is not a multiple of the rate process the remainder part padding a zero
    if rem != 0 {
        state1[0] += &input1[input_idx];
        state1[T::R] += &T::C2;
        state2[0] += &input2[input_idx];
        state2[T::R] += &T::C2;
        // apply permutation after adding the input vector
        poseidon_perm_dual::<T>(&mut state1, &mut state2);
    }

    // return the first element of the state vector as the hash digest
    // state

    let output:Vec<T::Fr> = vec![state1[0],state2[0]];

    output
}
