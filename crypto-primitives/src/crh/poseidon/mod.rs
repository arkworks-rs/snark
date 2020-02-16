extern crate hex;
extern crate rand;

use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::{Field, PrimeField, SquareRootField, UniformRand};
use std::ops::Mul;

use std::time::Instant;

use algebra::biginteger::BigInteger768;
use algebra::{to_bytes, ToBytes};
use algebra::field_new;

use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use crate::crh::poseidon::parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters};
use crate::crh::poseidon::poseidon_dual::poseidon_engine_dual;
use crate::crh::poseidon::poseidon_original::poseidon_engine;
use crate::crh::poseidon::mul_inv::bench_mul_inv;
use std::marker::PhantomData;
use crate::crh::Batched2to1CRH;
use std::error::Error;

pub mod poseidon_original;
pub mod poseidon_dual;
pub mod parameters;
pub mod mul_inv;

pub trait PoseidonParameters: 'static{

    type Fr: PrimeField + SquareRootField + Into<<Self::Fr as PrimeField>::BigInt>;

    const T: usize;  // Number of S-Boxesb
    const R_F:i32;   // Number of full rounds
    const R_P:i32;   // Number of partial rounds
    const R:usize;   // The rate of the hash function
    const ZERO:Self::Fr;   // The zero element in the field
    const C2:Self::Fr;     // The constant 3 to add in the position corresponding to the capacity
    const ROUND_CST: &'static[Self::Fr];  // Array of round constants
    const MDS_CST: &'static[Self::Fr];  // The MDS matrix

}

pub struct PoseidonCRH<P: PoseidonParameters> {
    parameters: PhantomData<P>

    //fn evaluate(parameters: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error> {}

}

impl<P: PoseidonParameters> Batched2to1CRH for PoseidonCRH<P> {
    const INPUT_NUM_PAIRS: usize = unimplemented!();
    type Output = P::Fr;
    type Parameters = P;

    fn evaluate(input: &[u8]) -> Result<Self::Output, Error> {}

}

//fn print_cst () {
//    let cst = Fr::from_str("3").map_err(|_| ()).unwrap();
//    println!("{:?}", cst);
//
//    let mut d_out = to_bytes!(cst).unwrap();
//    d_out.reverse();
//
//    println!("constant = {:?}", hex::encode(d_out));
//
//}

//fn print_cst () {
//    for i in 0..195 {
//        let cst = Fr::from_str(ROUND_CST[i]).map_err(|_| ()).unwrap();
//        println!("{:?}", cst);
//    }
//}
//
//fn print_mds () {
//    for i in 0..9 {
//        let cst = Fr::from_str(MDS[i]).map_err(|_| ()).unwrap();
//        println!("{:?}", cst);
//    }
//}

#[test]
fn test_cst() {

    //  the number of rounds to test
    let num_rounds = 1000;

    // the vectors that store random input data
    let mut vec_elem_4753:Vec<MNT4753Fr> = Vec::new();
    let mut vec_elem_6753:Vec<MNT6753Fr> = Vec::new();

    // the random number generator to generate random input data
    // let mut rng = &mut thread_rng();
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    // we need the double of number of rounds because we have two inputs
    for _ in 0..(2*num_rounds) {
        vec_elem_4753.push(MNT4753Fr::rand(&mut rng));
        vec_elem_6753.push(MNT6753Fr::rand(&mut rng));
    }

    // =============================================================================
    // Calculate Poseidon Hash for mnt4753 combining 2 hashes
    let now_4753_dual = Instant::now();

    let mut output_4753_dual = Vec::new();

    for i in 0..num_rounds/2 {
        let mut input1 = vec![vec_elem_4753[4*i], vec_elem_4753[4*i+1]];
        let mut input2 = vec![vec_elem_4753[4*i+2], vec_elem_4753[4*i+3]];

        // Call the poseidon hash
        let output = poseidon_engine_dual::<MNT4753PoseidonParameters>(&mut input1, &mut input2);
        output_4753_dual.push(output[0]);
        output_4753_dual.push(output[1]);
    }

    let new_now_4753_dual  = Instant::now();
    // =============================================================================
    // =============================================================================
    // Calculate Poseidon Hash for mnt4753 original
    let now_4753_single = Instant::now();

    let mut output_4753_single = Vec::new();

    for i in 0..num_rounds {
        let mut input1 = vec![vec_elem_4753[2*i], vec_elem_4753[2*i+1]];

        // Call the poseidon hash
        let output = poseidon_engine::<MNT4753PoseidonParameters>(&mut input1);
        output_4753_single.push(output);
    }

    let new_now_4753_single  = Instant::now();
    // =============================================================================


    // =============================================================================
    // Calculate Poseidon Hash for mnt6753 combining 2 hashes
    let now_6753_dual = Instant::now();

    let mut output_6753_dual = Vec::new();

    for i in 0..num_rounds/2 {
        //let mut input = Vec::new();
        let mut input1 = vec![vec_elem_6753[4*i],vec_elem_6753[4*i+1]];
        let mut input2 = vec![vec_elem_6753[4*i+2],vec_elem_6753[4*i+3]];

        // Call the poseidon hash
        let output = poseidon_engine_dual::<MNT6753PoseidonParameters>(&mut input1, &mut input2);
        output_6753_dual.push(output[0]);
        output_6753_dual.push(output[1]);
    }
    let new_now_6753_dual  = Instant::now();
    // =============================================================================
    // =============================================================================
    // Calculate Poseidon Hash for mnt6753 original
    let now_6753_single = Instant::now();

    let mut output_6753_single = Vec::new();

    for i in 0..num_rounds {
        //let mut input = Vec::new();
        let mut input1 = vec![vec_elem_6753[2*i],vec_elem_6753[2*i+1]];

        // Call the poseidon hash
        let output = poseidon_engine::<MNT6753PoseidonParameters>(&mut input1);
        output_6753_single.push(output);
    }
    let new_now_6753_single  = Instant::now();
    // =============================================================================

    // =============================================================================
    // Compare results
    for i in 0..num_rounds {
        if output_4753_dual[i] != output_4753_single[i] {
            println!("Hash outputs, position {}, for MNT4 are not equal.",i);
        }
    }
    println!("End comparison for MNT4.");
    for i in 0..num_rounds {
        if output_6753_dual[i] != output_6753_single[i] {
            println!("Hash outputs, position {}, for MNT6 are not equal.",i);
        }
    }
    println!("End comparison for MNT6.");


//    // =============================================================================
//    // Print the result of for mnt4753
//    for i in 0..num_rounds {
//
//        // Reverse order to output the data
//        let mut d_in_0 = to_bytes!(vec_elem_4753[2*i]).unwrap();
//        d_in_0.reverse();
//        let mut d_in_1 = to_bytes!(vec_elem_4753[2*i + 1]).unwrap();
//        d_in_1.reverse();
//
//        let mut d_out_dual = to_bytes!(output_4753_dual[i]).unwrap();
//        d_out_dual.reverse();
//        let mut d_out_single = to_bytes!(output_4753_single[i]).unwrap();
//        d_out_single.reverse();
//
//        println!("input[0] = {:?}", hex::encode(d_in_0));
//        println!("input[1] = {:?}", hex::encode(d_in_1));
//        println!("hash MNT4753 single = {:?}", hex::encode(d_out_single));
//        println!("hash MNT4753 dual   = {:?}", hex::encode(d_out_dual));
//
//    }
//    // =============================================================================
//
//    // =============================================================================
//    // Print the result for mnt6753
//    for i in 0..num_rounds {
//        // Reverse order to output the data
//        let mut d_in_0 = to_bytes!(vec_elem_6753[2*i]).unwrap();
//        d_in_0.reverse();
//        let mut d_in_1 = to_bytes!(vec_elem_6753[2*i + 1]).unwrap();
//        d_in_1.reverse();
//
//        let mut d_out_dual = to_bytes!(output_6753_dual[i]).unwrap();
//        d_out_dual.reverse();
//        let mut d_out_single = to_bytes!(output_6753_single[i]).unwrap();
//        d_out_single.reverse();
//
//        println!("input[0] = {:?}", hex::encode(d_in_0));
//        println!("input[1] = {:?}", hex::encode(d_in_1));
//        println!("hash MNT6753 single = {:?}", hex::encode(d_out_single));
//        println!("hash MNT6753 dual   = {:?}", hex::encode(d_out_dual));
//    }
//    // =============================================================================

    // =============================================================================
    // Report the timing results

    let duration_4753_single =  new_now_4753_single.duration_since(now_4753_single);
    println!("Time for {} rounds MNT4753 single = {:?}", num_rounds, duration_4753_single.as_millis());

    let duration_4753_dual =  new_now_4753_dual.duration_since(now_4753_dual);
    println!("Time for {} rounds MNT4753 dual   = {:?}", num_rounds, duration_4753_dual.as_millis());

    let duration_6753_single =  new_now_6753_single.duration_since(now_6753_single);
    println!("Time for {} rounds MNT6753 single = {:?}", num_rounds, duration_6753_single.as_millis());

    let duration_6753_dual =  new_now_6753_dual.duration_since(now_6753_dual);
    println!("Time for {} rounds MNT6753 dual   = {:?}", num_rounds, duration_6753_dual.as_millis());

    // =============================================================================

    //bench_mul_inv();

}
