use crate::crh::poseidon::PoseidonParameters;

use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::{Field, PrimeField, SquareRootField, UniformRand};
use std::ops::Mul;

use std::time::Instant;

use algebra::biginteger::BigInteger768;
use algebra::{to_bytes, ToBytes};
use algebra::field_new;
use rand_xorshift::XorShiftRng;
use super::rand::SeedableRng;

// Function that benchmarks computation time of multiplication and division
pub fn bench_mul_inv() {

    //  the number of rounds to test
    let num_rounds = 1000;

    // the vectors that store random input data
    let mut vec_elem_4753: Vec<MNT4753Fr> = Vec::new();
    let mut vec_elem_6753: Vec<MNT6753Fr> = Vec::new();

    // the random number generator to generate random input data
    // let mut rng = &mut thread_rng();
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    // we need the double of number of rounds because we have two inputs
    for _ in 0..(2 * num_rounds) {
        vec_elem_4753.push(MNT4753Fr::rand(&mut rng));
        vec_elem_6753.push(MNT6753Fr::rand(&mut rng));
    }

    // =============================================================================
    // Calculate multiplication time MNT4
    let mut vec_mul_out_4753 = Vec::new();
    let now_mul_4753 = Instant::now();

    for i in 0..num_rounds {
        let o1 = vec_elem_4753[2 * i] * &vec_elem_4753[2 * i + 1];
        vec_mul_out_4753.push(o1);
    }

    let new_now_mul_4753 = Instant::now();
    println!("Last mul MNT4 = {:?}", vec_mul_out_4753[num_rounds - 1]);
    // =============================================================================
    // =============================================================================
    // Calculate multiplication time MNT6
    let mut vec_mul_out_6753 = Vec::new();
    let now_mul_6753 = Instant::now();

    for i in 0..num_rounds {
        let o1 = vec_elem_6753[2 * i] * &vec_elem_6753[2 * i + 1];
        vec_mul_out_6753.push(o1);
    }

    let new_now_mul_6753 = Instant::now();
    println!("Last mul MNT6 = {:?}", vec_mul_out_6753[num_rounds - 1]);
    // =============================================================================


    // =============================================================================
    // Calculate inversion time MNT4
    let mut vec_inv_out_4753 = Vec::new();
    let now_inv_4753 = Instant::now();

    for i in 0..num_rounds {
        let o1 = vec_elem_4753[i].inverse();
        vec_inv_out_4753.push(o1);
    }

    let new_now_inv_4753 = Instant::now();
    println!("Last inv = {:?}", vec_inv_out_4753[num_rounds - 1]);
    // =============================================================================
    // =============================================================================
    // Calculate inversion time MNT6
    let mut vec_inv_out_6753 = Vec::new();
    let now_inv_6753 = Instant::now();

    for i in 0..num_rounds {
        let o1 = vec_elem_6753[i].inverse();
        vec_inv_out_6753.push(o1);
    }

    let new_now_inv_6753 = Instant::now();
    println!("Last inv = {:?}", vec_inv_out_6753[num_rounds - 1]);
    // =============================================================================

    // Report timing results
    let duration_mul_4753 = new_now_mul_4753.duration_since(now_mul_4753);
    println!("Time for {} mul  rounds MNT4753 = {:?}", num_rounds, duration_mul_4753.as_micros());

    let duration_mul_6753 = new_now_mul_6753.duration_since(now_mul_6753);
    println!("Time for {} mul  rounds MNT6753 = {:?}", num_rounds, duration_mul_6753.as_micros());

    let duration_inv_4753 = new_now_inv_4753.duration_since(now_inv_4753);
    println!("Time for {} inv  rounds MNT4753 = {:?}", num_rounds, duration_inv_4753.as_micros());

    let duration_inv_6753 = new_now_inv_6753.duration_since(now_inv_6753);
    println!("Time for {} inv  rounds MNT6753 = {:?}", num_rounds, duration_inv_6753.as_micros());
}