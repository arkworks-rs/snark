use accel::*;
mod helpers;
use crate::helpers::create_pseudo_uniform_random_elems;
use algebra::bw6_761::G1Projective;
use algebra_core::{
    curves::{AffineCurve, ProjectiveCurve},
    fields::PrimeField,
    BatchGroupArithmeticSlice, UniformRand,
};
use gpu::gpu_scalar_mul;
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use rayon::prelude::*;

const LOG2_N: usize = 16;
const CHUNK_SIZE: usize = 1024;
const CUDA_GROUP_SIZE: usize = 1 << 5;

pub type G1 = G1Projective;
pub type BigInt = <<G1Projective as ProjectiveCurve>::ScalarField as PrimeField>::BigInt;

fn main() -> error::Result<()> {
    let device = Device::nth(0)?;
    let ctx = device.create_context();

    let _pf = Profiler::start(&ctx);
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    // Allocate memories on GPU
    let n = 1 << LOG2_N;
    let mut exps_h = Vec::with_capacity(n);

    let now = std::time::Instant::now();
    let mut bases_h: Vec<<G1 as ProjectiveCurve>::Affine> =
        create_pseudo_uniform_random_elems(&mut rng, LOG2_N);
    for _ in 0..n {
        exps_h.push(<G1 as ProjectiveCurve>::ScalarField::rand(&mut rng).into_repr());
    }
    println!("Generated random elems: {}us", now.elapsed().as_micros());

    let bases_proj: Vec<_> = bases_h.iter().map(|p| p.into_projective()).collect();

    let now = std::time::Instant::now();
    let mut bases = gpu_scalar_mul(&ctx, &bases_proj[..], &exps_h[..], CUDA_GROUP_SIZE);
    println!("GPU mul: {}us", now.elapsed().as_micros());

    let mut exps_cpu = exps_h.to_vec();
    let now = std::time::Instant::now();
    bases_h
        .par_chunks_mut(CHUNK_SIZE)
        .zip(exps_cpu.par_chunks_mut(CHUNK_SIZE))
        .for_each(|(b, s)| b[..].batch_scalar_mul_in_place(&mut s[..], 4));
    println!("CPU mul: {}us", now.elapsed().as_micros());

    G1::batch_normalization(&mut bases);
    for (b_h, b) in bases_h.into_iter().zip(bases.into_iter()) {
        assert_eq!(b_h, b.into_affine());
    }
    Ok(())
}
