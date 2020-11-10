use algebra_core::{
    cuda::scalar_mul::{GPUScalarMul, GPUScalarMulSlice, MAX_GROUP_ELEM_BYTES},
    AffineCurve, BatchGroupArithmeticSlice, PrimeField, UniformRand, Zero,
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use crate::{cfg_chunks_mut, tests::helpers::create_pseudo_uniform_random_elems};

const CHUNK_SIZE: usize = 1 << 12;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

#[allow(unused)]
pub fn test_cuda_scalar_mul<G: AffineCurve>() {
    #[cfg(not(feature = "big_n"))]
    const MAX_LOGN: usize = 14;
    #[cfg(feature = "big_n")]
    const MAX_LOGN: usize = 20;

    let cuda_group_size = 1 << 5;
    if core::mem::size_of::<G>() >= MAX_GROUP_ELEM_BYTES {
        println!("Group size too large to run on GPU, defaulting to CPU-only implementation");
    }

    const SAMPLES: usize = 1 << MAX_LOGN;

    let _lol = G::Projective::zero();
    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let exps_h = (0..SAMPLES)
        .map(|_| G::ScalarField::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let mut bases_h = create_pseudo_uniform_random_elems::<G, XorShiftRng>(&mut rng, MAX_LOGN);

    let mut bases_d = bases_h.to_vec();
    let mut exps_cpu = exps_h.to_vec();

    let now = std::time::Instant::now();
    cfg_chunks_mut!(bases_h, CHUNK_SIZE)
        .zip(cfg_chunks_mut!(exps_cpu, CHUNK_SIZE))
        .for_each(|(b, s)| b[..].batch_scalar_mul_in_place(&mut s[..], 4));
    println!("CPU mul: {}us", now.elapsed().as_micros());

    <G as AffineCurve>::Projective::clear_gpu_profiling_data();

    let mut junk_data = bases_d.to_vec();
    for _ in 0..3 {
        let now = std::time::Instant::now();
        &mut junk_data[..].cpu_gpu_scalar_mul(&exps_h[..], cuda_group_size, CHUNK_SIZE);
        println!("CPU + GPU mul: {}us", now.elapsed().as_micros());
    }
    let now = std::time::Instant::now();
    &mut bases_d[..].cpu_gpu_scalar_mul(&exps_h[..], cuda_group_size, CHUNK_SIZE);
    println!("CPU + GPU mul: {}us", now.elapsed().as_micros());

    for (b_h, b_d) in bases_h.into_iter().zip(bases_d.into_iter()) {
        assert_eq!(b_h, b_d);
    }
}
