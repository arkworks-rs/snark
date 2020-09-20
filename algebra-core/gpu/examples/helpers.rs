use algebra_core::{
    cfg_chunks_mut, AffineCurve, BatchGroupArithmeticSlice, BigInteger64, ProjectiveCurve,
    UniformRand,
};
use rand::{distributions::Uniform, prelude::Distribution, Rng};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub fn create_pseudo_uniform_random_elems<C: AffineCurve, R: Rng>(
    rng: &mut R,
    max_logn: usize,
) -> Vec<C> {
    const AFFINE_BATCH_SIZE: usize = 4096;
    println!("Starting");
    let now = std::time::Instant::now();
    // Generate pseudorandom group elements
    let step = Uniform::new(0, 1 << (max_logn + 5));
    let elem = C::Projective::rand(rng).into_affine();
    let mut random_elems = vec![elem; 1 << max_logn];
    let mut scalars: Vec<BigInteger64> = (0..1 << max_logn)
        .map(|_| BigInteger64::from(step.sample(rng)))
        .collect();
    cfg_chunks_mut!(random_elems, AFFINE_BATCH_SIZE)
        .zip(cfg_chunks_mut!(scalars, AFFINE_BATCH_SIZE))
        .for_each(|(e, s)| {
            e[..].batch_scalar_mul_in_place::<BigInteger64>(&mut s[..], 1);
        });

    println!("Initial generation: {:?}", now.elapsed().as_micros());
    random_elems
}
