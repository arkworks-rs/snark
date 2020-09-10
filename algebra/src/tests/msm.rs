use algebra_core::{
    msm::VariableBaseMSM, AffineCurve, PrimeField, ProjectiveCurve, UniformRand, Zero,
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use crate::tests::helpers::create_pseudo_uniform_random_elems;

fn _naive_var_base_msm<G: AffineCurve>(
    bases: &[G],
    scalars: &[<G::ScalarField as PrimeField>::BigInt],
) -> G::Projective {
    let mut acc = G::Projective::zero();

    for (base, scalar) in bases.iter().zip(scalars.iter()) {
        acc += &base.mul(*scalar);
    }
    acc
}

#[allow(unused)]
pub fn test_msm<G: AffineCurve>() {
    const MAX_LOGN: usize = 14;
    const SAMPLES: usize = 1 << MAX_LOGN;

    let _lol = G::Projective::zero();
    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES)
        .map(|_| G::ScalarField::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = create_pseudo_uniform_random_elems::<G, XorShiftRng>(&mut rng, MAX_LOGN);

    // let naive = naive_var_base_msm(g.as_slice(), v.as_slice());

    let now = std::time::Instant::now();
    let even_faster = VariableBaseMSM::multi_scalar_mul_batched(
        g.as_slice(),
        v.as_slice(),
        <G::ScalarField as PrimeField>::size_in_bits(),
    );
    println!(
        "new MSM for {} elems: {:?}",
        SAMPLES,
        now.elapsed().as_micros()
    );

    let now = std::time::Instant::now();
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());
    println!(
        "old MSM for {} elems: {:?}",
        SAMPLES,
        now.elapsed().as_micros()
    );

    assert_eq!(even_faster.into_affine(), fast.into_affine());
}
