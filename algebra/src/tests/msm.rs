use algebra_core::{
    msm::VariableBaseMSM, AffineCurve, PrimeField, ProjectiveCurve, UniformRand, Zero,
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

fn naive_var_base_msm<G: AffineCurve>(
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
    #[cfg(not(feature = "big_n"))]
    const MAX_LOGN: usize = 14;
    #[cfg(feature = "big_n")]
    const MAX_LOGN: usize = 21;

    const SAMPLES: usize = 1 << MAX_LOGN;

    let _lol = G::Projective::zero();
    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES)
        .map(|_| G::ScalarField::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();

    let g = (0..SAMPLES)
        .map(|_| G::Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let naive = naive_var_base_msm(g.as_slice(), v.as_slice());

    let now = std::time::Instant::now();
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());
    println!(
        "old MSM for {} elems: {:?}",
        SAMPLES,
        now.elapsed().as_micros()
    );

    assert_eq!(naive.into_affine(), fast.into_affine());
}
