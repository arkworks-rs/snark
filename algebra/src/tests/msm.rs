#![cfg(feature = "bls12_381")]
use algebra_core::{UniformRand, msm::VariableBaseMSM, ProjectiveCurve, AffineCurve, PrimeField, Zero};
use crate::bls12_381::{G1Projective, Fr};
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

#[test]
fn test_with_bls12() {
    const SAMPLES: usize = 1 << 10;

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES)
        .map(|_| Fr::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = (0..SAMPLES)
        .map(|_| G1Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

    assert_eq!(naive.into_affine(), fast.into_affine());
}

#[test]
fn test_with_bls12_unequal_numbers() {
    const SAMPLES: usize = 1 << 10;

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES - 1)
        .map(|_| Fr::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = (0..SAMPLES)
        .map(|_| G1Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());

    assert_eq!(naive.into_affine(), fast.into_affine());
}
