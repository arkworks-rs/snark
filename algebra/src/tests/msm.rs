#![cfg(any(feature = "bls12_381", feature = "bw6_761", feature = "bn254"))]
#[cfg(feature = "bls12_381")]
use crate::bls12_381::{Fr, G1Affine, G1Projective};
#[cfg(all(feature = "bn254", not(feature = "bls12_381")))]
use crate::bn254::{Fr, G1Affine, G1Projective};
#[cfg(all(feature = "bw6_761", not(feature = "bls12_381")))]
use crate::bw6_761::{Fr, G1Affine, G1Projective};

use algebra_core::{
    msm::VariableBaseMSM, AffineCurve, PrimeField, ProjectiveCurve, UniformRand, Zero,
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use crate::tests::helpers::create_pseudo_uniform_random_elems;

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
fn test() {
    test_msm::<G1Affine>();
}

fn test_msm<G: AffineCurve>() {
    const MAX_LOGN: usize = 15;
    const SAMPLES: usize = 1 << MAX_LOGN;

    let _lol = G1Projective::zero();
    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES)
        .map(|_| Fr::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = create_pseudo_uniform_random_elems::<G1Affine>(&mut rng, MAX_LOGN);

    // let naive = naive_var_base_msm(g.as_slice(), v.as_slice());

    let now = std::time::Instant::now();
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());
    println!(
        "old MSM for {} elems: {:?}",
        SAMPLES,
        now.elapsed().as_micros()
    );
    let now = std::time::Instant::now();
    let even_faster = VariableBaseMSM::multi_scalar_mul_batched(
        g.as_slice(),
        v.as_slice(),
        <<G1Affine as AffineCurve>::ScalarField as PrimeField>::size_in_bits(),
    );
    println!(
        "new MSM for {} elems: {:?}",
        SAMPLES,
        now.elapsed().as_micros()
    );

    assert_eq!(even_faster.into_affine(), fast.into_affine());
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
