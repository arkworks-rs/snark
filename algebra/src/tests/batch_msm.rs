#![cfg(feature = "bn_382")]
use crate::bn_382::{Fp, G1Affine, G1Projective};
use algebra_core::{
    msm::{FixedBaseMSM, VariableBaseMSM}, AffineCurve, PrimeField, ProjectiveCurve, UniformRand, Zero,
};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use std::time::Instant;
/*
use super::*;
use crate::fields::bn_382::Fp;
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;
use crate::UniformRand;
use crate::FixedBaseMSM;
*/
use rand::Rng;

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
fn test_with_bn_382() {
    const SAMPLES: usize = 1 << 10;

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES)
        .map(|_| Fp::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = (0..SAMPLES)
        .map(|_| G1Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
    let fast = VariableBaseMSM::multi_scalar_mul(g.as_slice(), v.as_slice());
    let affine = VariableBaseMSM::multi_scalar_mul_affine(g.as_slice(), v.as_slice());

    assert_eq!(naive, fast);
    assert_eq!(naive, affine)
}

#[test]
fn test_with_bn_382_unequal_numbers() {
    const SAMPLES: usize = 1 << 10;

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let v = (0..SAMPLES-1)
        .map(|_| Fp::rand(&mut rng).into_repr())
        .collect::<Vec<_>>();
    let g = (0..SAMPLES)
        .map(|_| G1Projective::rand(&mut rng).into_affine())
        .collect::<Vec<_>>();

    let naive = naive_var_base_msm(g.as_slice(), v.as_slice());
    let fast = VariableBaseMSM::multi_scalar_mul_affine(g.as_slice(), v.as_slice());
    let affine = VariableBaseMSM::multi_scalar_mul_affine(g.as_slice(), v.as_slice());

    assert_eq!(naive, fast);
    assert_eq!(naive, affine)
}

#[test]
fn batch_addition()
{
    let mut length = 1000000;
    let rng = &mut XorShiftRng::seed_from_u64(234872845u64);

    let size_in_bits = <Fp as PrimeField>::size_in_bits();
    let window_size = FixedBaseMSM::get_mul_window_size(length);
    let mut v = FixedBaseMSM::multi_scalar_mul::<G1Projective>
    (
        size_in_bits,
        window_size,
        &FixedBaseMSM::get_window_table
        (
            size_in_bits,
            window_size,
            G1Projective::prime_subgroup_generator()
        ),
        &(0..length).map(|_| Fp::rand(rng)).collect::<Vec<Fp>>(),
    );
    ProjectiveCurve::batch_normalization(&mut v);

    let vector = v.iter().map(|e| e.into_affine()).collect::<Vec<_>>();

    println!();
    println!("{}", "*****BENCHMARKING OPTIMISED AFFINE VS SERIAL JACOBIAN MIXED ADDITION*****");
    loop
    {
        let mut vectors = (0..1000000/length).map(|_| vector[rng.gen_range(0, length/2)..rng.gen_range(length/2, length)].to_vec()).collect::<Vec<_>>();

        //println!("{}{:?}", "Lengths: ".magenta(), vectors.iter().map(|v| v.len()).collect::<Vec<_>>());
        println!();
        println!("{}{:?}", "Length: ", length);
        println!("{}{:?}", "Buckets: ", 1000000/length);

        let start = Instant::now();

        let sum = vectors.iter().map
        (
            |v|
            {
                let mut sum = G1Projective::zero();
                for point in v.iter()
                {
                    sum.add_assign_mixed(point);
                }
                sum.into_affine()
            }
        ).collect::<Vec<_>>();

        let serial = start.elapsed();
        println!("     {}{:?}", "serial: ", serial);
        let start = Instant::now();

        AffineCurve::add_points(&mut vectors);

        let batch = start.elapsed();
        println!("     {}{:?}", "batch: ", batch);
        println!("     {}{:?}", "batch/serial: ", batch.as_secs_f32()/serial.as_secs_f32());
        
        assert_eq!(vectors.iter().map(|v| if v.len() == 0 {G1Affine::zero()} else {v[0]}).collect::<Vec<_>>().iter().eq(sum.iter()), true);

        length = length/10;
        if length == 1 {break}
    }
}

#[test]
fn multiexp()
{
    let mut length = 1000000;
    let rng = &mut XorShiftRng::seed_from_u64(234872845u64);

    let size_in_bits = <Fp as PrimeField>::size_in_bits();
    let window_size = FixedBaseMSM::get_mul_window_size(length);
    let mut v = FixedBaseMSM::multi_scalar_mul::<G1Projective>
    (
        size_in_bits,
        window_size,
        &FixedBaseMSM::get_window_table
        (
            size_in_bits,
            window_size,
            G1Projective::prime_subgroup_generator()
        ),
        &(0..length).map(|_| Fp::rand(rng)).collect::<Vec<Fp>>(),
    );
    ProjectiveCurve::batch_normalization(&mut v);

    let bases = v.iter().map(|e| e.into_affine()).collect::<Vec<_>>();
    let scalars = (0..length).map(|_| Fp::rand(rng).into_repr()).collect::<Vec<_>>();

    println!();
    println!("{}", "*****BENCHMARKING MULTIEXP WITH OPTIMISED AFFINE VS*****");
    println!("{}", "******MULTIEXP WITH SERIAL JACOBIAN MIXED ADDITION******");
    length = 1000000;
    loop
    {
        let base = bases[0..length].to_vec();
        let scalar = scalars[0..length].to_vec();

        println!("{}{:?}", "Length: ", length);

        let start = Instant::now();

        let s1 = VariableBaseMSM::multi_scalar_mul_affine(&base, &scalar);

        let batch = start.elapsed();
        println!("     {}{:?}", "batch: ", batch);
        let start = Instant::now();

        let s2 = VariableBaseMSM::multi_scalar_mul(&base, &scalar);

        let serial = start.elapsed();
        println!("     {}{:?}", "serial: ", serial);
        println!("     {}{:?}", "batch/serial: ", batch.as_secs_f32()/serial.as_secs_f32());

        assert_eq!(s1, s2);
        if length == 1 {break}
        length = length/10;
    }
}

#[test]
fn multiexp_with_tweedledee()
{
    use crate::tweedle::{Fp, dee::Projective as G1Projective};

    let mut length = 1000000;
    let rng = &mut XorShiftRng::seed_from_u64(234872845u64);

    let size_in_bits = <Fp as PrimeField>::size_in_bits();
    let window_size = FixedBaseMSM::get_mul_window_size(length);
    let mut v = FixedBaseMSM::multi_scalar_mul::<G1Projective>
    (
        size_in_bits,
        window_size,
        &FixedBaseMSM::get_window_table
        (
            size_in_bits,
            window_size,
            G1Projective::prime_subgroup_generator()
        ),
        &(0..length).map(|_| Fp::rand(rng)).collect::<Vec<Fp>>(),
    );
    ProjectiveCurve::batch_normalization(&mut v);

    let bases = v.iter().map(|e| e.into_affine()).collect::<Vec<_>>();
    let scalars = (0..length).map(|_| Fp::rand(rng).into_repr()).collect::<Vec<_>>();

    println!();
    println!("{}", "*****BENCHMARKING MULTIEXP WITH OPTIMISED AFFINE VS*****");
    println!("{}", "******MULTIEXP WITH SERIAL JACOBIAN MIXED ADDITION******");
    length = 1000000;
    loop
    {
        let base = bases[0..length].to_vec();
        let scalar = scalars[0..length].to_vec();

        println!("{}{:?}", "Length: ", length);

        let start = Instant::now();

        let s1 = VariableBaseMSM::multi_scalar_mul_affine(&base, &scalar);

        let batch = start.elapsed();
        println!("     {}{:?}", "batch: ", batch);
        let start = Instant::now();

        let s2 = VariableBaseMSM::multi_scalar_mul(&base, &scalar);

        let serial = start.elapsed();
        println!("     {}{:?}", "serial: ", serial);
        println!("     {}{:?}", "batch/serial: ", batch.as_secs_f32()/serial.as_secs_f32());

        assert_eq!(s1, s2);
        if length == 1 {break}
        length = length/10;
    }
}
