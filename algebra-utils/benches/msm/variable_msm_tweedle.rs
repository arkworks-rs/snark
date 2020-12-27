#[macro_use]
extern crate criterion;

#[macro_use]
extern crate bench_utils;

use criterion::{BenchmarkId, Criterion, BatchSize};

use algebra::tweedle::{
    Fp, dee::Projective as G1Projective,
    dee::Affine as G1Affine
};

use algebra::{UniformRand, PrimeField, ProjectiveCurve};
use rand::rngs::OsRng;
use std::time::{SystemTime, UNIX_EPOCH};


use algebra_core::msm::VariableBaseMSM;
use rand_xorshift::XorShiftRng;
use std::fs::File;
use rand::SeedableRng;
use algebra_core::{ToBytes, BigInteger256, FromBytes};

const PARAM_C: usize = 16;

fn variable_msm_affine_sd(c: &mut Criterion) {

    let mut group = c.benchmark_group("variable_base_msm_affine_sd-tweedle_dee-variable number of bases = number of scalars");
    let samples = (14..=23).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    for &samples in samples.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(samples), &samples, |b, _samples| {
            b.iter_batched(|| {
                let (v, g) = load_data(samples);
                (v, g)
            },
                           |(v, g)| {
                               add_to_trace!(
             || format!("****************{}*******************", samples),
             || format!("--->START TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                               VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(), PARAM_C);
                               add_to_trace!(
             || format!("****************{}*******************", samples),
             || format!("--->END TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                           },
                           BatchSize::PerIteration);
        });
    }
}

fn variable_msm_affine(c: &mut Criterion) {
    let mut group = c.benchmark_group("variable_base_msm_affine-tweedle_dee-variable number of bases = number of scalars");
    let samples = (14..=23).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    for &samples in samples.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(samples), &samples, |b, _samples| {
            b.iter_batched(|| {
                let (v, g) = load_data(samples);
                (v, g)
            },
                           |(v, g)| {
                               add_to_trace!(
             || format!("****************{}*******************", samples),
             || format!("--->START TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                               VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(), PARAM_C);
                               add_to_trace!(
             || format!("****************{}*******************", samples),
             || format!("--->END TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
                     );
                           },
                           BatchSize::PerIteration);
        });
    }
}

fn variable_msm_fast(c: &mut Criterion) {

    let mut group = c.benchmark_group("variable_base_msm_fast-tweedle_dee-variable number of bases = number of scalars");
    let samples = (14..=23).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    for &samples in samples.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(samples), &samples, |b, _samples| {
            b.iter_batched(|| {
                let (v, g) = load_data(samples);
                (v, g)
            },
                           |(v, g)| {
                               add_to_trace!(
             || format!("****************{}*******************", samples),
             || format!("--->START TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                               VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(), PARAM_C);
                               add_to_trace!(
             || format!("****************{}*******************", samples),
             || format!("--->END TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                           },
                           BatchSize::PerIteration);
        });
    }
}

fn variable_msm_affine_sd_4(c: &mut Criterion) {

    let mut group = c.benchmark_group("variable_base_msm_affine_sd-tweedle_dee-variable number of bases = number of scalars");
    let samples = (17..=23).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    let c = 4;

    for &samples in samples.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(samples.to_string()+"_"+ &*c.to_string()), &samples, |b, _samples| {
            b.iter_batched(|| {
                let (v, g) = load_data(samples);
                (v, g)
            },
                           |(v, g)| {
                               add_to_trace!(
             || format!("****************{}-{}*******************", samples, c),
             || format!("--->START TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                               VariableBaseMSM::multi_scalar_mul_affine_sd_c(g.as_slice(), v.as_slice(), c);
                               add_to_trace!(
             || format!("****************{}-{}*******************", samples, c),
             || format!("--->END TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                           },
                           BatchSize::PerIteration);
        });
    }
}

fn variable_msm_affine_4(c: &mut Criterion) {

    let mut group = c.benchmark_group("variable_base_msm_affine-tweedle_dee-variable number of bases = number of scalars");
    let samples = (17..=23).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    let c = 4;

    for &samples in samples.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(samples.to_string()+"_"+ &*c.to_string()), &samples, |b, _samples| {
            b.iter_batched(|| {
                let (v, g) = load_data(samples);
                (v, g)
            },
                           |(v, g)| {
                               add_to_trace!(
             || format!("****************{}-{}*******************", samples, c),
             || format!("--->START TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                               VariableBaseMSM::multi_scalar_mul_affine_c(g.as_slice(), v.as_slice(), c);
                               add_to_trace!(
             || format!("****************{}-{}*******************", samples, c),
             || format!("--->END TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                           },
                           BatchSize::PerIteration);
        });
    }
}

fn variable_msm_fast_4(c: &mut Criterion) {

    let mut group = c.benchmark_group("variable_base_msm_fast-tweedle_dee-variable number of bases = number of scalars");
    let samples = (17..=23).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    let c = 4;

    for &samples in samples.iter() {
        group.bench_with_input(BenchmarkId::from_parameter(samples.to_string()+"_"+ &*c.to_string()), &samples, |b, _samples| {
        //group.bench_with_input(BenchmarkId::from_parameter(samples), &samples, |b, _samples| {
            b.iter_batched(|| {
                let (v, g) = load_data(samples);
                (v, g)
            },
                           |(v, g)| {
                               add_to_trace!(
             || format!("****************{}-{}*******************", samples, c),
             || format!("--->START TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                               VariableBaseMSM::multi_scalar_mul_c(g.as_slice(), v.as_slice(), c);
                               add_to_trace!(
             || format!("****************{}-{}*******************", samples, c),
             || format!("--->END TIMESTAMP: {:?}", SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_secs())
        );
                           },
                           BatchSize::PerIteration);
        });
    }
}

fn load_data(samples: usize) -> (Vec<BigInteger256>,Vec<G1Affine>) {

    let mut fs = File::open("./scalars_bases_tweedle").unwrap();
    let mut v = Vec::with_capacity(samples);
    let mut g = Vec::with_capacity(samples);

    for _i in 0..samples {
        let elem1 = BigInteger256::read(&mut fs).unwrap();
        let elem2 = G1Affine::read(&mut fs).unwrap();
        v.push(elem1);
        g.push(elem2);
    }
    (v, g)
}


#[allow(dead_code)]
fn generate_data() {

    const SAMPLES: usize = 1<<23;

    let mut rng = XorShiftRng::seed_from_u64(234872845u64);

    let mut fs = File::create("./scalars_bases_tweedle").unwrap();

    for _i in 0..SAMPLES {
        let elem = Fp::rand(&mut rng).into_repr();
        elem.write(&mut fs).unwrap();
        let elem = G1Projective::rand(&mut rng).into_affine();
        elem.write(&mut fs).unwrap();
    }
}


criterion_group! {
    name = variable_msm_eval_tweedle;
    config = Criterion::default().sample_size(10);
    targets = variable_msm_fast, variable_msm_affine, variable_msm_affine_sd
    // targets = variable_msm_fast_4, variable_msm_affine_4, variable_msm_affine_sd_4
}

criterion_main! (
    variable_msm_eval_tweedle
);

