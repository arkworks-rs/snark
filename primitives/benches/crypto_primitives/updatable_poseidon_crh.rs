use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};
use algebra::{
    fields::{
        mnt4753::Fr as MNT4753Fr,
        mnt6753::Fr as MNT6753Fr,
    }
};

use algebra::UniformRand;
use rand_xorshift::XorShiftRng;
use rand::SeedableRng;
use primitives::crh::{
    poseidon::{
        parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters},
        updatable::{UpdatablePoseidonHash, UpdatableBatchPoseidonHash},
    },
    UpdatableFieldBasedHash, UpdatableBatchFieldBasedHash,
};
use primitives::PoseidonParameters;

type UpdatableMNT4PoseidonHash = UpdatablePoseidonHash<MNT4753Fr, MNT4753PoseidonParameters>;
type UpdatableMNT6PoseidonHash = UpdatablePoseidonHash<MNT6753Fr, MNT6753PoseidonParameters>;
type UpdatableMNT4BatchPoseidonHash = UpdatableBatchPoseidonHash<MNT4753Fr, MNT4753PoseidonParameters>;
type UpdatableMNT6BatchPoseidonHash = UpdatableBatchPoseidonHash<MNT6753Fr, MNT6753PoseidonParameters>;

fn poseidon_crh_eval_mnt4(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let input = vec![MNT4753Fr::rand(&mut rng); 100];

    c.bench_function("Poseidon CRH Eval for MNT4", move |b| {
        b.iter(|| {
            poseidon_eval_mnt4(input.as_slice(), UpdatableMNT4PoseidonHash::new(None));
        })
    });
}

fn poseidon_eval_mnt4(inputs: &[MNT4753Fr], mut h: UpdatableMNT4PoseidonHash){
    for i in inputs {
        h.update(i);
    }
    h.finalize();
}

fn poseidon_crh_eval_mnt6(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let input = vec![MNT6753Fr::rand(&mut rng); 100];

    c.bench_function("Poseidon CRH Eval for MNT6", move |b| {
        b.iter(|| {
            poseidon_eval_mnt6(input.as_slice(), UpdatableMNT6PoseidonHash::new(None));
        })
    });
}

fn poseidon_eval_mnt6(inputs: &[MNT6753Fr], mut h: UpdatableMNT6PoseidonHash){
    for i in inputs {
        h.update(i);
    }
    h.finalize();
}

fn batch_poseidon_crh_eval_mnt4(c: &mut Criterion) {

    //  the number of rounds to test
    let num_hashes = 1000;

    // the random number generator to generate random input data
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let input = vec![MNT4753Fr::rand(&mut rng); MNT4753PoseidonParameters::R * num_hashes];

    c.bench_function("Batch Poseidon CRH Eval for MNT4 (1000 hashes)", move |b| {
        b.iter(|| {
            batch_poseidon_eval_mnt4(input.as_slice(), UpdatableMNT4BatchPoseidonHash::new(Some(128)));
        })
    });
}

fn batch_mnt4_hashes_per_core(c: &mut Criterion) {
    const NUM_HASHES: usize = 1000;

    // the random number generator to generate random input data
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let input = [MNT4753Fr::rand(&mut rng); MNT4753PoseidonParameters::R * NUM_HASHES];
    let max = (NUM_HASHES/rayon::current_num_threads()).next_power_of_two();
    let mut cpu_load = 1;
    while cpu_load <= max {
        c.bench_with_input(
            BenchmarkId::new("Batch Poseidon for MNT4 (1000 hashes) with cpu load", cpu_load),
           &cpu_load,
            move |b, cpu_load| {
                b.iter(|| {
                    batch_poseidon_eval_mnt4(
                        &input,
                        UpdatableMNT4BatchPoseidonHash::new(Some(*cpu_load))
                    );
                })
            });
        cpu_load *= 2;
    }
}

fn batch_poseidon_eval_mnt4(inputs: &[MNT4753Fr], mut h: UpdatableMNT4BatchPoseidonHash){
    for i in inputs.chunks(MNT4753PoseidonParameters::R){
        h.update(i);
    }
    h.finalize();
}

fn batch_poseidon_crh_eval_mnt6(c: &mut Criterion) {

    //  the number of rounds to test
    let num_hashes = 1000;

    // the random number generator to generate random input data
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let input = vec![MNT6753Fr::rand(&mut rng); MNT6753PoseidonParameters::R * num_hashes];

    c.bench_function("Batch Poseidon CRH Eval for MNT6 (1000 hashes)", move |b| {
        b.iter(|| {
            batch_poseidon_eval_mnt6(input.as_slice(), UpdatableMNT6BatchPoseidonHash::new(Some(128)));
        })
    });
}

fn batch_mnt6_hashes_per_core(c: &mut Criterion) {
    const NUM_HASHES: usize = 1000;

    // the random number generator to generate random input data
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let input = [MNT6753Fr::rand(&mut rng); MNT6753PoseidonParameters::R * NUM_HASHES];
    let max = (NUM_HASHES/rayon::current_num_threads()).next_power_of_two();
    let mut cpu_load = 1;
    while cpu_load <= max {
        c.bench_with_input(
            BenchmarkId::new("Batch Poseidon for MNT6 (1000 hashes) with cpu_load", cpu_load),
            &cpu_load,
            move |b, cpu_load| {
                b.iter(|| {
                    batch_poseidon_eval_mnt6(
                        &input,
                        UpdatableMNT6BatchPoseidonHash::new(Some(*cpu_load))
                    );
                })
            });
        cpu_load *= 2;
    }
}

fn batch_poseidon_eval_mnt6(inputs: &[MNT6753Fr], mut h: UpdatableMNT6BatchPoseidonHash){
    for i in inputs.chunks(MNT6753PoseidonParameters::R){
        h.update(i);
    }
    h.finalize();
}


criterion_group! {
    name = updatable_crh_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = poseidon_crh_eval_mnt4, poseidon_crh_eval_mnt6,
}

criterion_group! {
    name = updatable_batch_crh_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = batch_poseidon_crh_eval_mnt4, batch_poseidon_crh_eval_mnt6,
}

criterion_group! {
    name = updatable_batch_poseidon_hashes_per_core;
    config = Criterion::default().sample_size(10);
    targets = batch_mnt4_hashes_per_core, batch_mnt6_hashes_per_core
}

criterion_main! (
    updatable_batch_poseidon_hashes_per_core
);