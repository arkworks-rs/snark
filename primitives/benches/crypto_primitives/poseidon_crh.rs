#[macro_use]
extern crate criterion;

use criterion::Criterion;
use primitives::{MNT4PoseidonHash, MNT6PoseidonHash, FieldBasedHash, PoseidonBatchHash, BatchFieldBasedHash};
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::UniformRand;
use rand_xorshift::XorShiftRng;
use rand::SeedableRng;
use primitives::crh::poseidon::parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters};

fn poseidon_crh_eval_mnt4(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let input = vec![MNT4753Fr::rand(&mut rng); 100];
    c.bench_function("Poseidon CRH Eval for MNT4", move |b| {
        b.iter(|| {
            MNT4PoseidonHash::evaluate(input.as_slice()).unwrap();
        })
    });
}

fn poseidon_crh_eval_mnt6(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let input = vec![MNT6753Fr::rand(&mut rng); 100];

    c.bench_function("Poseidon CRH Eval for MNT6", move |b| {
        b.iter(|| {
            MNT6PoseidonHash::evaluate(input.as_slice()).unwrap();
        })
    });
}

fn batch_poseidon_crh_eval_mnt4(c: &mut Criterion) {

    type Mnt4BatchPoseidonHash = PoseidonBatchHash<MNT4753Fr, MNT4753PoseidonParameters>;
    //  the number of rounds to test
    let num_hashes = 1000;

    // the vectors that store random input data
    let mut input = Vec::new();

    // the random number generator to generate random input data
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    // we need the double of number of rounds because we have two inputs
    for _ in 0..num_hashes {
        input.push(MNT4753Fr::rand(&mut rng));
        input.push(MNT4753Fr::rand(&mut rng));
    }

    c.bench_function("Batch Poseidon CRH Eval for MNT4 (1000 hashes)", move |b| {
        b.iter(|| {
            Mnt4BatchPoseidonHash::batch_evaluate(&input).unwrap();
        })
    });
}

fn batch_poseidon_crh_eval_mnt6(c: &mut Criterion) {

    type Mnt6BatchPoseidonHash = PoseidonBatchHash<MNT6753Fr, MNT6753PoseidonParameters>;
    //  the number of rounds to test
    let num_hashes = 1000;

    // the vectors that store random input data
    let mut input = Vec::new();

    // the random number generator to generate random input data
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    // we need the double of number of rounds because we have two inputs
    for _ in 0..num_hashes {
        input.push(MNT6753Fr::rand(&mut rng));
        input.push(MNT6753Fr::rand(&mut rng));
    }

    c.bench_function("Batch Poseidon CRH Eval for MNT6 (1000 hashes)", move |b| {
        b.iter(|| {
            Mnt6BatchPoseidonHash::batch_evaluate(&input).unwrap();
        })
    });
}

criterion_group! {
    name = crh_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = poseidon_crh_eval_mnt4, poseidon_crh_eval_mnt6,
}

criterion_group! {
    name = batch_crh_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = batch_poseidon_crh_eval_mnt4, batch_poseidon_crh_eval_mnt6,
}

criterion_main! (
    crh_poseidon_eval, batch_crh_poseidon_eval
);

