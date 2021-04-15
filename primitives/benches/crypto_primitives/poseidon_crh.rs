use criterion::{criterion_group, criterion_main, Criterion};
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
        PoseidonHash
    },
    FieldBasedHash,
};

type MNT4PoseidonHash = PoseidonHash<MNT4753Fr, MNT4753PoseidonParameters>;
type MNT6PoseidonHash = PoseidonHash<MNT6753Fr, MNT6753PoseidonParameters>;

fn poseidon_crh_eval_mnt4(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let samples = 2000;
    let mut h = MNT4PoseidonHash::init_constant_length(samples, None);

    c.bench_function("Poseidon CRH Eval for MNT4", move |b| {
        b.iter(|| {
            for _ in 0..samples {
                h.update(MNT4753Fr::rand(&mut rng));
            }
            h.finalize();
        })
    });
}

fn poseidon_crh_eval_mnt6(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let samples = 2000;
    let mut h = MNT6PoseidonHash::init_constant_length(samples, None);

    c.bench_function("Poseidon CRH Eval for MNT6", move |b| {
        b.iter(|| {
            for _ in 0..samples {
                h.update(MNT6753Fr::rand(&mut rng));
            }
            h.finalize();
        })
    });
}

criterion_group! {
    name = crh_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = poseidon_crh_eval_mnt4, poseidon_crh_eval_mnt6,
}

criterion_main! (
    crh_poseidon_eval
);