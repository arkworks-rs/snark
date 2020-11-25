use criterion::{criterion_group, criterion_main, Criterion};
use algebra::{
    fields::{
        mnt4753::Fr as MNT4753Fr,
        mnt6753::Fr as MNT6753Fr,
        bn_382::Fr as BN382Fr,
        bn_382::Fq as BN382Fq,
    }
};

use algebra::UniformRand;
use rand_xorshift::XorShiftRng;
use rand::SeedableRng;
use primitives::crh::{
    poseidon::parameters::*,
    FieldBasedHash,
};

fn poseidon_crh_eval_mnt4(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let mut h = MNT4PoseidonHash::init(None);

    c.bench_function("Poseidon CRH Eval for MNT4", move |b| {
        b.iter(|| {
            for _ in 0..2000 {
                h.update(MNT4753Fr::rand(&mut rng));
            }
            h.finalize();
        })
    });
}

fn poseidon_crh_eval_mnt6(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let mut h = MNT6PoseidonHash::init(None);

    c.bench_function("Poseidon CRH Eval for MNT6", move |b| {
        b.iter(|| {
            for _ in 0..2000 {
                h.update(MNT6753Fr::rand(&mut rng));
            }
            h.finalize();
        })
    });
}

fn poseidon_crh_eval_bn382fr(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let mut h = BN382FrPoseidonHash::init(None);

    c.bench_function("Poseidon CRH Eval for BN382Fr", move |b| {
        b.iter(|| {
            for _ in 0..2000 {
                h.update(BN382Fr::rand(&mut rng));
            }
            h.finalize();
        })
    });
}

fn poseidon_crh_eval_bn382fq(c: &mut Criterion) {

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let mut h = BN382FqPoseidonHash::init(None);

    c.bench_function("Poseidon CRH Eval for BN382Fq", move |b| {
        b.iter(|| {
            for _ in 0..2000 {
                h.update(BN382Fq::rand(&mut rng));
            }
            h.finalize();
        })
    });
}

criterion_group! {
    name = crh_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = poseidon_crh_eval_mnt4, poseidon_crh_eval_mnt6,
              poseidon_crh_eval_bn382fq, poseidon_crh_eval_bn382fr,
}

criterion_main! (
    crh_poseidon_eval
);