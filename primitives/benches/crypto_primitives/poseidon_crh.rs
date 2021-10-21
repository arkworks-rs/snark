use algebra::fields::{
    bn_382::Fq as BN382Fq, bn_382::Fr as BN382Fr, mnt4753::Fr as MNT4753Fr,
    mnt6753::Fr as MNT6753Fr, tweedle::Fq as tweedleFq, tweedle::Fr as tweedleFr,
};
use criterion::{criterion_group, criterion_main, Criterion};

use algebra::UniformRand;
use primitives::crh::{poseidon::parameters::*, FieldBasedHash};
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

fn poseidon_crh_eval_mnt4(c: &mut Criterion) {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let samples = 2000;
    let mut h = MNT4PoseidonHash::init_constant_length(samples, None);

    c.bench_function("Poseidon CRH Eval for MNT4", move |b| {
        b.iter(|| {
            for _ in 0..samples {
                h.update(MNT4753Fr::rand(&mut rng));
            }
            h.finalize().unwrap();
            h.reset(None);
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
            h.finalize().unwrap();
            h.reset(None);
        })
    });
}

fn poseidon_crh_eval_bn382fr(c: &mut Criterion) {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let samples = 2000;
    let mut h = BN382FrPoseidonHash::init_constant_length(samples, None);

    c.bench_function("Poseidon CRH Eval for BN382Fr", move |b| {
        b.iter(|| {
            for _ in 0..samples {
                h.update(BN382Fr::rand(&mut rng));
            }
            h.finalize().unwrap();
            h.reset(None);
        })
    });
}

fn poseidon_crh_eval_bn382fq(c: &mut Criterion) {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let samples = 2000;
    let mut h = BN382FqPoseidonHash::init_constant_length(samples, None);

    c.bench_function("Poseidon CRH Eval for BN382Fq", move |b| {
        b.iter(|| {
            for _ in 0..samples {
                h.update(BN382Fq::rand(&mut rng));
            }
            h.finalize().unwrap();
            h.reset(None);
        })
    });
}

fn poseidon_crh_eval_tweedlefr(c: &mut Criterion) {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let samples = 2000;
    let mut h = TweedleFrPoseidonHash::init_constant_length(samples, None);

    c.bench_function("Poseidon CRH Eval for tweedleFr", move |b| {
        b.iter(|| {
            for _ in 0..samples {
                h.update(tweedleFr::rand(&mut rng));
            }
            h.finalize().unwrap();
            h.reset(None);
        })
    });
}

fn poseidon_crh_eval_tweedlefq(c: &mut Criterion) {
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    let samples = 2000;
    let mut h = TweedleFqPoseidonHash::init_constant_length(samples, None);

    c.bench_function("Poseidon CRH Eval for tweedleFq", move |b| {
        b.iter(|| {
            for _ in 0..samples {
                h.update(tweedleFq::rand(&mut rng));
            }
            h.finalize().unwrap();
            h.reset(None);
        })
    });
}

criterion_group! {
    name = crh_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = poseidon_crh_eval_mnt4, poseidon_crh_eval_mnt6,
              poseidon_crh_eval_bn382fq, poseidon_crh_eval_bn382fr,
              poseidon_crh_eval_tweedlefq, poseidon_crh_eval_tweedlefr,
}

criterion_main!(crh_poseidon_eval);
