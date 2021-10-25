#[macro_use]
extern crate criterion;

use algebra::{curves::mnt4753::G1Projective as MNT4Projective, UniformRand};
use criterion::Criterion;
use primitives::commitment::{pedersen::*, CommitmentScheme};

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct CommWindow;

impl PedersenWindow for CommWindow {
    const WINDOW_SIZE: usize = 250;
    const NUM_WINDOWS: usize = 8;
}

fn pedersen_comm_setup(c: &mut Criterion) {
    c.bench_function("Pedersen Commitment Setup", move |b| {
        b.iter(|| {
            let mut rng = &mut rand::thread_rng();
            PedersenCommitment::<MNT4Projective, CommWindow>::setup(&mut rng).unwrap()
        })
    });
}

fn pedersen_comm_eval(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let parameters = PedersenCommitment::<MNT4Projective, CommWindow>::setup(&mut rng).unwrap();
    let input = vec![5u8; 128];
    c.bench_function("Pedersen Commitment Eval", move |b| {
        b.iter(|| {
            let rng = &mut rand::thread_rng();
            let commitment_randomness = PedersenRandomness::rand(rng);
            PedersenCommitment::<MNT4Projective, CommWindow>::commit(
                &parameters,
                &input,
                &commitment_randomness,
            )
            .unwrap();
        })
    });
}

criterion_group! {
    name = comm_setup;
    config = Criterion::default().sample_size(20);
    targets = pedersen_comm_setup
}

criterion_group! {
    name = comm_eval;
    config = Criterion::default().sample_size(20);
    targets = pedersen_comm_eval
}

criterion_main!(comm_setup, comm_eval);
