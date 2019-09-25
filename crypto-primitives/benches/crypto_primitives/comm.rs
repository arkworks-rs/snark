use rand;

#[macro_use]
extern crate criterion;

use algebra::{UniformRand, curves::edwards_bls12::EdwardsProjective as Edwards};
use criterion::Criterion;
use crypto_primitives::commitment::{pedersen::*, CommitmentScheme};

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
            PedersenCommitment::<Edwards, CommWindow>::setup(&mut rng).unwrap()
        })
    });
}

fn pedersen_comm_eval(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let parameters = PedersenCommitment::<Edwards, CommWindow>::setup(&mut rng).unwrap();
    let input = vec![5u8; 128];
    c.bench_function("Pedersen Commitment Eval", move |b| {
        b.iter(|| {
            let rng = &mut rand::thread_rng();
            let commitment_randomness = PedersenRandomness::rand(rng);
            PedersenCommitment::<Edwards, CommWindow>::commit(
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
    config = Criterion::default().sample_size(5);
    targets = pedersen_comm_setup
}

criterion_group! {
    name = comm_eval;
    config = Criterion::default().sample_size(10);
    targets = pedersen_comm_eval
}

criterion_main!(comm_setup, comm_eval);
