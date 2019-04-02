use rand;

#[macro_use]
extern crate criterion;

use algebra::curves::edwards_bls12::EdwardsAffine as Edwards;
use criterion::Criterion;
use dpc::crypto_primitives::crh::{pedersen::*, FixedLengthCRH};

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct HashWindow;

impl PedersenWindow for HashWindow {
    const WINDOW_SIZE: usize = 250;
    const NUM_WINDOWS: usize = 8;
}

fn pedersen_crh_setup(c: &mut Criterion) {
    c.bench_function("Pedersen CRH Setup", move |b| {
        b.iter(|| {
            let mut rng = rand::OsRng::new().unwrap();
            PedersenCRH::<Edwards, HashWindow>::setup(&mut rng).unwrap()
        })
    });
}

fn pedersen_crh_eval(c: &mut Criterion) {
    let mut rng = rand::OsRng::new().unwrap();
    let parameters = PedersenCRH::<Edwards, HashWindow>::setup(&mut rng).unwrap();
    let input = vec![5u8; 128];
    c.bench_function("Pedersen CRH Eval", move |b| {
        b.iter(|| {
            PedersenCRH::<Edwards, HashWindow>::evaluate(&parameters, &input).unwrap();
        })
    });
}

criterion_group! {
    name = crh_setup;
    config = Criterion::default().sample_size(5);
    targets = pedersen_crh_setup
}

criterion_group! {
    name = crh_eval;
    config = Criterion::default().sample_size(10);
    targets = pedersen_crh_eval
}

criterion_main!(crh_setup, crh_eval);
