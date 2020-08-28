use rand;

#[macro_use]
extern crate criterion;

use algebra::ed_on_bls12_377::EdwardsProjective as Edwards;
use criterion::Criterion;
use crypto_primitives::crh::{pedersen::*, FixedLengthCRH};

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct HashWindow;

impl Window for HashWindow {
    const WINDOW_SIZE: usize = 250;
    const NUM_WINDOWS: usize = 8;
}

fn pedersen_crh_setup(c: &mut Criterion) {
    c.bench_function("Pedersen CRH Setup", move |b| {
        b.iter(|| {
            let mut rng = &mut rand::thread_rng();
            CRH::<Edwards, HashWindow>::setup(&mut rng).unwrap()
        })
    });
}

fn pedersen_crh_eval(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let parameters = CRH::<Edwards, HashWindow>::setup(&mut rng).unwrap();
    let input = vec![5u8; 128];
    c.bench_function("Pedersen CRH Eval", move |b| {
        b.iter(|| CRH::<Edwards, HashWindow>::evaluate(&parameters, &input).unwrap())
    });
}

criterion_group! {
    name = crh_setup;
    config = Criterion::default().sample_size(10);
    targets = pedersen_crh_setup
}

criterion_group! {
    name = crh_eval;
    config = Criterion::default().sample_size(10);
    targets = pedersen_crh_eval
}

criterion_main!(crh_setup, crh_eval);
