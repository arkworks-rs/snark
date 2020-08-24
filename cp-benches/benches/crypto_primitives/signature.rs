#[macro_use]
extern crate criterion;

use algebra::ed_on_bls12_377::EdwardsProjective as Edwards;
use blake2::Blake2s;
use criterion::Criterion;
use crypto_primitives::signature::{schnorr::*, SignatureScheme};
use rand::{self, Rng};

type SchnorrEdwards = Schnorr<Edwards, Blake2s>;
fn schnorr_signature_setup(c: &mut Criterion) {
    c.bench_function("SchnorrEdwards: Setup", move |b| {
        b.iter(|| {
            let mut rng = &mut rand::thread_rng();
            SchnorrEdwards::setup(&mut rng).unwrap()
        })
    });
}

fn schnorr_signature_keygen(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let parameters = SchnorrEdwards::setup(&mut rng).unwrap();

    c.bench_function("SchnorrEdwards: KeyGen", move |b| {
        b.iter(|| {
            let mut rng = &mut rand::thread_rng();
            SchnorrEdwards::keygen(&parameters, &mut rng).unwrap()
        })
    });
}

fn schnorr_signature_sign(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
    let (_, sk) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
    let message = [100u8; 128];

    c.bench_function("SchnorrEdwards: Sign", move |b| {
        b.iter(|| {
            let mut rng = &mut rand::thread_rng();
            SchnorrEdwards::sign(&parameters, &sk, &message, &mut rng).unwrap()
        })
    });
}

fn schnorr_signature_verify(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
    let (pk, sk) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
    let message = [100u8; 128];
    let signature = SchnorrEdwards::sign(&parameters, &sk, &message, &mut rng).unwrap();

    c.bench_function("SchnorrEdwards: Verify", move |b| {
        b.iter(|| SchnorrEdwards::verify(&parameters, &pk, &message, &signature).unwrap())
    });
}

fn schnorr_signature_randomize_pk(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
    let (pk, _) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
    let randomness: [u8; 32] = rng.gen();

    c.bench_function("SchnorrEdwards: Randomize PubKey", move |b| {
        b.iter(|| SchnorrEdwards::randomize_public_key(&parameters, &pk, &randomness).unwrap())
    });
}

fn schnorr_signature_randomize_signature(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
    let (_, sk) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
    let randomness: [u8; 32] = rng.gen();
    let message = [100u8; 128];
    let signature = SchnorrEdwards::sign(&parameters, &sk, &message, &mut rng).unwrap();

    c.bench_function("SchnorrEdwards: Randomize Signature", move |b| {
        b.iter(|| {
            SchnorrEdwards::randomize_signature(&parameters, &signature, &randomness).unwrap()
        })
    });
}
criterion_group! {
    name = schnorr_sig;
    config = Criterion::default().sample_size(20);
    targets = schnorr_signature_setup, schnorr_signature_keygen, schnorr_signature_sign,
                schnorr_signature_verify, schnorr_signature_randomize_pk, schnorr_signature_randomize_signature
}
criterion_main!(schnorr_sig);
