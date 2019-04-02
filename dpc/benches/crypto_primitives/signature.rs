#[macro_use]
extern crate criterion;

mod affine {
    use algebra::curves::edwards_bls12::EdwardsAffine as Edwards;
    use blake2::Blake2s;
    use criterion::Criterion;
    use dpc::crypto_primitives::signature::{schnorr::*, SignatureScheme};
    use rand::{self, Rng};

    type SchnorrEdwards = SchnorrSignature<Edwards, Blake2s>;
    fn schnorr_signature_setup(c: &mut Criterion) {
        c.bench_function("SchnorrEdwardsAffine: Setup", move |b| {
            b.iter(|| {
                let mut rng = rand::OsRng::new().unwrap();
                SchnorrEdwards::setup(&mut rng).unwrap()
            })
        });
    }

    fn schnorr_signature_keygen(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();

        c.bench_function("SchnorrEdwardsAffine: KeyGen", move |b| {
            b.iter(|| {
                let mut rng = rand::OsRng::new().unwrap();
                SchnorrEdwards::keygen(&parameters, &mut rng).unwrap()
            })
        });
    }

    fn schnorr_signature_sign(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
        let (_, sk) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
        let message = [100u8; 128];

        c.bench_function("SchnorrEdwardsAffine: Sign", move |b| {
            b.iter(|| {
                let mut rng = rand::OsRng::new().unwrap();
                SchnorrEdwards::sign(&parameters, &sk, &message, &mut rng).unwrap()
            })
        });
    }

    fn schnorr_signature_verify(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
        let (pk, sk) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
        let message = [100u8; 128];
        let signature = SchnorrEdwards::sign(&parameters, &sk, &message, &mut rng).unwrap();

        c.bench_function("SchnorrEdwardsAffine: Verify", move |b| {
            b.iter(|| SchnorrEdwards::verify(&parameters, &pk, &message, &signature).unwrap())
        });
    }

    fn schnorr_signature_randomize_pk(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
        let (pk, _) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
        let randomness: [u8; 32] = rng.gen();

        c.bench_function("SchnorrEdwardsAffine: Randomize PubKey", move |b| {
            b.iter(|| SchnorrEdwards::randomize_public_key(&parameters, &pk, &randomness).unwrap())
        });
    }

    fn schnorr_signature_randomize_signature(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
        let (_, sk) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
        let randomness: [u8; 32] = rng.gen();
        let message = [100u8; 128];
        let signature = SchnorrEdwards::sign(&parameters, &sk, &message, &mut rng).unwrap();

        c.bench_function("SchnorrEdwardsAffine: Randomize Signature", move |b| {
            b.iter(|| {
                SchnorrEdwards::randomize_signature(&parameters, &signature, &randomness).unwrap()
            })
        });
    }
    criterion_group! {
        name = schnorr_sig_affine;
        config = Criterion::default().sample_size(20);
        targets = schnorr_signature_setup, schnorr_signature_keygen, schnorr_signature_sign,
                  schnorr_signature_verify, schnorr_signature_randomize_pk, schnorr_signature_randomize_signature
    }
}

mod projective {
    use algebra::curves::edwards_bls12::EdwardsProjective as Edwards;
    use blake2::Blake2s;
    use criterion::Criterion;
    use dpc::crypto_primitives::signature::{schnorr::*, SignatureScheme};
    use rand::{self, Rng};

    type SchnorrEdwards = SchnorrSignature<Edwards, Blake2s>;
    fn schnorr_signature_setup(c: &mut Criterion) {
        c.bench_function("SchnorrEdwardsProjective: Setup", move |b| {
            b.iter(|| {
                let mut rng = rand::OsRng::new().unwrap();
                SchnorrEdwards::setup(&mut rng).unwrap()
            })
        });
    }

    fn schnorr_signature_keygen(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();

        c.bench_function("SchnorrEdwardsProjective: KeyGen", move |b| {
            b.iter(|| {
                let mut rng = rand::OsRng::new().unwrap();
                SchnorrEdwards::keygen(&parameters, &mut rng).unwrap()
            })
        });
    }

    fn schnorr_signature_sign(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
        let (_, sk) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
        let message = [100u8; 128];

        c.bench_function("SchnorrEdwardsProjective: Sign", move |b| {
            b.iter(|| {
                let mut rng = rand::OsRng::new().unwrap();
                SchnorrEdwards::sign(&parameters, &sk, &message, &mut rng).unwrap()
            })
        });
    }

    fn schnorr_signature_verify(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
        let (pk, sk) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
        let message = [100u8; 128];
        let signature = SchnorrEdwards::sign(&parameters, &sk, &message, &mut rng).unwrap();

        c.bench_function("SchnorrEdwardsProjective: Verify", move |b| {
            b.iter(|| SchnorrEdwards::verify(&parameters, &pk, &message, &signature).unwrap())
        });
    }

    fn schnorr_signature_randomize_pk(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
        let (pk, _) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
        let randomness: [u8; 32] = rng.gen();

        c.bench_function("SchnorrEdwardsProjective: Randomize PubKey", move |b| {
            b.iter(|| SchnorrEdwards::randomize_public_key(&parameters, &pk, &randomness).unwrap())
        });
    }

    fn schnorr_signature_randomize_signature(c: &mut Criterion) {
        let mut rng = rand::OsRng::new().unwrap();
        let parameters = SchnorrEdwards::setup(&mut rng).unwrap();
        let (_, sk) = SchnorrEdwards::keygen(&parameters, &mut rng).unwrap();
        let randomness: [u8; 32] = rng.gen();
        let message = [100u8; 128];
        let signature = SchnorrEdwards::sign(&parameters, &sk, &message, &mut rng).unwrap();

        c.bench_function("SchnorrEdwardsProjective: Randomize Signature", move |b| {
            b.iter(|| {
                SchnorrEdwards::randomize_signature(&parameters, &signature, &randomness).unwrap()
            })
        });
    }
    criterion_group! {
        name = schnorr_sig_projective;
        config = Criterion::default().sample_size(20);
        targets = schnorr_signature_setup, schnorr_signature_keygen, schnorr_signature_sign,
                  schnorr_signature_verify, schnorr_signature_randomize_pk, schnorr_signature_randomize_signature
    }
}
use crate::{affine::schnorr_sig_affine, projective::schnorr_sig_projective};
criterion_main!(schnorr_sig_affine, schnorr_sig_projective);
