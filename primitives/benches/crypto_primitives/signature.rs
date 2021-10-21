#[macro_use]
extern crate criterion;

mod projective {
    use algebra::curves::mnt4753::G1Projective as MNT4Projective;
    use blake2::Blake2s;
    use criterion::Criterion;
    use primitives::signature::{schnorr::*, SignatureScheme};
    use rand::{self, Rng};

    type SchnorrMNT4Affine = SchnorrSignature<MNT4Projective, Blake2s>;
    fn schnorr_signature_setup(c: &mut Criterion) {
        c.bench_function("SchnorrMNT4Projective: Setup", move |b| {
            b.iter(|| {
                let mut rng = &mut rand::thread_rng();
                SchnorrMNT4Affine::setup(&mut rng).unwrap()
            })
        });
    }

    fn schnorr_signature_keygen(c: &mut Criterion) {
        let mut rng = &mut rand::thread_rng();
        let parameters = SchnorrMNT4Affine::setup(&mut rng).unwrap();

        c.bench_function("SchnorrMNT4Projective: KeyGen", move |b| {
            b.iter(|| {
                let mut rng = &mut rand::thread_rng();
                SchnorrMNT4Affine::keygen(&parameters, &mut rng).unwrap()
            })
        });
    }

    fn schnorr_signature_sign(c: &mut Criterion) {
        let mut rng = &mut rand::thread_rng();
        let parameters = SchnorrMNT4Affine::setup(&mut rng).unwrap();
        let (_, sk) = SchnorrMNT4Affine::keygen(&parameters, &mut rng).unwrap();
        let message = [100u8; 128];

        c.bench_function("SchnorrMNT4Projective: Sign", move |b| {
            b.iter(|| {
                let mut rng = &mut rand::thread_rng();
                SchnorrMNT4Affine::sign(&parameters, &sk, &message, &mut rng).unwrap()
            })
        });
    }

    fn schnorr_signature_verify(c: &mut Criterion) {
        let mut rng = &mut rand::thread_rng();
        let parameters = SchnorrMNT4Affine::setup(&mut rng).unwrap();
        let (pk, sk) = SchnorrMNT4Affine::keygen(&parameters, &mut rng).unwrap();
        let message = [100u8; 128];
        let signature = SchnorrMNT4Affine::sign(&parameters, &sk, &message, &mut rng).unwrap();

        c.bench_function("SchnorrMNT4Projective: Verify", move |b| {
            b.iter(|| SchnorrMNT4Affine::verify(&parameters, &pk, &message, &signature).unwrap())
        });
    }

    fn schnorr_signature_randomize_pk(c: &mut Criterion) {
        let mut rng = &mut rand::thread_rng();
        let parameters = SchnorrMNT4Affine::setup(&mut rng).unwrap();
        let (pk, _) = SchnorrMNT4Affine::keygen(&parameters, &mut rng).unwrap();
        let randomness: [u8; 32] = rng.gen();

        c.bench_function("SchnorrMNT4Projective: Randomize PubKey", move |b| {
            b.iter(|| {
                SchnorrMNT4Affine::randomize_public_key(&parameters, &pk, &randomness).unwrap()
            })
        });
    }

    fn schnorr_signature_randomize_signature(c: &mut Criterion) {
        let mut rng = &mut rand::thread_rng();
        let parameters = SchnorrMNT4Affine::setup(&mut rng).unwrap();
        let (_, sk) = SchnorrMNT4Affine::keygen(&parameters, &mut rng).unwrap();
        let randomness: [u8; 32] = rng.gen();
        let message = [100u8; 128];
        let signature = SchnorrMNT4Affine::sign(&parameters, &sk, &message, &mut rng).unwrap();

        c.bench_function("SchnorrMNT4Projective: Randomize Signature", move |b| {
            b.iter(|| {
                SchnorrMNT4Affine::randomize_signature(&parameters, &signature, &randomness)
                    .unwrap()
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

mod field_impl {
    use algebra::{
        curves::mnt6753::G1Projective as MNT6G1Projective, fields::mnt4753::Fr as MNT4Fr,
        UniformRand,
    };
    use criterion::Criterion;
    use primitives::{
        crh::MNT4PoseidonHash,
        signature::{schnorr::field_based_schnorr::*, FieldBasedSignatureScheme},
    };

    type SchnorrMNT4Fr =
        FieldBasedSchnorrSignatureScheme<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash>;

    fn schnorr_signature_keygen(c: &mut Criterion) {
        c.bench_function("FieldSchnorrMNT4: KeyGen", move |b| {
            b.iter(|| {
                let mut rng = &mut rand::thread_rng();
                SchnorrMNT4Fr::keygen(&mut rng)
            })
        });
    }

    fn schnorr_signature_sign(c: &mut Criterion) {
        let mut rng = &mut rand::thread_rng();
        let (pk, sk) = SchnorrMNT4Fr::keygen(&mut rng);
        let message = MNT4Fr::rand(rng);

        c.bench_function("FieldSchnorrMNT4: Sign", move |b| {
            b.iter(|| {
                let mut rng = &mut rand::thread_rng();
                SchnorrMNT4Fr::sign(&mut rng, &pk, &sk, message).unwrap()
            })
        });
    }

    fn schnorr_signature_verify(c: &mut Criterion) {
        let mut rng = &mut rand::thread_rng();
        let (pk, sk) = SchnorrMNT4Fr::keygen(&mut rng);
        let message = MNT4Fr::rand(rng);
        let signature = SchnorrMNT4Fr::sign(&mut rng, &pk, &sk, message).unwrap();

        c.bench_function("FieldSchnorrMNT4: Verify", move |b| {
            b.iter(|| SchnorrMNT4Fr::verify(&pk, message, &signature).unwrap())
        });
    }

    criterion_group! {
        name = field_based_schnorr_sig;
        config = Criterion::default().sample_size(20);
        targets = schnorr_signature_keygen,
                  schnorr_signature_sign, schnorr_signature_verify,
    }
}

use crate::{field_impl::field_based_schnorr_sig, projective::schnorr_sig_projective};
criterion_main!(schnorr_sig_projective, field_based_schnorr_sig);
