use algebra::curves::mnt6753::G1Projective as MNT6G1Projective;
use algebra::fields::mnt4753::Fr as MNT4Fr;
use algebra::UniformRand;
use criterion::Criterion;
use primitives::{
    crh::{bowe_hopwood::BoweHopwoodPedersenCRH, pedersen::PedersenWindow, MNT4PoseidonHash},
    vrf::{ecvrf::FieldBasedEcVrf, FieldBasedVrf},
    FixedLengthCRH,
};

#[macro_use]
extern crate criterion;

#[derive(Clone)]
struct TestWindow {}
impl PedersenWindow for TestWindow {
    const WINDOW_SIZE: usize = 128;
    const NUM_WINDOWS: usize = 2;
}

type BHMNT6 = BoweHopwoodPedersenCRH<MNT6G1Projective, TestWindow>;
type EcVrfMNT4 = FieldBasedEcVrf<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash, BHMNT6>;

fn ecvrf_keygen(c: &mut Criterion) {
    c.bench_function("FieldSchnorrMNT4: KeyGen", move |b| {
        b.iter(|| {
            let mut rng = &mut rand::thread_rng();
            EcVrfMNT4::keygen(&mut rng)
        })
    });
}

fn ecvrf_prove(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let pp = <BHMNT6 as FixedLengthCRH>::setup(rng).unwrap();
    let (pk, sk) = EcVrfMNT4::keygen(&mut rng);
    let message = MNT4Fr::rand(rng);

    c.bench_function("FieldSchnorrMNT4: Sign", move |b| {
        b.iter(|| {
            let mut rng = &mut rand::thread_rng();
            EcVrfMNT4::prove(&mut rng, &pp, &pk, &sk, message).unwrap()
        })
    });
}

fn ecvrf_verify(c: &mut Criterion) {
    let mut rng = &mut rand::thread_rng();
    let pp = <BHMNT6 as FixedLengthCRH>::setup(rng).unwrap();
    let (pk, sk) = EcVrfMNT4::keygen(&mut rng);
    let message = MNT4Fr::rand(rng);
    let proof = EcVrfMNT4::prove(&mut rng, &pp, &pk, &sk, message).unwrap();

    c.bench_function("FieldSchnorrMNT4: Proof To Hash", move |b| {
        b.iter(|| EcVrfMNT4::proof_to_hash(&pp, &pk, message, &proof).unwrap())
    });
}

criterion_group! {
    name = ecvrf;
    config = Criterion::default().sample_size(20);
    targets = ecvrf_keygen, ecvrf_prove, ecvrf_verify
}

criterion_main!(ecvrf);
