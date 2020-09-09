#[macro_use]
extern crate criterion;

use criterion::Criterion;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::Field;
use algebra::UniformRand;
use rand_xorshift::XorShiftRng;
use rand::SeedableRng;
use primitives::{
    crh::{
        MNT4PoseidonHash, MNT6PoseidonHash,
        batched_crh::{
            MNT4BatchPoseidonHash, MNT6BatchPoseidonHash
        },
    },
    merkle_tree::{
        field_based_mht::{
            poseidon::{
                MNT4753_MHT_POSEIDON_PARAMETERS, MNT6753_MHT_POSEIDON_PARAMETERS,
            },
            FieldBasedMerkleTree,
            FieldBasedOptimizedMHT, FieldBasedMerkleTreePrecomputedEmptyConstants,
            FieldBasedMerkleTreeParameters, BatchFieldBasedMerkleTreeParameters,
        },
    }
};

#[derive(Clone, Debug)]
struct MNT4753FieldBasedMerkleTreeParams;
impl FieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
    type Data = MNT4753Fr;
    type H = MNT4PoseidonHash;
    const MERKLE_ARITY: usize = 2;
    const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> = Some(MNT4753_MHT_POSEIDON_PARAMETERS);
}

impl BatchFieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
    type BH = MNT4BatchPoseidonHash;
}

type MNT4PoseidonMHT = FieldBasedOptimizedMHT<MNT4753FieldBasedMerkleTreeParams>;

#[derive(Clone, Debug)]
struct MNT6753FieldBasedMerkleTreeParams;
impl FieldBasedMerkleTreeParameters for MNT6753FieldBasedMerkleTreeParams {
    type Data = MNT6753Fr;
    type H = MNT6PoseidonHash;
    const MERKLE_ARITY: usize = 2;
    const EMPTY_HASH_CST: Option<FieldBasedMerkleTreePrecomputedEmptyConstants<'static, Self::H>> = Some(MNT6753_MHT_POSEIDON_PARAMETERS);
}

impl BatchFieldBasedMerkleTreeParameters for MNT6753FieldBasedMerkleTreeParams {
    type BH = MNT6BatchPoseidonHash;
}

type MNT6PoseidonMHT = FieldBasedOptimizedMHT<MNT6753FieldBasedMerkleTreeParams>;

const BENCH_HEIGHT: usize = 11;

fn batch_poseidon_mht_eval_mnt4_full(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch Full Poseidon MHT Eval for MNT4 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..num_leaves {
                tree.append(MNT4753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_full(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch Full Poseidon MHT Eval for MNT6 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..num_leaves {
                tree.append(MNT6753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt4_3_4(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch 3/4 Poseidon MHT Eval for MNT4 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..(num_leaves * 3)/4 {
                tree.append(MNT4753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_3_4(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch 3/4 Poseidon MHT Eval for MNT6 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..(num_leaves * 3)/4 {
                tree.append(MNT6753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt4_half(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch half Poseidon MHT Eval for MNT4 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/2 {
                tree.append(MNT4753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_half(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch half Poseidon MHT Eval for MNT6 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/2 {
                tree.append(MNT6753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt4_1_4(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch 1/4 Poseidon MHT Eval for MNT4 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/4 {
                tree.append(MNT4753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_1_4(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch 1/4 Poseidon MHT Eval for MNT6 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/4 {
                tree.append(MNT6753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt4_interleaved(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch interleaved Poseidon MHT Eval for MNT4 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/3 {
                tree.append(MNT4753Fr::zero());
            }
            for _ in 0..(num_leaves * 2)/3 {
                tree.append(MNT4753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_interleaved(c: &mut Criterion) {

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT);

    c.bench_function("Batch interleaved Poseidon MHT Eval for MNT6 (2^10 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/3 {
                tree.append(MNT6753Fr::zero());
            }
            for _ in 0..(num_leaves * 2)/3 {
                tree.append(MNT6753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}


criterion_group! {
    name = mht_poseidon_eval;
    config = Criterion::default().sample_size(100);
    targets = batch_poseidon_mht_eval_mnt4_full, batch_poseidon_mht_eval_mnt6_full, batch_poseidon_mht_eval_mnt4_3_4, batch_poseidon_mht_eval_mnt6_3_4, batch_poseidon_mht_eval_mnt4_half, batch_poseidon_mht_eval_mnt6_half, batch_poseidon_mht_eval_mnt4_1_4, batch_poseidon_mht_eval_mnt6_1_4, batch_poseidon_mht_eval_mnt4_interleaved, batch_poseidon_mht_eval_mnt6_interleaved
}

criterion_main! (
    mht_poseidon_eval
);
