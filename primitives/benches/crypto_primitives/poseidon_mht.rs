#[macro_use]
extern crate criterion;

use algebra::{
    fields::{
        mnt4753::Fr as MNT4753Fr, mnt6753::Fr as MNT6753Fr
    },
    Field, UniformRand
};

use primitives::{
    crh::parameters::{
        MNT4PoseidonHash, MNT6PoseidonHash,
        MNT4BatchPoseidonHash, MNT6BatchPoseidonHash
    },
    merkle_tree::{
        field_based_mht::{
            parameters::{
                MNT4753_MHT_POSEIDON_PARAMETERS, MNT6753_MHT_POSEIDON_PARAMETERS,
            },
            FieldBasedMerkleTree,
            FieldBasedOptimizedMHT, FieldBasedMerkleTreePrecomputedZeroConstants,
            FieldBasedMerkleTreeParameters, BatchFieldBasedMerkleTreeParameters,
        },
    }
};

use criterion::{
    Criterion, BenchmarkId
};

use rand_xorshift::XorShiftRng;
use rand::SeedableRng;

#[derive(Clone, Debug)]
struct MNT4753FieldBasedMerkleTreeParams;
impl FieldBasedMerkleTreeParameters for MNT4753FieldBasedMerkleTreeParams {
    type Data = MNT4753Fr;
    type H = MNT4PoseidonHash;
    const MERKLE_ARITY: usize = 2;
    const ZERO_NODE_CST: Option<FieldBasedMerkleTreePrecomputedZeroConstants<'static, Self::H>> = Some(MNT4753_MHT_POSEIDON_PARAMETERS);
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
    const ZERO_NODE_CST: Option<FieldBasedMerkleTreePrecomputedZeroConstants<'static, Self::H>> = Some(MNT6753_MHT_POSEIDON_PARAMETERS);
}

impl BatchFieldBasedMerkleTreeParameters for MNT6753FieldBasedMerkleTreeParams {
    type BH = MNT6BatchPoseidonHash;
}

type MNT6PoseidonMHT = FieldBasedOptimizedMHT<MNT6753FieldBasedMerkleTreeParams>;

const BENCH_HEIGHT: usize = 11;

fn batch_poseidon_mht_eval_mnt4_full(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch Full Poseidon MHT Eval for MNT4 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..num_leaves {
                tree.append(MNT4753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_full(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch Full Poseidon MHT Eval for MNT6 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..num_leaves {
                tree.append(MNT6753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt4_3_4(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch 3/4 Poseidon MHT Eval for MNT4 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..(num_leaves * 3)/4 {
                tree.append(MNT4753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_3_4(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch 3/4 Poseidon MHT Eval for MNT6 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..(num_leaves * 3)/4 {
                tree.append(MNT6753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt4_half(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch half Poseidon MHT Eval for MNT4 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/2 {
                tree.append(MNT4753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_half(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch half Poseidon MHT Eval for MNT6 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/2 {
                tree.append(MNT6753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt4_1_4(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch 1/4 Poseidon MHT Eval for MNT4 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/4 {
                tree.append(MNT4753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_1_4(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch 1/4 Poseidon MHT Eval for MNT6 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/4 {
                tree.append(MNT6753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt4_interleaved(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch interleaved Poseidon MHT Eval for MNT4 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/3 {
                tree.append(MNT4753Fr::zero()).unwrap();
            }
            for _ in 0..(num_leaves * 2)/3 {
                tree.append(MNT4753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_interleaved(c: &mut Criterion) {

    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT, num_leaves).unwrap();

    c.bench_function(format!("Batch interleaved Poseidon MHT Eval for MNT6 ({} leaves)", num_leaves).as_str(), move |b| {
        b.iter(|| {
            for _ in 0..num_leaves/3 {
                tree.append(MNT6753Fr::zero()).unwrap();
            }
            for _ in 0..(num_leaves * 2)/3 {
                tree.append(MNT6753Fr::rand(&mut rng)).unwrap();
            }
            tree.finalize_in_place();
            tree.reset();
        })
    });
}

/// Let's create a full tree with different processing_step sizes and bench the total time
fn batch_poseidon_mht_tune_processing_step_mnt4(c: &mut Criterion) {
    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut processing_steps = Vec::with_capacity(BENCH_HEIGHT - 1);
    for i in 0..BENCH_HEIGHT {
        processing_steps.push(2usize.pow(i as u32));
    }

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut group = c.benchmark_group(
        format!("tune processing_step_size for MNT4 with a tree of height {}", BENCH_HEIGHT)
    );

    for processing_step in processing_steps.iter() {
        let mut tree = MNT4PoseidonMHT::init(BENCH_HEIGHT, *processing_step).unwrap();
        group.bench_with_input(BenchmarkId::from_parameter(processing_step), processing_step, |b, _processing_step| {
            b.iter(|| {
                for _ in 0..num_leaves {
                    tree.append(MNT4753Fr::rand(&mut rng)).unwrap();
                }
                tree.finalize_in_place();
                tree.reset();
            });
        });
    }
}

/// Let's create a full tree with different processing_step sizes and bench the total time
fn batch_poseidon_mht_tune_processing_step_mnt6(c: &mut Criterion) {
    let num_leaves = 2usize.pow(BENCH_HEIGHT as u32 - 1);

    let mut processing_steps = Vec::with_capacity(BENCH_HEIGHT - 1);
    for i in 0..BENCH_HEIGHT {
        processing_steps.push(2usize.pow(i as u32));
    }

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut group = c.benchmark_group(
        format!("tune processing_step_size for MNT6 with a tree of height {}", BENCH_HEIGHT)
    );

    for processing_step in processing_steps.iter() {
        let mut tree = MNT6PoseidonMHT::init(BENCH_HEIGHT, *processing_step).unwrap();

        group.bench_with_input(BenchmarkId::from_parameter(processing_step), processing_step, |b, _processing_step| {
            b.iter(|| {
                for _ in 0..num_leaves {
                    tree.append(MNT6753Fr::rand(&mut rng)).unwrap();
                }
                tree.finalize_in_place();
                tree.reset();
            });
        });
    }
}

criterion_group! {
    name = mht_poseidon_eval;
    config = Criterion::default().sample_size(100);
    targets =
        batch_poseidon_mht_eval_mnt4_full, batch_poseidon_mht_eval_mnt6_full,
        batch_poseidon_mht_eval_mnt4_3_4, batch_poseidon_mht_eval_mnt6_3_4,
        batch_poseidon_mht_eval_mnt4_half, batch_poseidon_mht_eval_mnt6_half,
        batch_poseidon_mht_eval_mnt4_1_4, batch_poseidon_mht_eval_mnt6_1_4,
        batch_poseidon_mht_eval_mnt4_interleaved, batch_poseidon_mht_eval_mnt6_interleaved
}

criterion_group!{
    name = mht_poseidon_tuning;
    config = Criterion::default().sample_size(10);
    targets =
        batch_poseidon_mht_tune_processing_step_mnt4, batch_poseidon_mht_tune_processing_step_mnt6
}

criterion_main! (
    mht_poseidon_tuning
);