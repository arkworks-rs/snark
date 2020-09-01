#[macro_use]
extern crate criterion;

use criterion::Criterion;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::fields::Field;
use algebra::UniformRand;
use rand_xorshift::XorShiftRng;
use rand::SeedableRng;
use primitives::crh::poseidon::parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters};
use primitives::merkle_tree::{
    field_based_mht::poseidon::{
        PoseidonMerkleTree, MNT4753MHTPoseidonParameters, MNT6753MHTPoseidonParameters,
    },
    FieldBasedMerkleTree
};

fn batch_poseidon_mht_eval_mnt4_full(c: &mut Criterion) {

    type MNT4PoseidonMHT = PoseidonMerkleTree<MNT4753Fr, MNT4753MHTPoseidonParameters, MNT4753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(num_leaves);

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

    type MNT6PoseidonMHT = PoseidonMerkleTree<MNT6753Fr, MNT6753MHTPoseidonParameters, MNT6753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(num_leaves);

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

    type MNT4PoseidonMHT = PoseidonMerkleTree<MNT4753Fr, MNT4753MHTPoseidonParameters, MNT4753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(num_leaves);

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

    type MNT6PoseidonMHT = PoseidonMerkleTree<MNT6753Fr, MNT6753MHTPoseidonParameters, MNT6753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(num_leaves);

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

    type MNT4PoseidonMHT = PoseidonMerkleTree<MNT4753Fr, MNT4753MHTPoseidonParameters, MNT4753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(num_leaves);

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

    type MNT6PoseidonMHT = PoseidonMerkleTree<MNT6753Fr, MNT6753MHTPoseidonParameters, MNT6753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(num_leaves);

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

    type MNT4PoseidonMHT = PoseidonMerkleTree<MNT4753Fr, MNT4753MHTPoseidonParameters, MNT4753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(num_leaves);

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

    type MNT6PoseidonMHT = PoseidonMerkleTree<MNT6753Fr, MNT6753MHTPoseidonParameters, MNT6753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(num_leaves);

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

    type MNT4PoseidonMHT = PoseidonMerkleTree<MNT4753Fr, MNT4753MHTPoseidonParameters, MNT4753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(num_leaves);

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

    type MNT6PoseidonMHT = PoseidonMerkleTree<MNT6753Fr, MNT6753MHTPoseidonParameters, MNT6753PoseidonParameters>;

    let num_leaves = 2usize.pow(10);

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(num_leaves);

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

