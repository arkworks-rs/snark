#[macro_use]
extern crate criterion;

use criterion::Criterion;
use primitives::{MNT4PoseidonHash, MNT6PoseidonHash, FieldBasedHash, BatchFieldBasedHash};
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::UniformRand;
use rand_xorshift::XorShiftRng;
use rand::SeedableRng;
use primitives::crh::poseidon::parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters};
use primitives::crh::poseidon::batched_crh::PoseidonBatchHash;
use primitives::merkle_tree::field_based_mht::batch_mht::poseidon::{PoseidonBatchMerkleTree, PoseidonBatchMerkleTreeMem};
use primitives::merkle_tree::field_based_mht::batch_mht::BatchMerkleTree;

fn batch_poseidon_mht_eval_mnt4(c: &mut Criterion) {

    type MNT4BatchedMerkleTree = PoseidonBatchMerkleTree<MNT4753Fr, MNT4753PoseidonParameters>;

    let num_leaves = 64;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4BatchedMerkleTree::new(num_leaves, 64);

    c.bench_function("Batch Poseidon MHT Eval for MNT4 (64 leaves)", move |b| {
        b.iter(|| {
            tree.reset();
            for _ in 0..num_leaves {
                tree.update(MNT4753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6(c: &mut Criterion) {

    type MNT6BatchedMerkleTree = PoseidonBatchMerkleTree<MNT6753Fr, MNT6753PoseidonParameters>;

    let num_leaves = 64;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6BatchedMerkleTree::new(num_leaves, 64);

    c.bench_function("Batch Poseidon MHT Eval for MNT6 (64 leaves)", move |b| {
        b.iter(|| {
            tree.reset();
            for _ in 0..num_leaves {
                tree.update(MNT6753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
        })
    });
}

fn batch_poseidon_mht_eval_mnt4_mem(c: &mut Criterion) {

    type MNT4BatchedMerkleTree = PoseidonBatchMerkleTreeMem<MNT4753Fr, MNT4753PoseidonParameters>;

    let num_leaves = 64;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4BatchedMerkleTree::new(num_leaves, 64);

    c.bench_function("Batch Poseidon MHT Eval for MNT4 opt mem (64 leaves)", move |b| {
        b.iter(|| {
            tree.reset();
            for _ in 0..num_leaves {
                tree.update(MNT4753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6_mem(c: &mut Criterion) {

    type MNT6BatchedMerkleTree = PoseidonBatchMerkleTreeMem<MNT6753Fr, MNT6753PoseidonParameters>;

    let num_leaves = 64;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6BatchedMerkleTree::new(num_leaves, 64);

    c.bench_function("Batch Poseidon MHT Eval for MNT6 opt mem (64 leaves)", move |b| {
        b.iter(|| {
            tree.reset();
            for _ in 0..num_leaves {
                tree.update(MNT6753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
        })
    });
}


criterion_group! {
    name = mht_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = batch_poseidon_mht_eval_mnt4, batch_poseidon_mht_eval_mnt6, batch_poseidon_mht_eval_mnt4_mem, batch_poseidon_mht_eval_mnt6_mem
}

criterion_main! (
    mht_poseidon_eval
);

