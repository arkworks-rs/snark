#[macro_use]
extern crate criterion;

use criterion::Criterion;
use primitives::{MNT4PoseidonHash, MNT6PoseidonHash, FieldBasedHash, PoseidonBatchHash, BatchFieldBasedHash, BatchMerkleTree};
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::fields::mnt6753::Fr as MNT6753Fr;
use algebra::UniformRand;
use rand_xorshift::XorShiftRng;
use rand::SeedableRng;
use primitives::crh::poseidon::parameters::{MNT4753PoseidonParameters, MNT6753PoseidonParameters};
use primitives::merkle_tree::field_based_mht::batched_mht::BatchMerkleTree;

fn batch_poseidon_mht_eval_mnt4(c: &mut Criterion) {

    type Mnt4BatchPoseidonHash = PoseidonBatchHash<MNT4753Fr, MNT4753PoseidonParameters>;
    type MNT4BatchedMerkleTree = BatchMerkleTree<MNT4PoseidonBatchHash>;

    let num_leaves = 1024;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4BatchedMerkleTree::new(num_leaves, 1024);

    c.bench_function("Batch Poseidon MHT Eval for MNT4 (1024 leaves)", move |b| {
        b.iter(|| {
            for _ in 0..num_leaves {
                tree.push(MNT4753Fr::rand(&mut rng));
                tree.finalize();
            }
        })
    });
}


criterion_group! {
    name = mht_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = batch_poseidon_mht_eval_mnt4
}

criterion_main! (
    mht_poseidon_eval
);

