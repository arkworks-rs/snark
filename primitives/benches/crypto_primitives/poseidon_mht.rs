#[macro_use]
extern crate criterion;

use criterion::Criterion;
use algebra::fields::mnt4753::Fr as MNT4753Fr;
use algebra::fields::mnt6753::Fr as MNT6753Fr;
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

//TODO: Bench new added functions ?
fn batch_poseidon_mht_eval_mnt4(c: &mut Criterion) {

    type MNT4PoseidonMHT = PoseidonMerkleTree<MNT4753Fr, MNT4753MHTPoseidonParameters, MNT4753PoseidonParameters>;

    let num_leaves = 64;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT4PoseidonMHT::init(num_leaves);

    c.bench_function("Batch Poseidon MHT Eval for MNT4 (64 leaves)", move |b| {
        b.iter(|| {
            tree.reset();
            for _ in 0..num_leaves {
                tree.append(MNT4753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
        })
    });
}

fn batch_poseidon_mht_eval_mnt6(c: &mut Criterion) {

    type MNT6PoseidonMHT = PoseidonMerkleTree<MNT6753Fr, MNT6753MHTPoseidonParameters, MNT6753PoseidonParameters>;

    let num_leaves = 64;

    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

    let mut tree = MNT6PoseidonMHT::init(num_leaves);

    c.bench_function("Batch Poseidon MHT Eval for MNT6 (64 leaves)", move |b| {
        b.iter(|| {
            tree.reset();
            for _ in 0..num_leaves {
                tree.append(MNT6753Fr::rand(&mut rng));
            }
            tree.finalize_in_place();
        })
    });
}


criterion_group! {
    name = mht_poseidon_eval;
    config = Criterion::default().sample_size(20);
    targets = batch_poseidon_mht_eval_mnt4, batch_poseidon_mht_eval_mnt6
}

criterion_main! (
    mht_poseidon_eval
);

