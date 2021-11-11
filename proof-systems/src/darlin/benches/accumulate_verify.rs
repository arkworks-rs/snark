use algebra::{AffineCurve, ToConstraintField};
use blake2::Blake2s;
use criterion::*;
use digest::Digest;
use poly_commit::{ipa_pc::InnerProductArgPC, PolynomialCommitment};
use proof_systems::darlin::pcd::GeneralPCD;
use proof_systems::darlin::{
    proof_aggregator::{accumulate_proofs, verify_aggregated_proofs},
    tests::{final_darlin::generate_test_data as generate_final_darlin_test_data, get_keys},
};
use rand::{thread_rng, SeedableRng};
use rand_xorshift::XorShiftRng;

fn bench_verify<G1: AffineCurve, G2: AffineCurve, D: Digest>(
    c: &mut Criterion,
    bench_name: &str,
    segment_size: usize,
    max_proofs: Vec<usize>,
) where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>
        + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>
        + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let rng = &mut XorShiftRng::seed_from_u64(1234567890u64);
    let mut group = c.benchmark_group(bench_name);
    let num_constraints = 1 << 19;

    //Generate DLOG keys
    let params_g1 = InnerProductArgPC::<G1, D>::setup(segment_size - 1).unwrap();
    let params_g2 = InnerProductArgPC::<G2, D>::setup(segment_size - 1).unwrap();

    let (committer_key_g1, verifier_key_g1, committer_key_g2, verifier_key_g2) =
        get_keys::<_, _, D>(&params_g1, &params_g2);

    let (final_darlin_pcd, index_vk) = generate_final_darlin_test_data::<G1, G2, D, _>(
        num_constraints,
        segment_size,
        &params_g1,
        &params_g2,
        1,
        rng,
    );

    // Generate proofs and bench
    for num_proofs in max_proofs.into_iter() {
        // Collect PCDs and vks
        let pcds = vec![GeneralPCD::FinalDarlin(final_darlin_pcd[0].clone()); num_proofs];
        let vks = vec![index_vk[0].clone(); num_proofs];

        // Accumulate PCDs
        let (proof_g1, proof_g2) = accumulate_proofs::<G1, G2, D>(
            pcds.as_slice(),
            vks.as_slice(),
            &committer_key_g1,
            &committer_key_g2,
        )
        .unwrap();

        group.bench_with_input(
            BenchmarkId::from_parameter(num_proofs),
            &num_proofs,
            |bn, _num_proofs| {
                bn.iter(|| {
                    // Verify accumulation
                    assert!(verify_aggregated_proofs::<G1, G2, D, _>(
                        pcds.as_slice(),
                        vks.as_slice(),
                        &proof_g1,
                        &proof_g2,
                        &verifier_key_g1,
                        &verifier_key_g2,
                        &mut thread_rng(),
                    )
                    .unwrap())
                });
            },
        );
    }
    group.finish();
}

fn bench_accumulate<G1: AffineCurve, G2: AffineCurve, D: Digest>(
    c: &mut Criterion,
    bench_name: &str,
    segment_size: usize,
    max_proofs: Vec<usize>,
) where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>
        + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>
        + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let rng = &mut XorShiftRng::seed_from_u64(1234567890u64);
    let mut group = c.benchmark_group(bench_name);
    let num_constraints = 1 << 19;

    //Generate DLOG keys
    let params_g1 = InnerProductArgPC::<G1, D>::setup(segment_size - 1).unwrap();
    let params_g2 = InnerProductArgPC::<G2, D>::setup(segment_size - 1).unwrap();

    let (committer_key_g1, _, committer_key_g2, _) = get_keys::<_, _, D>(&params_g1, &params_g2);

    let (final_darlin_pcd, index_vk) = generate_final_darlin_test_data::<G1, G2, D, _>(
        num_constraints,
        segment_size,
        &params_g1,
        &params_g2,
        1,
        rng,
    );

    // Generate proofs and bench
    for num_proofs in max_proofs.into_iter() {
        // Collect PCDs and vks
        let pcds = vec![GeneralPCD::FinalDarlin(final_darlin_pcd[0].clone()); num_proofs];
        let vks = vec![index_vk[0].clone(); num_proofs];

        group.bench_with_input(
            BenchmarkId::from_parameter(num_proofs),
            &num_proofs,
            |bn, _num_proofs| {
                bn.iter(|| {
                    accumulate_proofs::<G1, G2, D>(
                        pcds.as_slice(),
                        vks.as_slice(),
                        &committer_key_g1,
                        &committer_key_g2,
                    )
                    .unwrap()
                });
            },
        );
    }
    group.finish();
}

// the maximum degree we expect to handle is 2^19, maybe even below, e.g. 2^18
// Segment size |H| => 42, segment size |H|/2 => 84

fn bench_verify_tweedle(c: &mut Criterion) {
    use algebra::curves::tweedle::{dee::Affine as TweedleDee, dum::Affine as TweedleDum};

    bench_verify::<TweedleDee, TweedleDum, Blake2s>(
        c,
        "tweedle-dee, |H| = segment_size = 1 << 19, proofs",
        1 << 19,
        vec![10, 50, 100, 200],
    );

    bench_verify::<TweedleDee, TweedleDum, Blake2s>(
        c,
        "tweedle-dee, |H| = 1 << 19, segment_size = |H|/2, proofs",
        1 << 18,
        vec![10, 50, 100, 200],
    );

    bench_verify::<TweedleDee, TweedleDum, Blake2s>(
        c,
        "tweedle-dee, |H| = 1 << 19, segment_size = |H|/4, proofs",
        1 << 17,
        vec![10, 50, 100, 200],
    );
}

fn bench_accumulate_tweedle(c: &mut Criterion) {
    use algebra::curves::tweedle::{dee::Affine as TweedleDee, dum::Affine as TweedleDum};

    bench_accumulate::<TweedleDee, TweedleDum, Blake2s>(
        c,
        "tweedle-dee, |H| = segment_size = 1 << 19, proofs",
        1 << 19,
        vec![10, 50, 100, 200],
    );

    bench_accumulate::<TweedleDee, TweedleDum, Blake2s>(
        c,
        "tweedle-dee, |H| = 1 << 19, segment_size = |H|/2, proofs",
        1 << 18,
        vec![10, 50, 100, 200],
    );

    bench_accumulate::<TweedleDee, TweedleDum, Blake2s>(
        c,
        "tweedle-dee, |H| = 1 << 19, segment_size = |H|/4, proofs",
        1 << 17,
        vec![10, 50, 100, 200],
    );
}

criterion_group!(
name = accumulate_verify;
config = Criterion::default().sample_size(10);
targets = bench_verify_tweedle, bench_accumulate_tweedle
);

criterion_main!(accumulate_verify);
