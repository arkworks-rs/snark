//! Test suite for PCD post processing (batch-verification, aggregation)
use algebra::AffineCurve;
use poly_commit::{
    PCUniversalParams, PolynomialCommitment,
    ipa_pc::{
        InnerProductArgPC,
        UniversalParams,
        CommitterKey as DLogCommitterKey, VerifierKey as DLogVerifierKey,
    }
};
use digest::Digest;

pub mod simple_marlin;
pub mod final_darlin;

#[allow(dead_code)]
/// Extract DLogCommitterKey and DLogVerifierKey from UniversalParams struct
pub fn get_keys<G1: AffineCurve, G2: AffineCurve, D: Digest>(
    params_g1: &UniversalParams<G1>,
    params_g2: &UniversalParams<G2>,
) -> (DLogCommitterKey<G1>, DLogVerifierKey<G1>, DLogCommitterKey<G2>, DLogVerifierKey<G2>)
{
    let (ck_g1, vk_g1) = InnerProductArgPC::<G1, D>::trim(
        params_g1,
        params_g1.max_degree(),
    ).unwrap();

    let (ck_g2, vk_g2) = InnerProductArgPC::<G2, D>::trim(
        params_g2,
        params_g2.max_degree(),
    ).unwrap();

    (ck_g1, vk_g1, ck_g2, vk_g2)
}

#[cfg(test)]
mod test {
    use super::*;
    use algebra::{curves::tweedle::{
        dee::Affine as DeeAffine, dum::Affine as DumAffine,
    }, UniformRand, ToConstraintField, serialize::test_canonical_serialize_deserialize, SemanticallyValid, CanonicalSerialize, CanonicalDeserialize};
    use marlin::VerifierKey as MarlinVerifierKey;
    use crate::darlin::{
        pcd::GeneralPCD,
        proof_aggregator::{accumulate_proofs, verify_aggregated_proofs, batch_verify_proofs},
        tests::{
            simple_marlin::generate_test_data as generate_simple_marlin_test_data,
            final_darlin::generate_test_data as generate_final_darlin_test_data,
        }
    };
    use blake2::Blake2s;
    use rand::{Rng, RngCore, SeedableRng, thread_rng};
    use rand_xorshift::XorShiftRng;
    use std::path::Path;
    use std::fs::File;
    use crate::darlin::data_structures::FinalDarlinProof;

    /// Generic test for `accumulate_proofs` and `verify_aggregated_proofs`
    fn test_accumulation<G1: AffineCurve, G2: AffineCurve, D: Digest, R: RngCore>(
        pcds: &mut [GeneralPCD<G1, G2, D>],
        vks: &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
        committer_key_g1: &DLogCommitterKey<G1>,
        committer_key_g2: &DLogCommitterKey<G2>,
        verifier_key_g1: &DLogVerifierKey<G1>,
        verifier_key_g2: &DLogVerifierKey<G2>,
        rng: &mut R
    )
        where
            G1: AffineCurve<BaseField=<G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
            G2: AffineCurve<BaseField=<G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
    {
        // Accumulate PCDs
        let (proof_g1, proof_g2) = accumulate_proofs::<G1, G2, D>(
            pcds,
            vks,
            committer_key_g1,
            committer_key_g2
        ).unwrap();

        // Verify accumulation
        assert!(verify_aggregated_proofs::<G1, G2, D, R>(
            pcds,
            vks,
            &proof_g1,
            &proof_g2,
            verifier_key_g1,
            verifier_key_g2,
            rng
        ).unwrap());

        // Pass wrong accumulation proof and check verification fails
        // Change one element in proof_g1
        let mut wrong_proof_g1 = proof_g1.clone().unwrap();
        wrong_proof_g1.pc_proof.c = G1::ScalarField::rand(rng);
        assert!(!verify_aggregated_proofs::<G1, G2, D, R>(
            pcds,
            vks,
            &Some(wrong_proof_g1),
            &proof_g2,
            verifier_key_g1,
            verifier_key_g2,
            rng
        ).unwrap());

        // Randomize usr_ins for one marlin PCD and assert AHP verification fails
        let idx: usize = rng.gen_range(0, pcds.len());
        let original_pcd = pcds[idx].clone(); // Save correct pcd
        pcds[idx].randomize_usr_ins(rng);

        let result = verify_aggregated_proofs::<G1, G2, D, R>(
            pcds,
            vks,
            &proof_g1,
            &proof_g2,
            verifier_key_g1,
            verifier_key_g2,
            rng
        );

        // Check AHP failed
        assert!(result.is_err());

        // Since the AHP failed, we are able to determine which proof verification has failed
        assert_eq!(result.unwrap_err().unwrap(), idx);

        // Restore correct PCD
        pcds[idx] = original_pcd;

        // Randomize sys_ins for one marlin PCD and assert AHP verification fails
        let idx: usize = rng.gen_range(0, pcds.len());
        pcds[idx].randomize_sys_ins(committer_key_g1, committer_key_g2, rng);

        let result = verify_aggregated_proofs::<G1, G2, D, R>(
            pcds,
            vks,
            &proof_g1,
            &proof_g2,
            verifier_key_g1,
            verifier_key_g2,
            rng
        );

        // Check AHP failed
        assert!(result.is_err());

        // Since the AHP failed, we are able to determine which proof verification has failed
        assert_eq!(result.unwrap_err().unwrap(), idx);
    }

    /// Generic test for `batch_verify_proofs`
    fn test_batch_verification<G1: AffineCurve, G2: AffineCurve, D: Digest, R: RngCore>(
        pcds: &mut [GeneralPCD<G1, G2, D>],
        vks: &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
        verifier_key_g1: &DLogVerifierKey<G1>,
        verifier_key_g2: &DLogVerifierKey<G2>,
        rng: &mut R
    )
        where
            G1: AffineCurve<BaseField=<G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
            G2: AffineCurve<BaseField=<G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
    {
        // Batch Verify
        assert!(batch_verify_proofs::<G1, G2, D, R>(
            pcds,
            vks,
            verifier_key_g1,
            verifier_key_g2,
            rng
        ).unwrap());

        // Pass wrong public inputs for one marlin PCD and check AHP verification fails
        // Select randomly one pcds to which changing the input
        let idx: usize = rng.gen_range(0, pcds.len());
        let original_pcd = pcds[idx].clone(); // Save correct pcd
        pcds[idx].randomize_usr_ins(rng);

        let result = batch_verify_proofs::<G1, G2, D, R>(
            pcds,
            vks,
            verifier_key_g1,
            verifier_key_g2,
            rng
        );

        // Check AHP failed
        assert!(result.is_err());

        // Since the AHP failed, we are able to determine which proof verification has failed
        assert_eq!(result.unwrap_err().unwrap(), idx);

        // Restore correct PCD
        pcds[idx] = original_pcd;

        // Randomize sys_ins for one marlin PCD and assert AHP verification fails
        let idx: usize = rng.gen_range(0, pcds.len());
        pcds[idx].randomize_sys_ins(verifier_key_g1, verifier_key_g2, rng);

        let result = batch_verify_proofs::<G1, G2, D, R>(
            pcds,
            vks,
            verifier_key_g1,
            verifier_key_g2,
            rng
        );

        // Check AHP failed
        assert!(result.is_err());

        // Since the AHP failed, we are able to determine which proof verification has failed
        assert_eq!(result.unwrap_err().unwrap(), idx);
    }

    type TestIPAPCDee = InnerProductArgPC<DeeAffine, Blake2s>;
    type TestIPAPCDum = InnerProductArgPC<DumAffine, Blake2s>;

    #[test]
    fn test_simple_marlin_proof_aggregator() {
        let rng = &mut XorShiftRng::seed_from_u64(1234567890u64);

        // Set params
        let max_proofs = 100;
        let max_pow = 7usize;
        let num_constraints = 1 << max_pow;
        let segment_size = num_constraints;

        //Generate keys
        let params_g1 = TestIPAPCDee::setup(segment_size - 1).unwrap();
        let params_g2 = TestIPAPCDum::setup(segment_size - 1).unwrap();

        let (
            committer_key_g1, verifier_key_g1,
            committer_key_g2, verifier_key_g2
        ) = get_keys::<_, _, Blake2s>(&params_g1, &params_g2);

        test_canonical_serialize_deserialize(true, &committer_key_g1);
        test_canonical_serialize_deserialize(true, &committer_key_g2);
        test_canonical_serialize_deserialize(true, &verifier_key_g1);
        test_canonical_serialize_deserialize(true, &verifier_key_g2);

        // Generate pcds and index vks: we want to generate PCDs with different segment
        // sizes up at least to max_proofs, so we are going to randomly sample the segment size
        // and the number of proofs with that specific segment size to generate (to save
        // time without lacking in expressivity we just generate one and clone it
        // iteration_num_proofs time.
        let mut generated_proofs = 0;
        let mut pcds = Vec::new();
        let mut simple_marlin_vks = Vec::new();
        let generation_rng = &mut thread_rng();
        while generated_proofs < max_proofs {
            let iteration_num_proofs: usize = generation_rng.gen_range(1, max_proofs);
            generated_proofs += iteration_num_proofs;
            let iteration_segment_size = 1 << (generation_rng.gen_range(2, max_pow));
            let (mut iteration_pcds, mut iteration_vks) = generate_simple_marlin_test_data(
                num_constraints - 1,
                iteration_segment_size,
                &params_g1,
                iteration_num_proofs,
                generation_rng
            );

            assert!(&iteration_pcds[0].proof.is_valid());
            test_canonical_serialize_deserialize(true, &iteration_pcds[0].proof);
            test_canonical_serialize_deserialize(true, &iteration_vks[0]);

            pcds.append(&mut iteration_pcds);
            simple_marlin_vks.append(&mut iteration_vks);
        }

        // Collect PCDs
        let mut simple_marlin_pcds = pcds
            .into_iter()
            .map(|simple_marlin_pcd| GeneralPCD::SimpleMarlin(simple_marlin_pcd))
            .collect::<Vec<_>>();

        println!("Test accumulation");
        test_accumulation::<DeeAffine, DumAffine, Blake2s, _>(
            simple_marlin_pcds.clone().as_mut_slice(),
            simple_marlin_vks.as_slice(),
            &committer_key_g1,
            &committer_key_g2,
            &verifier_key_g1,
            &verifier_key_g2,
            rng
        );

        println!("Test batch verification");
        test_batch_verification::<DeeAffine, DumAffine, Blake2s, _>(
            simple_marlin_pcds.as_mut_slice(),
            simple_marlin_vks.as_slice(),
            &verifier_key_g1,
            &verifier_key_g2,
            rng
        );
    }

    #[test]
    fn test_final_darlin_proof_aggregator() {
        let rng = &mut XorShiftRng::seed_from_u64(1234567890u64);

        // Set params
        let max_proofs = 100;
        let max_pow = 7usize;
        let num_constraints = 1 << max_pow;
        let segment_size = num_constraints;

        //Generate keys
        let params_g1 = TestIPAPCDee::setup(segment_size - 1).unwrap();
        let params_g2 = TestIPAPCDum::setup(segment_size - 1).unwrap();

        let (
            committer_key_g1, verifier_key_g1,
            committer_key_g2, verifier_key_g2
        ) = get_keys::<_, _, Blake2s>(&params_g1, &params_g2);

        test_canonical_serialize_deserialize(true, &committer_key_g1);
        test_canonical_serialize_deserialize(true, &committer_key_g2);
        test_canonical_serialize_deserialize(true, &verifier_key_g1);
        test_canonical_serialize_deserialize(true, &verifier_key_g2);

        // Generate pcds and index vks: we want to generate PCDs with different segment
        // sizes up to at least max_proofs, so we are going to randomly sample the segment size
        // and the number of proofs with that specific segment size to generate (to save
        // time without lacking in expressivity we just generate one and clone it
        // iteration_num_proofs time.
        let mut generated_proofs = 0;
        let mut pcds = Vec::new();
        let mut final_darlin_vks = Vec::new();
        let generation_rng = &mut thread_rng();
        while generated_proofs < max_proofs {
            let iteration_num_proofs: usize = generation_rng.gen_range(1, max_proofs);
            generated_proofs += iteration_num_proofs;
            let iteration_segment_size = 1 << (generation_rng.gen_range(2, max_pow));
            let (mut iteration_pcds, mut iteration_vks) = generate_final_darlin_test_data(
                num_constraints - 1,
                iteration_segment_size,
                &params_g1,
                &params_g2,
                iteration_num_proofs,
                generation_rng
            );

            assert!(&iteration_pcds[0].final_darlin_proof.is_valid());
            test_canonical_serialize_deserialize(true, &iteration_pcds[0].final_darlin_proof);
            test_canonical_serialize_deserialize(true, &iteration_vks[0]);

            pcds.append(&mut iteration_pcds);
            final_darlin_vks.append(&mut iteration_vks);
        }

        // Collect PCDs
        let mut final_darlin_pcds = pcds
            .into_iter()
            .map(|final_darlin_pcd| GeneralPCD::FinalDarlin(final_darlin_pcd))
            .collect::<Vec<_>>();

        println!("Test accumulation");
        test_accumulation::<DeeAffine, DumAffine, Blake2s, _>(
            final_darlin_pcds.clone().as_mut_slice(),
            final_darlin_vks.as_slice(),
            &committer_key_g1,
            &committer_key_g2,
            &verifier_key_g1,
            &verifier_key_g2,
            rng
        );

        println!("Test batch verification");
        test_batch_verification::<DeeAffine, DumAffine, Blake2s, _>(
            final_darlin_pcds.as_mut_slice(),
            final_darlin_vks.as_slice(),
            &verifier_key_g1,
            &verifier_key_g2,
            rng
        );
    }

    #[test]
    fn test_mixed_proof_aggregator() {
        let rng = &mut XorShiftRng::seed_from_u64(1234567890u64);

        // Set params
        let max_proofs = 100;
        let max_pow = 7usize;
        let num_constraints = 1 << max_pow;
        let segment_size = num_constraints;

        //Generate keys
        let params_g1 = TestIPAPCDee::setup(segment_size - 1).unwrap();
        let params_g2 = TestIPAPCDum::setup(segment_size - 1).unwrap();

        let (
            committer_key_g1, verifier_key_g1,
            committer_key_g2, verifier_key_g2
        ) = get_keys::<_, _, Blake2s>(&params_g1, &params_g2);

        test_canonical_serialize_deserialize(true, &committer_key_g1);
        test_canonical_serialize_deserialize(true, &committer_key_g2);
        test_canonical_serialize_deserialize(true, &verifier_key_g1);
        test_canonical_serialize_deserialize(true, &verifier_key_g2);

        // Generate pcds and index vks: we want to generate PCDs with different segment
        // sizes up to at least max_proofs, so we are going to randomly sample the segment size
        // and the number of proofs with that specific segment size to generate (to save
        // time without lacking in expressivity we just generate one and clone it
        // iteration_num_proofs time.
        let generation_rng = &mut thread_rng();
        let mut generated_proofs = 0;
        let mut pcds = Vec::new();
        let mut vks = Vec::new();
        while generated_proofs < max_proofs {
            let iteration_num_proofs: usize = generation_rng.gen_range(1, max_proofs);
            generated_proofs += iteration_num_proofs;
            let iteration_segment_size = 1 << (generation_rng.gen_range(2, max_pow));

            // Randomly choose if to generate a SimpleMarlinProof or a FinalDarlinProof
            let simple: bool = generation_rng.gen();
            if simple {
                let (iteration_pcds, mut iteration_vks) = generate_simple_marlin_test_data(
                    num_constraints - 1,
                    iteration_segment_size,
                    &params_g1,
                    iteration_num_proofs,
                    generation_rng
                );

                assert!(&iteration_pcds[0].proof.is_valid());
                test_canonical_serialize_deserialize(true, &iteration_pcds[0].proof);
                test_canonical_serialize_deserialize(true, &iteration_vks[0]);

                let mut iteration_pcds = iteration_pcds.into_iter().map(|pcd| GeneralPCD::SimpleMarlin(pcd)).collect::<Vec<_>>();

                pcds.append(&mut iteration_pcds);
                vks.append(&mut iteration_vks);
            } else {
                let (iteration_pcds, mut iteration_vks) = generate_final_darlin_test_data(
                    num_constraints - 1,
                    iteration_segment_size,
                    &params_g1,
                    &params_g2,
                    iteration_num_proofs,
                    generation_rng
                );

                assert!(&iteration_pcds[0].final_darlin_proof.is_valid());
                test_canonical_serialize_deserialize(true, &iteration_pcds[0].final_darlin_proof);
                test_canonical_serialize_deserialize(true, &iteration_vks[0]);

                let mut iteration_pcds = iteration_pcds.into_iter().map(|pcd| GeneralPCD::FinalDarlin(pcd)).collect::<Vec<_>>();

                pcds.append(&mut iteration_pcds);
                vks.append(&mut iteration_vks);
            }
        }

        println!("Test accumulation");
        test_accumulation::<DeeAffine, DumAffine, Blake2s, _>(
            pcds.clone().as_mut_slice(),
            vks.as_slice(),
            &committer_key_g1,
            &committer_key_g2,
            &verifier_key_g1,
            &verifier_key_g2,
            rng
        );

        println!("Test batch verification");
        test_batch_verification::<DeeAffine, DumAffine, Blake2s, _>(
            pcds.as_mut_slice(),
            vks.as_slice(),
            &verifier_key_g1,
            &verifier_key_g2,
            rng
        );
    }

    #[ignore]
    #[test]
    fn test_final_darlin_size() {

        // Set params
        let num_constraints = 1 << 19;
        let segment_size = 1 << 17;

        //Generate keys
        let params_g1 = TestIPAPCDee::setup(segment_size - 1).unwrap();
        let params_g2 = TestIPAPCDum::setup(segment_size - 1).unwrap();

        let generation_rng = &mut thread_rng();

        let file_path = "./size_test_proof";
        let proof;

        if Path::new(file_path).exists() {

            let fs = File::open(file_path).unwrap();
            proof = FinalDarlinProof::deserialize(fs).unwrap();

        } else {

            let (iteration_pcds, _) = generate_final_darlin_test_data::<_, _, Blake2s, _>(
                num_constraints - 1,
                segment_size,
                &params_g1,
                &params_g2,
                1,
                generation_rng
            );

            proof = iteration_pcds[0].final_darlin_proof.clone();

            let fs = File::create(file_path).unwrap();
            proof.serialize(fs).unwrap();
            let _ = std::fs::remove_file(file_path);
        }

        test_canonical_serialize_deserialize(true, &proof);

        println!("{} - FinalDarlinProof", proof.serialized_size());
        println!("-- {} - MarlinProof", proof.proof.serialized_size());
        println!("---- {} - commitments ({})",
                 proof.proof.commitments.serialized_size() -
                    (proof.proof.commitments.iter().flatten().collect::<Vec<_>>().len() * 4),
                 proof.proof.commitments.iter().flatten().collect::<Vec<_>>().len()
        );
        println!("---- {} - evaluations ({})",
                 proof.proof.evaluations.serialized_size() - 8,
                 proof.proof.evaluations.len()
        );
        println!("---- {} - pc_proof", proof.proof.pc_proof.serialized_size());
        println!("-- {} - FinalDarlinDeferredData", proof.deferred.serialized_size());
        println!("---- {} - DLogAccumulatorG1", proof.deferred.pre_previous_acc.serialized_size());
        println!("------ {} - G_final", proof.deferred.pre_previous_acc.g_final.serialized_size());
        println!("------ {} - xi_s", proof.deferred.pre_previous_acc.xi_s.serialized_size());
        println!("---- {} - DLogAccumulatorG2", proof.deferred.previous_acc.serialized_size());
        println!("------ {} - G_final", proof.deferred.previous_acc.g_final.serialized_size());
        println!("------ {} - xi_s", proof.deferred.previous_acc.xi_s.serialized_size());
    }
}