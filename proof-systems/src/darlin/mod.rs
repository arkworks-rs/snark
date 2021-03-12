use algebra::{
    AffineCurve, ToConstraintField
};
use marlin::{MarlinConfig, VerifierKey as MarlinVerifierKey};
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC,
        CommitterKey as DLogCommitterKey, VerifierKey as DLogVerifierKey,
    },
    Error as PCError
};
use crate::darlin::{
    accumulators::{
        dlog::DLogAccumulator,
        Accumulator, AccumulationProof
    },
    pcd::{
        PCD,
        simple_marlin::{SimpleMarlinPCD, SimpleMarlinPCDVerifierKey},
        finalized_darlin::{FinalDarlinPCD, FinalDarlinPCDVerifierKey}
    }
};
use rand::{RngCore, thread_rng};
use digest::Digest;
use rayon::prelude::*;

pub mod pcd;
pub mod accumulators;

//TODO: Remove dependency from R: RngCore when not needed

//TODO: Get rid of this MarlinConfig template as it obliges all the proofs to have the same
//      MarlinConfig. Either we keep it inside the proof or the PCD (by converting it to a struct)
//      or we remove it (regarding LC_OPT, we never use it; regarding ZK, only the prover needs to
//      know about it; regarding SegmentSize it can be done in a transparent way)

//TODO: Do the same with Digest template
pub fn accumulate_proofs<G1, G2, D: Digest, R: RngCore, MC: MarlinConfig>(
    final_darlin_pcds:      &[FinalDarlinPCD<G1, G2, D, MC>],
    final_darlin_vks:       &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    marlin_pcds:            &[SimpleMarlinPCD<G1, D, MC>],
    marlin_vks:             &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    g1_ck:                  &DLogCommitterKey<G1>,
    g2_ck:                  &DLogCommitterKey<G2>,
    rng:                    &mut R
) -> Result<
    (
        AccumulationProof<G1>,
        AccumulationProof<G2>,
    ), PCError>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    // Final Darlin
    let darlin_accs = final_darlin_pcds
        .par_iter()
        .zip(final_darlin_vks)
        .map(|(final_darlin_pcd, final_darlin_vk)|
            {
                let vk = FinalDarlinPCDVerifierKey::<G1, G2, D>{
                    marlin_vk: final_darlin_vk,
                    dlog_vks: (g1_ck, g2_ck)
                };
                final_darlin_pcd.succinct_verify(&vk, &mut thread_rng())
            }
        ).collect::<Result<Vec<_>, PCError>>()?;

    let mut accs_g1 = darlin_accs.iter().flat_map(|acc| acc.0.clone()).collect::<Vec<_>>();
    let accs_g2 = darlin_accs.into_iter().flat_map(|acc| acc.1).collect::<Vec<_>>();

    // Marlin
    let mut marlin_accs_g1 = marlin_pcds
        .into_par_iter()
        .zip(marlin_vks)
        .map(|(marlin_pcd, marlin_vk)|
            {
                let vk = SimpleMarlinPCDVerifierKey::<G1, D>(marlin_vk, g1_ck);
                marlin_pcd.succinct_verify(&vk, &mut thread_rng())
            }
        ).collect::<Result<Vec<_>, PCError>>()?;

    accs_g1.append(&mut marlin_accs_g1);

    // Create accumulation proofs
    let (_, acc_proof_g1) = DLogAccumulator::<G1>::accumulate::<R, D>(g1_ck, accs_g1, rng)?;
    let (_, acc_proof_g2) = DLogAccumulator::<G2>::accumulate::<R, D>(g2_ck, accs_g2, rng)?;

    Ok((acc_proof_g1, acc_proof_g2))
}

pub fn verify_aggregated_proofs<G1, G2, D: Digest, R: RngCore, MC: MarlinConfig>(
    final_darlin_pcds:      &[FinalDarlinPCD<G1, G2, D, MC>],
    final_darlin_vks:       &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    marlin_pcds:            &[SimpleMarlinPCD<G1, D, MC>],
    marlin_vks:             &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    accumulation_proof_g1:  &AccumulationProof<G1>,
    accumulation_proof_g2:  &AccumulationProof<G2>,
    g1_vk:                  &DLogVerifierKey<G1>,
    g2_vk:                  &DLogVerifierKey<G2>,
    rng:                    &mut R
) -> Result<bool, PCError>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    // Final Darlin
    let darlin_accs = final_darlin_pcds
        .into_par_iter()
        .zip(final_darlin_vks)
        .map(|(final_darlin_pcd, final_darlin_vk)|
            {
                let vk = FinalDarlinPCDVerifierKey::<G1, G2, D>{
                    marlin_vk: final_darlin_vk,
                    dlog_vks: (g1_vk, g2_vk)
                };
                final_darlin_pcd.succinct_verify(&vk, &mut thread_rng())
            }
        ).collect::<Result<Vec<_>, PCError>>()?;

    let mut accs_g1 = darlin_accs.iter().flat_map(|acc| acc.0.clone()).collect::<Vec<_>>();
    let accs_g2 = darlin_accs.into_iter().flat_map(|acc| acc.1).collect::<Vec<_>>();

    // Marlin
    let mut marlin_accs_g1 = marlin_pcds
        .into_par_iter()
        .zip(marlin_vks)
        .map(|(marlin_pcd, marlin_vk)|
            {
                let vk = SimpleMarlinPCDVerifierKey::<G1, D>(marlin_vk, g1_vk);
                marlin_pcd.succinct_verify(&vk, &mut thread_rng())
            }
        ).collect::<Result<Vec<_>, PCError>>()?;

    accs_g1.append(&mut marlin_accs_g1);

    // Create dummy DLOGAccumulators on which to call verify accumulate
    let dummy_g1 = DLogAccumulator::<G1>::default();
    let dummy_g2 = DLogAccumulator::<G2>::default();

    // Verify accumulators and accumulation proofs
    let result = dummy_g1.verify_accumulate::<R, D>(g1_vk, accs_g1, accumulation_proof_g1, rng)? &&
        dummy_g2.verify_accumulate::<R, D>(g2_vk, accs_g2, accumulation_proof_g2, rng)?;

    Ok(result)
}