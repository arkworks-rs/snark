use algebra::{
    AffineCurve, ToConstraintField
};
use marlin::VerifierKey as MarlinVerifierKey;
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
        final_darlin::{FinalDarlinPCD, FinalDarlinPCDVerifierKey}
    }
};
use rand::{RngCore, thread_rng};
use digest::Digest;
use rayon::prelude::*;

pub mod pcd;
pub mod accumulators;

#[cfg(test)]
mod tests;

//TODO: Add tracing
//TODO: Remove dependency from R: RngCore when not needed
//TODO: Do the same with Digest template
pub fn accumulate_proofs<G1, G2, D: Digest>(
    final_darlin_pcds:      &[FinalDarlinPCD<G1, G2, D>],
    final_darlin_vks:       &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    marlin_pcds:            &[SimpleMarlinPCD<G1, D>],
    marlin_vks:             &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    g1_ck:                  &DLogCommitterKey<G1>,
    g2_ck:                  &DLogCommitterKey<G2>,
) -> Result<
    (
        Option<AccumulationProof<G1>>,
        Option<AccumulationProof<G2>>,
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
    let acc_proof_g1 = if accs_g1.is_empty() {
        None
    } else {
        Some(DLogAccumulator::<G1>::accumulate::<D>(g1_ck, accs_g1)?.1)
    };

    let acc_proof_g2 = if accs_g2.is_empty() {
        None
    } else {
        Some(DLogAccumulator::<G2>::accumulate::<D>(g2_ck, accs_g2)?.1)
    };

    Ok((acc_proof_g1, acc_proof_g2))
}

pub fn verify_aggregated_proofs<G1, G2, D: Digest, R: RngCore>(
    final_darlin_pcds:      &[FinalDarlinPCD<G1, G2, D>],
    final_darlin_vks:       &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    marlin_pcds:            &[SimpleMarlinPCD<G1, D>],
    marlin_vks:             &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    accumulation_proof_g1:  Option<&AccumulationProof<G1>>,
    accumulation_proof_g2:  Option<&AccumulationProof<G2>>,
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

    // Verify accumulators and accumulation proofs
    let result_accumulate_g1 = if accumulation_proof_g1.is_some() {
        let dummy_g1 = DLogAccumulator::<G1>::default();
        dummy_g1.verify_accumulate::<R, D>(g1_vk, accs_g1, accumulation_proof_g1.unwrap(), rng)?
    } else {
        true
    };

    let result_accumulate_g2 = if accumulation_proof_g2.is_some() {
        let dummy_g2 = DLogAccumulator::<G2>::default();
        dummy_g2.verify_accumulate::<R, D>(g2_vk, accs_g2, accumulation_proof_g2.unwrap(), rng)?
    } else {
        true
    };

    Ok(result_accumulate_g1 && result_accumulate_g2)
}