//! Utilities for proof post-processing of `GeneralPCD`, i.e. SimpleMarlin and
//! FinalDarlin PCD, using batch verification and aggregation of their dlog hard parts.
use crate::darlin::{
    accumulators::{
        dlog::{DLogItem, DLogItemAccumulator},
        AccumulationProof, ItemAccumulator,
    },
    pcd::{DualPCDVerifierKey, GeneralPCD, PCD},
};
use algebra::{AffineCurve, ToConstraintField};
use digest::Digest;
use marlin::VerifierKey as MarlinVerifierKey;
use poly_commit::ipa_pc::{
    CommitterKey as DLogCommitterKey, InnerProductArgPC, VerifierKey as DLogVerifierKey,
};
use rand::RngCore;
use rayon::prelude::*;

/// Given a set of PCDs, their corresponding Marlin verification keys, and the DLogCommitterKey(s)
/// over two groups of a curve cycle, compute and return the associated accumulators via the
/// succinct verification of them.
/// In case of failure, return the indices of the proofs that have caused the failure (if it's possible
/// to establish it).
/// The PCDs are allowed to use different size restrictions of the DLogCommitterKey `g1_ck` and `g2_ck`.
pub(crate) fn get_accumulators<G1, G2, D: Digest>(
    pcds: &[GeneralPCD<G1, G2, D>],
    vks: &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    g1_ck: &DLogCommitterKey<G1>,
    g2_ck: &DLogCommitterKey<G2>,
) -> Result<(Vec<DLogItem<G1>>, Vec<DLogItem<G2>>), Option<Vec<usize>>>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>
        + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>
        + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let accumulators_time = start_timer!(|| "Compute accumulators");

    if pcds.is_empty() || vks.is_empty() || pcds.len() != vks.len() {
        return Err(None);
    }

    let (accs, failing_indices): (Vec<_>, Vec<_>) = pcds
        .into_par_iter()
        .zip(vks)
        .enumerate()
        .map(|(i, (pcd, vk))| {
            // recall that we use FinalDarlinVerifierKeys to handle
            // polymorphic verification of final Darlin/simpleM arlin PCDs
            let vk = DualPCDVerifierKey::<G1, G2, D> {
                final_darlin_vk: vk,
                dlog_vks: (g1_ck, g2_ck),
            };
            // No need to trim the vk here to the specific segment size used
            // to generate the proof for this pcd, as the IPA succinct_check
            // function doesn't use vk.comm_key at all.
            pcd.succinct_verify(&vk).map_err(|_| i)
        })
        .partition(Result::is_ok);
    end_timer!(accumulators_time);

    let accs = accs.into_iter().map(Result::unwrap).collect::<Vec<_>>();
    let mut failing_indices = failing_indices
        .into_iter()
        .map(Result::unwrap_err)
        .collect::<Vec<_>>();

    if failing_indices.is_empty() {
        // All succinct verifications passed: collect and return the accumulators
        let accs_g1 = accs
            .iter()
            .flat_map(|acc| acc.0.clone())
            .collect::<Vec<_>>();
        let accs_g2 = accs.into_iter().flat_map(|acc| acc.1).collect::<Vec<_>>();
        Ok((accs_g1, accs_g2))
    } else {
        // Otherwise, collect and return as error the indices of all the failing proofs
        // sorted in ascending order
        failing_indices.sort_unstable();
        Err(Some(failing_indices))
    }
}

/// Given a set of PCDs, their corresponding Marlin verification keys, and the DLogCommitterKey(s)
/// from both groups of our EC cycle, compute and return an accumulation proof(s) for
/// the dlog accumulators/"items".
/// In case of failure, returns the indices of the proofs which caused it (if possible).
/// The PCDs are allowed to use different size restrictions of the DLogCommitterKey
/// `g1_ck` and `g2_ck`.
pub fn accumulate_proofs<G1, G2, D: Digest>(
    pcds: &[GeneralPCD<G1, G2, D>],
    vks: &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    g1_ck: &DLogCommitterKey<G1>,
    g2_ck: &DLogCommitterKey<G2>,
) -> Result<(Option<AccumulationProof<G1>>, Option<AccumulationProof<G2>>), Option<Vec<usize>>>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>
        + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>
        + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let accumulation_time = start_timer!(|| "Accumulate proofs");

    // Get accumulators from pcds
    let (accs_g1, accs_g2) =
        get_accumulators::<G1, G2, D>(pcds, vks, g1_ck, g2_ck).map_err(|e| {
            end_timer!(accumulation_time);
            e
        })?;

    // Create accumulation proofs
    let acc_proof_g1 = if accs_g1.is_empty() {
        None
    } else {
        Some(
            DLogItemAccumulator::<G1, D>::accumulate_items(g1_ck, accs_g1)
                .map_err(|_| {
                    end_timer!(accumulation_time);
                    None
                })?
                .1,
        )
    };

    let acc_proof_g2 = if accs_g2.is_empty() {
        None
    } else {
        Some(
            DLogItemAccumulator::<G2, D>::accumulate_items(g2_ck, accs_g2)
                .map_err(|_| {
                    end_timer!(accumulation_time);
                    None
                })?
                .1,
        )
    };

    end_timer!(accumulation_time);

    Ok((acc_proof_g1, acc_proof_g2))
}

/// Verifies a set of PCDs which is augmented by an accumulation proof for their
/// dlog items. (This is cheaper than batch verification, as it doesn't need to
/// do any batching of witnesses.)
/// In case of failure, returns the indices of the proofs which caused it (if possible).
/// The PCDs are allowed to use different size restrictions of the DLogCommitterKey
/// `g1_ck` and `g2_ck`.
pub fn verify_aggregated_proofs<G1, G2, D: Digest, R: RngCore>(
    pcds: &[GeneralPCD<G1, G2, D>],
    vks: &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    accumulation_proof_g1: &Option<AccumulationProof<G1>>,
    accumulation_proof_g2: &Option<AccumulationProof<G2>>,
    g1_vk: &DLogVerifierKey<G1>,
    g2_vk: &DLogVerifierKey<G2>,
    rng: &mut R,
) -> Result<bool, Option<Vec<usize>>>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>
        + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>
        + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let verification_time = start_timer!(|| "Verify aggregated proofs");

    // Do the succinct verification of the PCDs and get their accumulators
    let (accs_g1, accs_g2) =
        get_accumulators::<G1, G2, D>(pcds, vks, g1_vk, g2_vk).map_err(|e| {
            end_timer!(verification_time);
            e
        })?;

    // fully verify the dlog aggregation proof in G1, if present.
    let result_accumulate_g1 = if accumulation_proof_g1.is_some() {
        let dummy_g1 = DLogItem::<G1>::default();
        DLogItemAccumulator::<G1, D>::verify_accumulated_items::<R>(
            &dummy_g1,
            g1_vk,
            accs_g1,
            accumulation_proof_g1.as_ref().unwrap(),
            rng,
        )
        .map_err(|_| {
            end_timer!(verification_time);
            None
        })?
    } else {
        true
    };

    // fully verify the dlog aggregation proof in G2, if present.
    let result_accumulate_g2 = if accumulation_proof_g2.is_some() {
        let dummy_g2 = DLogItem::<G2>::default();
        DLogItemAccumulator::<G2, D>::verify_accumulated_items::<R>(
            &dummy_g2,
            g2_vk,
            accs_g2,
            accumulation_proof_g2.as_ref().unwrap(),
            rng,
        )
        .map_err(|_| {
            end_timer!(verification_time);
            None
        })?
    } else {
        true
    };

    end_timer!(verification_time);

    Ok(result_accumulate_g1 && result_accumulate_g2)
}

/// Batch verification of PCDs consisting of FinalDarlin/SimpleMarlin PCDs.
/// The succinct parts are processed in serial, the dlog items (in both of the groups G1
/// and G2) are verified in batch.
/// In case of failure, returns the indices of the proofs which caused it (if possible).
/// The PCDs are allowed to use different size restrictions of the DLogCommitterKey
/// `g1_ck` and `g2_ck`.
pub fn batch_verify_proofs<G1, G2, D: Digest, R: RngCore>(
    pcds: &[GeneralPCD<G1, G2, D>],
    vks: &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    g1_vk: &DLogVerifierKey<G1>,
    g2_vk: &DLogVerifierKey<G2>,
    rng: &mut R,
) -> Result<bool, Option<Vec<usize>>>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>
        + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>
        + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let verification_time = start_timer!(|| "Batch verify proofs");

    // Do the succinct verification of the PCDs and get their accumulators
    let (accs_g1, accs_g2) =
        get_accumulators::<G1, G2, D>(pcds, vks, g1_vk, g2_vk).map_err(|e| {
            end_timer!(verification_time);
            e
        })?;

    // Verify accumulators (hard part)
    let result_g1 = if accs_g1.is_empty() {
        true
    } else {
        DLogItemAccumulator::<G1, D>::check_items::<R>(g1_vk, &accs_g1, rng).map_err(|_| {
            end_timer!(verification_time);
            None
        })?
    };

    let result_g2 = if accs_g2.is_empty() {
        true
    } else {
        DLogItemAccumulator::<G2, D>::check_items::<R>(g2_vk, &accs_g2, rng).map_err(|_| {
            end_timer!(verification_time);
            None
        })?
    };

    end_timer!(verification_time);

    Ok(result_g1 && result_g2)
}
