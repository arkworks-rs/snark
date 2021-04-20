use algebra::{
    AffineCurve, ToConstraintField
};
use marlin::VerifierKey as MarlinVerifierKey;
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC,
        CommitterKey as DLogCommitterKey, VerifierKey as DLogVerifierKey,
    },
};
use crate::darlin::{
    accumulators::{
        dlog::DLogItem,
        ItemAccumulator, AccumulationProof
    },
    pcd::{
        PCD,
        final_darlin::FinalDarlinPCDVerifierKey,
    },
};
use rand::RngCore;
use digest::Digest;
use rayon::prelude::*;
use crate::darlin::pcd::GeneralPCD;
use crate::darlin::accumulators::dlog::DLogItemAccumulator;

/// Given a set of PCDs, their corresponding Marlin verification keys, and the DLogCommitterKey(s)
/// over two groups of a curve cycle, compute and return the associated accumulators via the
/// succinct verification of them; in case of failure, if it's possible to estabilish it, return the
/// index of the proof that has caused the failure. It is implicitly allowed for the PCDs to be
/// produced (thus verified) using DLogCommitterKey of different sizes, as long as they are smaller
/// equal than `g1_ck` and `g2_ck`.
pub(crate) fn get_accumulators<G1, G2, D: Digest>(
    pcds:      &[GeneralPCD<G1, G2, D>],
    vks:       &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    g1_ck:     &DLogCommitterKey<G1>,
    g2_ck:     &DLogCommitterKey<G2>,
) -> Result<(Vec<DLogItem<G1>>, Vec<DLogItem<G2>>), Option<usize>>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let accumulators_time = start_timer!(|| "Compute accumulators");
    let accs = pcds
        .into_par_iter()
        .zip(vks)
        .enumerate()
        .map(|(i, (pcd, vk))|
            {
                let vk = FinalDarlinPCDVerifierKey::<G1, G2, D>{
                    final_darlin_vk: vk,
                    dlog_vks: (g1_ck, g2_ck)
                };
                // No need to trim the vk here to the specific segment size used
                // to generate the proof for this pcd, as the IPA succinct_check
                // function doesn't use vk.comm_key at all.
                pcd.succinct_verify(&vk).map_err(|_| Some(i))
            }
        ).collect::<Result<Vec<_>, _>>().map_err(|e| {
            end_timer!(accumulators_time);
            e
        })?;

    let accs_g1 = accs.iter().flat_map(|acc| acc.0.clone()).collect::<Vec<_>>();
    let accs_g2 = accs.into_iter().flat_map(|acc| acc.1).collect::<Vec<_>>();

    end_timer!(accumulators_time);

    Ok((accs_g1, accs_g2))
}

/// Given a set of PCDs, their corresponding Marlin verification keys, and the DLogCommitterKey(s)
/// over two groups of a curve cycle, compute and return the accumulation proofs for the PCDs
/// on both groups (if needed); in case of failure, if it's possible to estabilish it, return the
/// index of the proof that has caused the failure. It is implicitly allowed for the PCDs to be
/// produced using DLogCommitterKey of different sizes, as long as they are smaller equal than
/// `g1_ck` and `g2_ck`.
pub fn accumulate_proofs<G1, G2, D: Digest>(
    pcds:      &[GeneralPCD<G1, G2, D>],
    vks:       &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    g1_ck:     &DLogCommitterKey<G1>,
    g2_ck:     &DLogCommitterKey<G2>,
) -> Result<
    (
        Option<AccumulationProof<G1>>,
        Option<AccumulationProof<G2>>,
    ), Option<usize>>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let accumulation_time = start_timer!(|| "Accumulate proofs");

    // Get accumulators from pcds
    let (accs_g1, accs_g2) = get_accumulators::<G1, G2, D>(pcds, vks, g1_ck, g2_ck)
        .map_err(|e| {
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
            .1
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
                .1
        )
    };

    end_timer!(accumulation_time);

    Ok((acc_proof_g1, acc_proof_g2))
}

/// Given a set of PCDs, their corresponding Marlin verification keys, the DLogVerifierKey(s)
/// over two groups of a curve cycle, and an accumulation proof in both (if needed) of the groups
/// verify the latters; in case of failure, if it's possible to estabilish it, return the
/// index of the proof that has caused the failure. It is implicitly allowed for the PCDs to be
/// produced (thus verified) using DLogVerifierKey of different sizes, as long as they are
/// smaller equal than `g1_vk` and `g2_vk`.
pub fn verify_aggregated_proofs<G1, G2, D: Digest, R: RngCore>(
    pcds:                   &[GeneralPCD<G1, G2, D>],
    vks:                    &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    accumulation_proof_g1:  &Option<AccumulationProof<G1>>,
    accumulation_proof_g2:  &Option<AccumulationProof<G2>>,
    g1_vk:                  &DLogVerifierKey<G1>,
    g2_vk:                  &DLogVerifierKey<G2>,
    rng:                    &mut R
) -> Result<bool, Option<usize>>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let verification_time = start_timer!(|| "Verify aggregated proofs");

    // Get accumulators from pcds
    let (accs_g1, accs_g2) = get_accumulators::<G1, G2, D>(pcds, vks, g1_vk, g2_vk)
        .map_err(|e| {
            end_timer!(verification_time);
            e
        })?;

    // Verify accumulators and accumulation proofs
    let result_accumulate_g1 = if accumulation_proof_g1.is_some() {
        let dummy_g1 = DLogItem::<G1>::default();
        DLogItemAccumulator::<G1, D>::verify_accumulated_items::<R>(
            &dummy_g1, g1_vk, accs_g1, accumulation_proof_g1.as_ref().unwrap(), rng
        ).map_err(|_| {
            end_timer!(verification_time);
            None
        })?
    } else {
        true
    };

    let result_accumulate_g2 = if accumulation_proof_g2.is_some() {
        let dummy_g2 = DLogItem::<G2>::default();
        DLogItemAccumulator::<G2, D>::verify_accumulated_items::<R>(
            &dummy_g2, g2_vk, accs_g2, accumulation_proof_g2.as_ref().unwrap(), rng
        ).map_err(|_| {
            end_timer!(verification_time);
            None
        })?
    } else {
        true
    };

    end_timer!(verification_time);

    Ok(result_accumulate_g1 && result_accumulate_g2)
}

/// Given a set of PCDs, their corresponding Marlin verification keys, the DLogVerifierKey(s)
/// over two groups of a curve cycle, verify the proofs; in case of failure, if it's possible
/// to estabilish it, return the index of the proof that has caused the failure. It is
/// implicitly allowed for the PCDs to be produced (thus verified) using DLogVerifierKey of
/// different sizes, as long as they are smaller equal than `g1_vk` and `g2_vk`.
pub fn batch_verify_proofs<G1, G2, D: Digest, R: RngCore>(
    pcds:                   &[GeneralPCD<G1, G2, D>],
    vks:                    &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    g1_vk:                  &DLogVerifierKey<G1>,
    g2_vk:                  &DLogVerifierKey<G2>,
    rng:                    &mut R
) -> Result<bool, Option<usize>>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let verification_time = start_timer!(|| "Batch verify proofs");

    // Get accumulators from pcds (perform succinct verification)
    let (accs_g1, accs_g2) = get_accumulators::<G1, G2, D>(pcds, vks, g1_vk, g2_vk)
        .map_err(|e| {
            end_timer!(verification_time);
            e
        })?;

    // Verify accumulators (hard part)
    let result_g1 = if accs_g1.is_empty() {
        true
    } else {
        DLogItemAccumulator::<G1, D>::check_items::<R>(
            g1_vk, &accs_g1, rng
        ).map_err(|_| {
            end_timer!(verification_time);
            None
        })?
    };

    let result_g2 = if accs_g2.is_empty() {
        true
    } else {
        DLogItemAccumulator::<G2, D>::check_items::<R>(
            g2_vk, &accs_g2, rng
        ).map_err(|_| {
            end_timer!(verification_time);
            None
        })?
    };

    end_timer!(verification_time);

    Ok(result_g1 && result_g2)
}