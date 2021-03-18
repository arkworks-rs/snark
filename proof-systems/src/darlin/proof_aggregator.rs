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
        dlog::DLogAccumulator,
        Accumulator, AccumulationProof
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

pub(crate) fn get_accumulators<G1, G2, D: Digest>(
    pcds:      &[GeneralPCD<G1, G2, D>],
    vks:       &[MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>],
    g1_ck:     &DLogCommitterKey<G1>,
    g2_ck:     &DLogCommitterKey<G2>,
) -> Result<(Vec<DLogAccumulator<G1>>, Vec<DLogAccumulator<G2>>), Option<usize>>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let accs = pcds
        .into_par_iter()
        .zip(vks)
        .enumerate()
        .map(|(i, (pcd, vk))|
            {
                let vk = FinalDarlinPCDVerifierKey::<G1, G2, D>{
                    marlin_vk: vk,
                    dlog_vks: (g1_ck, g2_ck)
                };
                // No need to trim the vk here to the specific segment size used
                // to generate the proof for this pcd, as the IPA succinct_check
                // function doesn't use vk.comm_key at all.
                pcd.succinct_verify(&vk).map_err(|_| Some(i))
            }
        ).collect::<Result<Vec<_>, _>>()?;

    let accs_g1 = accs.iter().flat_map(|acc| acc.0.clone()).collect::<Vec<_>>();
    let accs_g2 = accs.into_iter().flat_map(|acc| acc.1).collect::<Vec<_>>();

    Ok((accs_g1, accs_g2))
}

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
    // Get accumulators from pcds
    let (accs_g1, accs_g2) = get_accumulators::<G1, G2, D>(pcds, vks, g1_ck, g2_ck)?;

    // Create accumulation proofs
    let acc_proof_g1 = if accs_g1.is_empty() {
        None
    } else {
        Some(
            DLogAccumulator::<G1>::accumulate::<D>(g1_ck, accs_g1)
            .map_err(|_| None)?
            .1
        )
    };

    let acc_proof_g2 = if accs_g2.is_empty() {
        None
    } else {
        Some(
            DLogAccumulator::<G2>::accumulate::<D>(g2_ck, accs_g2)
                .map_err(|_| None)?
                .1
        )
    };

    Ok((acc_proof_g1, acc_proof_g2))
}

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
    // Get accumulators from pcds
    let (accs_g1, accs_g2) = get_accumulators::<G1, G2, D>(pcds, vks, g1_vk, g2_vk)?;

    // Verify accumulators and accumulation proofs
    let result_accumulate_g1 = if accumulation_proof_g1.is_some() {
        let dummy_g1 = DLogAccumulator::<G1>::default();
        dummy_g1.verify_accumulate::<R, D>(
            g1_vk, accs_g1, accumulation_proof_g1.as_ref().unwrap(), rng
        ).map_err(|_| None)?
    } else {
        true
    };

    let result_accumulate_g2 = if accumulation_proof_g2.is_some() {
        let dummy_g2 = DLogAccumulator::<G2>::default();
        dummy_g2.verify_accumulate::<R, D>(
            g2_vk, accs_g2, accumulation_proof_g2.as_ref().unwrap(), rng
        ).map_err(|_| None)?
    } else {
        true
    };

    Ok(result_accumulate_g1 && result_accumulate_g2)
}

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
    // Get accumulators from pcds (perform succinct verification)
    let (accs_g1, accs_g2) = get_accumulators::<G1, G2, D>(pcds, vks, g1_vk, g2_vk)?;

    // Verify accumulators (hard part)
    let result_g1 = if accs_g1.is_empty() {
        true
    } else {
        DLogAccumulator::<G1>::check_accumulators::<R, D>(
            g1_vk, &accs_g1, rng
        ).map_err(|_| None)?
    };

    let result_g2 = if accs_g2.is_empty() {
        true
    } else {
        DLogAccumulator::<G2>::check_accumulators::<R, D>(
            g2_vk, &accs_g2, rng
        ).map_err(|_| None)?
    };

    Ok(result_g1 && result_g2)
}