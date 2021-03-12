use algebra::AffineCurve;
use r1cs_core::ConstraintSynthesizer;
use marlin::{MarlinConfig, ProverKey as MarlinProverKey, VerifierKey as MarlinVerifierKey, Error as MarlinError, AHPForR1CS};
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC, UniversalParams,
        CommitterKey as DLogCommitterKey, VerifierKey as DLogVerifierKey,
    },
    PolynomialCommitment, Error as PCError
};
use crate::darlin::{
    accumulators::{
        dlog::DLogAccumulator,
        Accumulator, AccumulationProof
    },
    simple_marlin::{SimpleMarlinPCD, SimpleMarlinPCDVerifierKey},
    finalized_darlin::{FinalDarlinPCD, FinalDarlinPCDVerifierKey}
};
use rand::{RngCore, thread_rng};
use digest::Digest;
use rayon::prelude::*;

//TODO: Remove dependency from R: RngCore when not needed

pub mod accumulators;
pub mod simple_marlin;
pub mod finalized_darlin;

pub struct PCDParameters {
    segment_size: usize
}

impl PCDParameters {
    pub fn universal_setup<G: AffineCurve, D: Digest>(
        &self,
        params: &UniversalParams<G>
    ) -> Result<(DLogCommitterKey<G>, DLogVerifierKey<G>), PCError>
    {
        InnerProductArgPC::<G, D>::trim(
            params,
            self.segment_size,
            0,
            None
        )
    }

    pub fn circuit_specific_setup<G: AffineCurve, C: ConstraintSynthesizer<G::ScalarField>, D: Digest>(
        &self,
        circuit: C,
        ck: &DLogCommitterKey<G>,
    ) -> Result<
        (
            MarlinProverKey<G::ScalarField, InnerProductArgPC<G, D>>,
            MarlinVerifierKey<G::ScalarField, InnerProductArgPC<G, D>>
        ), MarlinError<PCError>
        >
    {
        let index = AHPForR1CS::<G::ScalarField>::index::<C>(circuit)?;

        let (index_comms, index_comm_rands): (_, _) =
            InnerProductArgPC::<G, D>::commit(ck, index.iter(), None).map_err(MarlinError::from_pc_err)?;

        let index_comms = index_comms
            .into_iter()
            .map(|c| c.commitment().clone())
            .collect();
        let index_vk = MarlinVerifierKey {
            index_info: index.index_info,
            index_comms,
        };

        let index_pk = MarlinProverKey {
            index,
            index_comm_rands,
            index_vk: index_vk.clone(),
        };

        Ok((index_pk, index_vk))
    }
}

/// this trait expresses the functions for proof carrying data, in which the PCD is assumed
/// to be a set of data consisting of a statement, some deferred elements and a proof.
pub trait PCD<'a>: Sized + Send + Sync {
    type PCDAccumulator: Accumulator<'a>;
    type PCDVerifierKey: AsRef<<Self::PCDAccumulator as Accumulator<'a>>::AccumulatorVerifierKey>;

    //TODO: IN ORDER TO ALLOW SIDECHAINS TO CHOOSE ARBITRARILY THE SEGMENT SIZE:
    //      In the IPA succinct verification code we should remove any reference
    //      to the vk length/supported degree and derive the info we need directly
    //      from the length of the L and R vector (length = segment size)
    fn succinct_verify<R: RngCore>(
        &self,
        vk:         &Self::PCDVerifierKey,
        rng:        &mut R,
    ) -> Result<Self::PCDAccumulator, PCError>;

    fn hard_verify<R: RngCore, D: Digest>(
        &self,
        acc:    Self::PCDAccumulator,
        vk:     &Self::PCDVerifierKey,
        rng:    &mut R,
    ) -> Result<bool, PCError>
    { <Self::PCDAccumulator as Accumulator>::check_accumulators::<R, D>(vk.as_ref(), &[acc], rng) }

    fn verify<R: RngCore, D: Digest>(
        &self,
        vk:         &Self::PCDVerifierKey,
        rng:        &mut R,
    ) -> Result<bool, PCError>
    {
        let acc = self.succinct_verify::<R>(vk, rng)?;
        self.hard_verify::<R, D>(acc, vk, rng)
    }
}

//TODO: Make SimpleMarlinPCD and FinalDarlinPCD "polymorhpic" ?

//TODO: Get rid of this MarlinConfig template as it obliges all the proofs to have the same
//      MarlinConfig. Either we keep it inside the proof or the PCD (by converting it to a struct)
//      or we remove it (regarding LC_OPT, we never use it; regarding ZK, only the prover needs to
//      know about it; regarding SegmentSize it can be done in a transparent way)

//TODO: Do the same with Digest template

//TODO: Consider using RecursiveDLogAccumulator (maybe renaming it into DualDLogAccumulator)
//      to recycle its code in here
pub fn accumulate_proofs<'a, G1, G2, D: Digest, R: RngCore, MC: MarlinConfig + Sync + Send>(
    final_darlin_pcds:      impl IntoIterator<Item = &'a FinalDarlinPCD<G1, G2, D, MC>>,
    final_darlin_vks:       impl IntoIterator<Item = &'a MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>>,
    marlin_pcds:            impl IntoIterator<Item = &'a SimpleMarlinPCD<G1, D, MC>>,
    marlin_vks:             impl IntoIterator<Item = &'a MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>>,
    g1_ck:                  &DLogCommitterKey<G1>,
    g2_ck:                  &DLogCommitterKey<G2>,
    rng:                    &mut R
) -> Result<
    (
        AccumulationProof<G1>,
        AccumulationProof<G2>,
    ), PCError>
where
    SimpleMarlinPCD<G1, D, MC>: 'a,
    FinalDarlinPCD<G1, G2, D, MC>: 'a,
    MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>: 'a,
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>,
{
    // Final Darlin
    let final_darlin_pcds = final_darlin_pcds.into_iter().collect::<Vec<_>>();
    let final_darlin_vks = final_darlin_vks.into_iter().collect::<Vec<_>>();

    let darlin_accs = final_darlin_pcds
        .into_par_iter()
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
    let accs_g2 = darlin_accs.into_iter().flat_map(|acc| acc.1.clone()).collect::<Vec<_>>();

    // Marlin
    let marlin_pcds = marlin_pcds.into_iter().collect::<Vec<_>>();
    let marlin_vks = marlin_vks.into_iter().collect::<Vec<_>>();

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

pub fn verify_aggregated_proofs<'a, G1, G2, D: Digest, R: RngCore, MC: MarlinConfig>(
    final_darlin_pcds:      impl IntoIterator<Item = &'a FinalDarlinPCD<G1, G2, D, MC>>,
    final_darlin_vks:       impl IntoIterator<Item = &'a MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>>,
    marlin_pcds:            impl IntoIterator<Item = &'a SimpleMarlinPCD<G1, D, MC>>,
    marlin_vks:             impl IntoIterator<Item = &'a MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>>,
    accumulation_proof_g1:  &AccumulationProof<G1>,
    accumulation_proof_g2:  &AccumulationProof<G2>,
    g1_vk:                  &DLogVerifierKey<G1>,
    g2_vk:                  &DLogVerifierKey<G2>,
    rng:                    &mut R
) -> Result<bool, PCError>
    where
        SimpleMarlinPCD<G1, D, MC>: 'a,
        FinalDarlinPCD<G1, G2, D, MC>: 'a,
        MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>: 'a,
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>,
{
    // Final Darlin
    let final_darlin_pcds = final_darlin_pcds.into_iter().collect::<Vec<_>>();
    let final_darlin_vks = final_darlin_vks.into_iter().collect::<Vec<_>>();

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
    let accs_g2 = darlin_accs.into_iter().flat_map(|acc| acc.1.clone()).collect::<Vec<_>>();

    // Marlin
    let marlin_pcds = marlin_pcds.into_iter().collect::<Vec<_>>();
    let marlin_vks = marlin_vks.into_iter().collect::<Vec<_>>();

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