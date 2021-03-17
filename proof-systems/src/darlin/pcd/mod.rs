use algebra::AffineCurve;
use r1cs_core::ConstraintSynthesizer;
use marlin::{ProverKey as MarlinProverKey, VerifierKey as MarlinVerifierKey, Error as MarlinError, AHPForR1CS};
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC, UniversalParams,
        CommitterKey as DLogCommitterKey, VerifierKey as DLogVerifierKey,
    },
    PolynomialCommitment, Error as PCError
};
use crate::darlin::accumulators::Accumulator;
use rand::RngCore;
use digest::Digest;

pub mod simple_marlin;
pub mod final_darlin;

pub struct PCDParameters {
    pub segment_size: usize
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
        let index = AHPForR1CS::index(circuit)?;

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