pub mod pcd;
pub mod accumulators;
pub mod proof_aggregator;
pub mod data_structures;
pub mod error;

pub mod tests;

use algebra::{AffineCurve, ToConstraintField};
use poly_commit::{ ipa_pc::{
    UniversalParams, InnerProductArgPC,
    CommitterKey as DLogProverKey,
    VerifierKey as DLogVerifierKey,
    Commitment
}, PolynomialCommitment, QuerySet, LabeledCommitment, Evaluations};
use marlin::{
    Marlin,
    ProverKey as MarlinProverKey, VerifierKey as MarlinVerifierKey,
};
use crate::darlin::{
    data_structures::*,
    pcd::{
        PCD, PCDCircuit,
        simple_marlin::MarlinProof,
        final_darlin::{FinalDarlinPCD, FinalDarlinPCDVerifierKey}
    },
    error::FinalDarlinError,
};
use rand::RngCore;
use digest::Digest;
use std::marker::PhantomData;


/// FinalDarlin proving system. For now, mainly a wrapper for SimpleMarlin calls.
pub type FinalDarlinProverKey<F, PC> = MarlinProverKey<F, PC>;
pub type FinalDarlinVerifierKey<F, PC> = MarlinVerifierKey<F, PC>;

pub struct FinalDarlin<'a, G1: AffineCurve, G2: AffineCurve, D: Digest>(
    #[doc(hidden)] PhantomData<G1>,
    #[doc(hidden)] PhantomData<G2>,
    #[doc(hidden)] PhantomData<D>,
    #[doc(hidden)] PhantomData<&'a ()>
);

impl<'a, G1, G2, D>FinalDarlin<'a, G1, G2, D>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
        D:  Digest + 'a,
{
    /// Generate the universal prover and verifier keys for the argument system.
    pub fn universal_setup(
        num_constraints: usize,
        num_variables: usize,
        num_non_zero: usize,
        zk:  bool,
    ) -> Result<(
            UniversalParams<G1>,
            UniversalParams<G2>,
        ), FinalDarlinError>
    {
        let srs_g1 = Marlin::<G1::ScalarField, InnerProductArgPC<G1, D>, D>::universal_setup(
            num_constraints,
            num_variables,
            num_non_zero,
            zk
        )?;

        let srs_g2 = Marlin::<G2::ScalarField, InnerProductArgPC<G2, D>, D>::universal_setup(
            num_constraints,
            num_variables,
            num_non_zero,
            zk
        )?;

        Ok((srs_g1, srs_g2))
    }

    /// Generate the index-specific (i.e., circuit-specific) prover and verifier
    /// keys. This is a deterministic algorithm that anyone can rerun.
    pub fn index<C: PCDCircuit<G1>>(
        committer_key: &DLogProverKey<G1>,
        config:   C::SetupData,
    ) -> Result<(
            FinalDarlinProverKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
            FinalDarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
        ), FinalDarlinError>
    {
        let c = C::init(config);
        let res = Marlin::<G1::ScalarField, InnerProductArgPC<G1, D>, D>::index(committer_key, c)?;

        Ok(res)
    }

    /// Create and return a FinalDarlinPCD, given previous PCDs and a PCDCircuit that verify them
    /// along with some incremental data.
    pub fn prove<C>(
        index_pk:         &FinalDarlinProverKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
        pc_pk:            &DLogProverKey<G1>,
        config:           C::SetupData,
        // In future, this will be explicitly a RainbowDarlinPCD
        previous:         Vec<C::PreviousPCD>,
        previous_vks:     Vec<<C::PreviousPCD as PCD>::PCDVerifierKey>,
        incremental_data: C::IncrementalData,
        zk:               bool,
        zk_rng:           Option<&mut dyn RngCore>,
    ) -> Result<FinalDarlinPCD<'a, G1, G2, D>, FinalDarlinError>
        where
            C: PCDCircuit<G1, SystemInputs = FinalDarlinDeferredData<G1, G2>>,
    {
        let c = C::init_state(
            config,
            previous,
            previous_vks,
            incremental_data
        );

        let sys_ins = c.get_sys_ins()?.clone();

        let usr_ins = c.get_usr_ins()?;

        let proof = Marlin::<G1::ScalarField, InnerProductArgPC<G1, D>, D>::prove(
            index_pk, pc_pk, c, zk, zk_rng
        )?;

        let proof = FinalDarlinProof::<G1, G2, D> { proof: MarlinProof(proof), deferred: sys_ins };
        let usr_ins = usr_ins
            .to_field_elements()
            .map_err(|_| FinalDarlinError::Other("Failed to convert usr ins to field elements".to_owned()))?;

        Ok(FinalDarlinPCD::<G1, G2, D>::new(proof, usr_ins))
    }

    /// Verify that a proof for the constrain system defined by `C` asserts that
    /// all constraints are satisfied.
    pub fn verify<R: RngCore>(
        index_vk:     &FinalDarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
        pc_vk_g1:     &DLogVerifierKey<G1>,
        pc_vk_g2:     &DLogVerifierKey<G2>,
        usr_ins:      &[G1::ScalarField],
        proof:        &FinalDarlinProof<G1, G2, D>,
        rng:          &mut R,
    ) -> Result<bool, FinalDarlinError>
    {
        let final_darlin_pcd = FinalDarlinPCD::<G1, G2, D>::new(
            proof.clone(), usr_ins.to_vec()
        );

        let final_darlin_pcd_vk = FinalDarlinPCDVerifierKey::<G1, G2, D>{
            final_darlin_vk: index_vk,
            dlog_vks: (pc_vk_g1, pc_vk_g2)
        };

        let res = final_darlin_pcd.verify(&final_darlin_pcd_vk, rng)?;

        Ok(res)
    }

    /// Verify that a proof for the constrain system defined by `C` asserts that
    /// all constraints are satisfied. Checks only that the sumcheck equations
    /// are satisfied.
    pub fn verify_ahp(
        index_vk:       &FinalDarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
        usr_ins:        &[G1::ScalarField],
        proof:          &FinalDarlinProof<G1, G2, D>,
    )  -> Result<(
        QuerySet<'a, G1::ScalarField>,
        Evaluations<'a, G1::ScalarField>,
        Vec<LabeledCommitment<Commitment<G1>>>,
        <InnerProductArgPC<G1, D> as PolynomialCommitment<G1::ScalarField>>::RandomOracle,
    ), FinalDarlinError>
    {
        // Get "system inputs"
        let mut public_inputs = proof.deferred.to_field_elements().unwrap();

        // Append user inputs
        public_inputs.extend_from_slice(usr_ins);

        // Verify AHP
        let res = Marlin::<G1::ScalarField, InnerProductArgPC<G1, D>, D>::verify_ahp(
            index_vk, public_inputs.as_slice(), &proof.proof
        )?;

        Ok(res)
    }

    /// Verify that a proof for the constrain system defined by `C` asserts that
    /// all constraints are satisfied. Checks only that the opening proof is
    /// satisfied.
    pub fn verify_opening(
        pc_vk:          &DLogVerifierKey<G1>,
        proof:          &FinalDarlinProof<G1, G2, D>,
        labeled_comms:  Vec<LabeledCommitment<Commitment<G1>>>,
        query_set:      QuerySet<'a, G1::ScalarField>,
        evaluations:    Evaluations<'a, G1::ScalarField>,
        fs_rng:         &mut <InnerProductArgPC<G1, D> as PolynomialCommitment<G1::ScalarField>>::RandomOracle,
    ) -> Result<bool, FinalDarlinError>
    {
        let res = Marlin::<G1::ScalarField, InnerProductArgPC<G1, D>, D>::verify_opening(
            pc_vk, &proof.proof, labeled_comms, query_set, evaluations, fs_rng
        )?;

        Ok(res)
    }
}