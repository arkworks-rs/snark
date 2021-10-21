//! Final Darlin proof carrying data. The final Darlin is the last node of the
//! exiting/conversion chain of our Darlin PCD scheme, and provides a (coboundary)
//! Marlin proof plus the dlog accumulators of the previous and pre-previous node.
use crate::darlin::{
    accumulators::dlog::{DLogItem, DualDLogItem, DualDLogItemAccumulator},
    accumulators::ItemAccumulator,
    data_structures::*,
    pcd::{error::PCDError, PCD},
    FinalDarlin, FinalDarlinVerifierKey,
};
use algebra::{AffineCurve, ToConstraintField};
use digest::Digest;
use poly_commit::{
    ipa_pc::{Commitment, InnerProductArgPC, VerifierKey as DLogVerifierKey},
    rng::FiatShamirRng,
};
use std::marker::PhantomData;

/// As every PCD, the `FinalDarlinPCD` comes as a proof plus "statement".
#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
pub struct FinalDarlinPCD<'a, G1: AffineCurve, G2: AffineCurve, D: Digest> {
    /// A `FinalDarlinProof` is a Marlin proof plus deferred dlog accumulators
    pub final_darlin_proof: FinalDarlinProof<G1, G2, D>,
    /// The user inputs form essentially the "statement" of the recursive proof.
    pub usr_ins: Vec<G1::ScalarField>,
    _lifetime: PhantomData<&'a ()>,
}

impl<'a, G1, G2, D> FinalDarlinPCD<'a, G1, G2, D>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>
        + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>
        + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
    D: Digest + 'a,
{
    pub fn new(
        final_darlin_proof: FinalDarlinProof<G1, G2, D>,
        usr_ins: Vec<G1::ScalarField>,
    ) -> Self {
        Self {
            final_darlin_proof,
            usr_ins,
            _lifetime: PhantomData,
        }
    }
}

/// To verify the PCD of a final Darlin we only need the `FinalDarlinVerifierKey` (or, the
/// IOP verifier key) of the final circuit and the two dlog committer keys for G1 and G2.
pub struct FinalDarlinPCDVerifierKey<'a, G1: AffineCurve, G2: AffineCurve, D: Digest> {
    pub final_darlin_vk: &'a FinalDarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
    pub dlog_vks: (&'a DLogVerifierKey<G1>, &'a DLogVerifierKey<G2>),
}

impl<'a, G1: AffineCurve, G2: AffineCurve, D: Digest>
    AsRef<(&'a DLogVerifierKey<G1>, &'a DLogVerifierKey<G2>)>
    for FinalDarlinPCDVerifierKey<'a, G1, G2, D>
{
    fn as_ref(&self) -> &(&'a DLogVerifierKey<G1>, &'a DLogVerifierKey<G2>) {
        &self.dlog_vks
    }
}

impl<'a, G1, G2, D> PCD for FinalDarlinPCD<'a, G1, G2, D>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>
        + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>
        + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
    D: Digest + 'a,
{
    type PCDAccumulator = DualDLogItemAccumulator<'a, G1, G2, D>;
    type PCDVerifierKey = FinalDarlinPCDVerifierKey<'a, G1, G2, D>;

    fn succinct_verify(
        &self,
        vk: &Self::PCDVerifierKey,
    ) -> Result<<Self::PCDAccumulator as ItemAccumulator>::Item, PCDError> {
        let succinct_time = start_timer!(|| "Finalized Darlin succinct verifier");

        // Verify sumchecks
        let (query_set, evaluations, labeled_comms, mut fs_rng) =
            FinalDarlin::<G1, G2, D>::verify_ahp(
                vk.dlog_vks.0,
                vk.final_darlin_vk,
                self.usr_ins.as_slice(),
                &self.final_darlin_proof,
            )
            .map_err(|e| {
                end_timer!(succinct_time);
                PCDError::FailedSuccinctVerification(format!("{:?}", e))
            })?;

        // Absorb evaluations and sample new challenge
        fs_rng.absorb(&self.final_darlin_proof.proof.evaluations);

        // Succinct verify DLOG proof
        let (xi_s, g_final) =
            InnerProductArgPC::<G1, D>::succinct_batch_check_individual_opening_challenges(
                vk.dlog_vks.0,
                &labeled_comms,
                &query_set,
                &evaluations,
                &self.final_darlin_proof.proof.pc_proof,
                &mut fs_rng,
            )
            .map_err(|e| {
                end_timer!(succinct_time);
                PCDError::FailedSuccinctVerification(e.to_string())
            })?;

        // Verification successfull: return new accumulator
        let acc = DLogItem::<G1> {
            g_final: Commitment::<G1> {
                comm: vec![g_final],
                shifted_comm: None,
            },
            xi_s,
        };

        end_timer!(succinct_time);
        Ok(DualDLogItem::<G1, G2>(
            vec![
                acc,
                self.final_darlin_proof.deferred.pre_previous_acc.clone(),
            ],
            vec![self.final_darlin_proof.deferred.previous_acc.clone()],
        ))
    }
}
