use algebra::AffineCurve;
use digest::Digest;
use marlin::{
    VerifierKey as MarlinVerifierKey,
    Proof as MarlinProof, Marlin
};
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC, VerifierKey as DLogVerifierKey
    },
    rng::FiatShamirRng,
};
use crate::darlin::{
    pcd::{PCD, error::PCDError},
    accumulators::{
        dlog::{DLogItem, DLogItemAccumulator}, ItemAccumulator
    },
};
use poly_commit::ipa_pc::Commitment;
use std::marker::PhantomData;

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
pub struct SimpleMarlinPCD<'a, G: AffineCurve, D: Digest> {
    pub proof:                     MarlinProof<G::ScalarField, InnerProductArgPC<G, D>>,
    pub usr_ins:                   Vec<G::ScalarField>,
    _lifetime:                     PhantomData<&'a ()>,
}

impl<'a, G, D> SimpleMarlinPCD<'a, G, D>
    where
        G: AffineCurve,
        D: Digest + 'a,
{
    pub fn new(
        proof:   MarlinProof<G::ScalarField, InnerProductArgPC<G, D>>,
        usr_ins: Vec<G::ScalarField>
    ) -> Self
    {
        Self { proof, usr_ins, _lifetime: PhantomData }
    }
}

pub struct SimpleMarlinPCDVerifierKey<'a, G: AffineCurve, D: Digest>(
    pub &'a MarlinVerifierKey<G::ScalarField, InnerProductArgPC<G, D>>,
    pub &'a DLogVerifierKey<G>
);

impl<'a, G: AffineCurve, D: Digest> AsRef<DLogVerifierKey<G>> for SimpleMarlinPCDVerifierKey<'a, G, D> {
    fn as_ref(&self) -> &DLogVerifierKey<G> {
        &self.1
    }
}

impl<'a, G, D> PCD for SimpleMarlinPCD<'a, G, D>
    where
        G: AffineCurve,
        D: Digest + 'a,
{
    type PCDAccumulator = DLogItemAccumulator<G, D>;
    type PCDVerifierKey = SimpleMarlinPCDVerifierKey<'a, G, D>;

    fn succinct_verify(
        &self,
        vk: &Self::PCDVerifierKey,
    ) -> Result<<Self::PCDAccumulator as ItemAccumulator>::Item, PCDError>
    {
        let succinct_time = start_timer!(|| "Marlin succinct verifier");

        // Verify sumchecks
        let (query_set, evaluations, labeled_comms, mut fs_rng) = Marlin::<G::ScalarField, InnerProductArgPC<G, D>, D>::verify_ahp(
            &vk.0,
            self.usr_ins.as_slice(),
            &self.proof,
        ).map_err(|e| {
            end_timer!(succinct_time);
            PCDError::FailedSuccinctVerification(format!("{:?}", e))
        })?;

        // Absorb evaluations and sample new challenge
        fs_rng.absorb(&self.proof.evaluations);

        // Succinct verify DLOG proof
        let (xi_s, g_final) = InnerProductArgPC::<G, D>::succinct_batch_check_individual_opening_challenges(
            &vk.1,
            &labeled_comms,
            &query_set,
            &evaluations,
            &self.proof.pc_proof,
            &mut fs_rng,
        ).map_err(|e| {
            end_timer!(succinct_time);
            PCDError::FailedSuccinctVerification(e.to_string())
        })?;

        // Successfull verification: return current accumulator
        let acc = DLogItem::<G> {
            g_final: Commitment::<G> {  comm: vec![g_final], shifted_comm: None  },
            xi_s,
        };

        end_timer!(succinct_time);
        Ok(acc)
    }
}