use algebra::{AffineCurve, Field, ToBytes, to_bytes, UniformRand, ToConstraintField};
use digest::Digest;
use marlin::{
    MarlinConfig,
    VerifierKey as MarlinVerifierKey,
    Marlin, Proof as MarlinProof
};
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC, VerifierKey as DLogVerifierKey, Commitment
    },
    Error
};
use crate::darlin::accumulators::dlog::{DLogAccumulator, RecursiveDLogAccumulator};
use crate::darlin::PCD;
use rand::RngCore;
use std::marker::PhantomData;

// Maybe later we will deferr algebraic checks over G1::BaseField
pub struct FinalDarlinDeferredData<G1: AffineCurve, G2: AffineCurve> {
    pub(crate) previous_acc:       DLogAccumulator<G2>,
    pub(crate) pre_previous_acc:   DLogAccumulator<G1>,
}

/// FinalDarlinPCD with two deferred DLOG accumulators.
pub struct FinalDarlinPCD<G1: AffineCurve, G2: AffineCurve, D: Digest, MC: MarlinConfig + Send + Sync> {
    /// Full Marlin proof without deferred arithmetics in G1.
    marlin_proof:       MarlinProof<G1::ScalarField, InnerProductArgPC<G1, D>>,
    deferred:           FinalDarlinDeferredData<G1, G2>,
    usr_ins:            Vec<G1::ScalarField>,
    _config:            PhantomData<MC>,
}

//TODO: Use references
pub struct FinalDarlinPCDVerifierKey<G1: AffineCurve, G2: AffineCurve, D: Digest>(
    pub MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
    pub DLogVerifierKey<G1>,
    pub DLogVerifierKey<G2>,
);

impl<
    G1: AffineCurve,
    G2: AffineCurve,
    D: Digest
> AsRef<(DLogVerifierKey<G1>, DLogVerifierKey<G2>)> for FinalDarlinPCDVerifierKey<G1, G2, D> {
    fn as_ref(&self) -> &(DLogVerifierKey<G1>, DLogVerifierKey<G2>) {
        &(self.1.clone(), self.2.clone())
    }
}

impl<G1, G2, D, MC> PCD for FinalDarlinPCD<G1, G2, D, MC>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>,
    D: Digest,
    MC: MarlinConfig + Send + Sync,
{
    type PCDAccumulator = RecursiveDLogAccumulator<G1, G2>;
    type PCDVerifierKey = FinalDarlinPCDVerifierKey<G1, G2, D>;

    fn succinct_verify<R: RngCore>(
        &self,
        vk: &Self::PCDVerifierKey,
        rng: &mut R
    ) -> Result<Self::PCDAccumulator, Error>
    {
        let succinct_time = start_timer!(|| "Finalized Darlin succinct verifier");

        // Verify sumchecks

        // Get "system inputs"
        // TODO: Take into account the specific structure of the accumulators to save constraints
        // Convert accumulator to bytes
        let mut input_bytes = Vec::new();
        input_bytes.append(&mut to_bytes!(self.deferred.pre_previous_acc).unwrap());
        input_bytes.append(&mut to_bytes!(self.deferred.previous_acc).unwrap());

        // Extract C::F1 field element from input_bytes
        let mut public_inputs = input_bytes.as_slice().to_field_elements().unwrap();

        // Append user inputs
        public_inputs.append(&mut self.usr_ins.clone());

        let ahp_result = Marlin::<G1::ScalarField, InnerProductArgPC<G1, D>, D, MC>::verify_ahp(
            &vk.0,
            public_inputs.as_slice(),
            &self.marlin_proof,
        );

        if ahp_result.is_err() {
            println!("AHP verification failed");
            return Err(Error::FailedSuccinctCheck)
        }
        let (query_set, evaluations, labeled_comms, mut fs_rng) = ahp_result.unwrap();

        fs_rng.absorb(&self.marlin_proof.evaluations);
        let opening_challenge: G1::ScalarField = u128::rand(&mut fs_rng).into();
        let opening_challenges = |pow| opening_challenge.pow(&[pow]);

        // Succinct verify DLOG proof
        let succinct_result = InnerProductArgPC::<G1, D>::succinct_batch_check_individual_opening_challenges(
            &vk.1,
            &labeled_comms,
            &query_set,
            &evaluations,
            &self.marlin_proof.pc_proof.proof,
            &opening_challenges,
            rng
        );

        if succinct_result.is_err() {
            println!("Succinct verification failed: {:?}", succinct_result.err());
            return Err(Error::FailedSuccinctCheck)
        }

        let (xi_s, g_final) = succinct_result.unwrap();
        let acc = DLogAccumulator::<G1> {
            g_final: Commitment::<G1> { comm: g_final, shifted_comm: None},
            xi_s,
        };

        end_timer!(succinct_time);
        Ok(RecursiveDLogAccumulator::<G1, G2>(vec![acc, self.deferred.pre_previous_acc.clone()], vec![self.deferred.previous_acc.clone()]))
    }
}

