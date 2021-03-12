use algebra::{AffineCurve, Field, UniformRand, ToBits, ToConstraintField};
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
use crate::darlin::pcd::PCD;
use rand::RngCore;
use std::marker::PhantomData;

// Maybe later we will deferr algebraic checks over G1::BaseField
pub struct FinalDarlinDeferredData<G1: AffineCurve, G2: AffineCurve> {
    pub(crate) previous_acc:       DLogAccumulator<G2>,
    pub(crate) pre_previous_acc:   DLogAccumulator<G1>,
}

impl<G1, G2> ToConstraintField<G1::ScalarField> for FinalDarlinDeferredData<G1, G2>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    fn to_field_elements(&self) -> Result<Vec<G1::ScalarField>, Box<dyn std::error::Error>> {
        let mut fes = Vec::new();

        // Convert previous_acc into G1::ScalarField field elements
        let g_final_g2 = self.previous_acc.g_final.comm.clone();
        for c in g_final_g2.into_iter() {
            fes.append(&mut c.to_field_elements()?);
        }

        for fe in self.previous_acc.xi_s.0.clone().into_iter() {
            fes.append(&mut fe.write_bits().to_field_elements()?);
        }

        // Convert pre_previous_acc into G1::ScalarField field elements
        fes.append(&mut self.pre_previous_acc.xi_s.0.clone());

        let g_final_g1 = self.pre_previous_acc.g_final.comm.clone();
        for c in g_final_g1 {
            let c_fes = c.to_field_elements()?;
            for fe in c_fes {
                fes.append(&mut fe.write_bits().to_field_elements()?);
            }
        }

        Ok(fes)
    }
}

/// FinalDarlinPCD with two deferred DLOG accumulators.
pub struct FinalDarlinPCD<G1: AffineCurve, G2: AffineCurve, D: Digest, MC: MarlinConfig> {
    /// Full Marlin proof without deferred arithmetics in G1.
    marlin_proof:       MarlinProof<G1::ScalarField, InnerProductArgPC<G1, D>>,
    deferred:           FinalDarlinDeferredData<G1, G2>,
    usr_ins:            Vec<G1::ScalarField>,
    _config:            PhantomData<MC>,
}

pub struct FinalDarlinPCDVerifierKey<'a, G1: AffineCurve, G2: AffineCurve, D: Digest> {
    pub marlin_vk: &'a MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
    pub dlog_vks: (&'a DLogVerifierKey<G1>, &'a DLogVerifierKey<G2>)
}

impl<
    'a,
    G1: AffineCurve,
    G2: AffineCurve,
    D: Digest
> AsRef<(&'a DLogVerifierKey<G1>, &'a DLogVerifierKey<G2>)> for FinalDarlinPCDVerifierKey<'a, G1, G2, D> {
    fn as_ref(&self) -> &(&'a DLogVerifierKey<G1>, &'a DLogVerifierKey<G2>) {
        &self.dlog_vks
    }
}

impl<'a, G1, G2, D, MC> PCD<'a> for FinalDarlinPCD<G1, G2, D, MC>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
    D: Digest + 'a,
    MC: MarlinConfig,
{
    type PCDAccumulator = RecursiveDLogAccumulator<G1, G2>;
    type PCDVerifierKey = FinalDarlinPCDVerifierKey<'a, G1, G2, D>;

    fn succinct_verify<R: RngCore>(
        &self,
        vk: &Self::PCDVerifierKey,
        rng: &mut R
    ) -> Result<Self::PCDAccumulator, Error>
    {
        let succinct_time = start_timer!(|| "Finalized Darlin succinct verifier");

        // Verify sumchecks

        // Get "system inputs"
        let mut public_inputs = self.deferred.to_field_elements().unwrap();

        // Append user inputs
        public_inputs.append(&mut self.usr_ins.clone());

        let ahp_result = Marlin::<G1::ScalarField, InnerProductArgPC<G1, D>, D, MC>::verify_ahp(
            vk.marlin_vk,
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
            vk.dlog_vks.0,
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
            g_final: Commitment::<G1> { comm: vec![g_final], shifted_comm: None},
            xi_s,
        };

        end_timer!(succinct_time);
        Ok(RecursiveDLogAccumulator::<G1, G2>(vec![acc, self.deferred.pre_previous_acc.clone()], vec![self.deferred.previous_acc.clone()]))
    }
}

