use algebra::{AffineCurve, ProjectiveCurve, Field, UniformRand, ToBits, ToConstraintField};
use digest::Digest;
use marlin::{
    VerifierKey as MarlinVerifierKey,
    Marlin, Proof as MarlinProof
};
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC,
        CommitterKey as DLogCommitterKey,
        VerifierKey as DLogVerifierKey,
        Commitment, SuccinctCheckPolynomial
    },
    Error
};
use crate::darlin::accumulators::dlog::{DLogAccumulator, DualDLogAccumulator};
use crate::darlin::pcd::PCD;
use rand::RngCore;

// Maybe later this will include more element as we will deferr algebraic checks over G1::BaseField
#[derive(Clone)]
pub struct FinalDarlinDeferredData<G1: AffineCurve, G2: AffineCurve> {
    pub(crate) previous_acc:       DLogAccumulator<G2>,
    pub(crate) pre_previous_acc:   DLogAccumulator<G1>,
}

impl<G1, G2> FinalDarlinDeferredData<G1, G2>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    pub fn generate_random<R: RngCore, D: Digest>(
        rng: &mut R,
        committer_key_g1: &DLogCommitterKey<G1>,
        committer_key_g2: &DLogCommitterKey<G2>
    ) -> Self
    {
        // Generate valid accumulator over G1 starting from random xi_s
        let log_key_len_g1 = algebra::log2(committer_key_g1.comm_key.len());
        let random_xi_s_g1 = SuccinctCheckPolynomial::<G1::ScalarField>(vec![G1::ScalarField::rand(rng); log_key_len_g1 as usize]);
        let g_final_g1 = InnerProductArgPC::<G1, D>::cm_commit(
            committer_key_g1.comm_key.as_slice(),
            random_xi_s_g1.compute_coeffs().as_slice(),
            None,
            None,
        );

        let acc_g1 = DLogAccumulator::<G1> {
            g_final: Commitment::<G1> {comm: vec![g_final_g1.into_affine()], shifted_comm: None },
            xi_s: random_xi_s_g1
        };

        // Generate valid accumulator over G2 starting from random xi_s
        let log_key_len_g2 = algebra::log2(committer_key_g2.comm_key.len());
        let random_xi_s_g2 = SuccinctCheckPolynomial::<G2::ScalarField>(vec![G2::ScalarField::rand(rng); log_key_len_g2 as usize]);

        let g_final_g2 = InnerProductArgPC::<G2, D>::cm_commit(
            committer_key_g2.comm_key.as_slice(),
            random_xi_s_g2.compute_coeffs().as_slice(),
            None,
            None,
        );

        let acc_g2 = DLogAccumulator::<G2> {
            g_final: Commitment::<G2> {comm: vec![g_final_g2.into_affine()], shifted_comm: None },
            xi_s: random_xi_s_g2
        };

        // Return accumulators in deferred struct
        Self {
            previous_acc: acc_g2,
            pre_previous_acc: acc_g1
        }
    }
}

impl<G1, G2> ToConstraintField<G1::ScalarField> for FinalDarlinDeferredData<G1, G2>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    fn to_field_elements(&self) -> Result<Vec<G1::ScalarField>, Box<dyn std::error::Error>> {
        let mut fes = Vec::new();

        // Convert previous_acc into G1::ScalarField field elements (matching field)

        let g_final_g2 = self.previous_acc.g_final.comm.clone();
        for c in g_final_g2.into_iter() {
            fes.append(&mut c.to_field_elements()?);
        }

        // Convert xi_s to G1::ScalarField elements, the field is not matching:
        // serialize them all to bits and pack them safely into native field elements
        // TODO: Later the xi_s will be only 128 bits long, so no need to take all
        //       the bits from each one.
        let mut xi_s_bits = Vec::new();
        for fe in self.previous_acc.xi_s.0.clone().into_iter() {
            xi_s_bits.append(&mut fe.write_bits());
        }
        fes.append(&mut xi_s_bits.to_field_elements()?);

        // Convert g_final_g1 to G1::ScalarField elements, the field is not matching:
        // serialize them all to bits and pack them safely into native field elements

        let g_final_g1 = self.pre_previous_acc.g_final.comm.clone();
        let mut g_final_g1_bits = Vec::new();
        for c in g_final_g1 {
            let c_fes = c.to_field_elements()?;
            for fe in c_fes {
                g_final_g1_bits.append(&mut fe.write_bits());
            }
        }
        fes.append(&mut g_final_g1_bits.to_field_elements()?);

        // Convert pre_previous_acc into G1::ScalarField field elements (matching field)
        fes.append(&mut self.pre_previous_acc.xi_s.0.clone());

        Ok(fes)
    }
}

/// FinalDarlinPCD with two deferred DLOG accumulators.
#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
pub struct FinalDarlinPCD<G1: AffineCurve, G2: AffineCurve, D: Digest> {
    /// Full Marlin proof without deferred arithmetics in G1.
    pub marlin_proof:       MarlinProof<G1::ScalarField, InnerProductArgPC<G1, D>>,
    pub deferred:           FinalDarlinDeferredData<G1, G2>,
    pub usr_ins:            Vec<G1::ScalarField>,
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

impl<'a, G1, G2, D> PCD<'a> for FinalDarlinPCD<G1, G2, D>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
    D: Digest + 'a,
{
    type PCDAccumulator = DualDLogAccumulator<G1, G2>;
    type PCDVerifierKey = FinalDarlinPCDVerifierKey<'a, G1, G2, D>;

    fn succinct_verify(
        &self,
        vk: &Self::PCDVerifierKey,
    ) -> Result<Self::PCDAccumulator, Error>
    {
        let succinct_time = start_timer!(|| "Finalized Darlin succinct verifier");

        // Verify sumchecks

        // Get "system inputs"
        let mut public_inputs = self.deferred.to_field_elements().unwrap();

        // Append user inputs
        public_inputs.append(&mut self.usr_ins.clone());

        let ahp_result = Marlin::<G1::ScalarField, InnerProductArgPC<G1, D>, D>::verify_ahp(
            vk.marlin_vk,
            public_inputs.as_slice(),
            &self.marlin_proof,
        );

        if ahp_result.is_err() {
            println!("AHP verification failed");
            end_timer!(succinct_time);
            return Err(Error::FailedSuccinctCheck)
        }
        let (query_set, evaluations, labeled_comms, mut fs_rng) = ahp_result.unwrap();

        // Absorb evaluations and sample new challenge
        fs_rng.absorb(&self.marlin_proof.evaluations);
        let opening_challenge: G1::ScalarField = u128::rand(&mut fs_rng).into();
        let opening_challenges = |pow| opening_challenge.pow(&[pow]);

        // Succinct verify DLOG proof
        let succinct_result = InnerProductArgPC::<G1, D>::succinct_batch_check_individual_opening_challenges(
            vk.dlog_vks.0,
            &labeled_comms,
            &query_set,
            &evaluations,
            &self.marlin_proof.pc_proof,
            &opening_challenges,
        );

        if succinct_result.is_err() {
            println!("Succinct verification failed: {:?}", succinct_result.err());
            end_timer!(succinct_time);
            return Err(Error::FailedSuccinctCheck)
        }

        // Verification successfull: return new accumulator
        let (xi_s, g_final) = succinct_result.unwrap();
        let acc = DLogAccumulator::<G1> {
            g_final: Commitment::<G1> { comm: vec![g_final], shifted_comm: None},
            xi_s,
        };

        end_timer!(succinct_time);
        Ok(DualDLogAccumulator::<G1, G2>(vec![acc, self.deferred.pre_previous_acc.clone()], vec![self.deferred.previous_acc.clone()]))
    }
}

