use algebra::{AffineCurve, Field, UniformRand};
use digest::Digest;
use marlin::{
    VerifierKey as MarlinVerifierKey,
    Proof as MarlinProof, MarlinConfig, Marlin
};
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC, VerifierKey as DLogVerifierKey
    },
    Error
};
use crate::darlin::PCD;
use rand::RngCore;
use crate::darlin::accumulators::dlog::DLogAccumulator;
use std::marker::PhantomData;
use poly_commit::ipa_pc::Commitment;

pub struct SimpleMarlinPCD<G: AffineCurve, D: Digest, MC: MarlinConfig> {
    proof:                     MarlinProof<G::ScalarField, InnerProductArgPC<G, D>>,
    usr_ins:                   Vec<G::ScalarField>,
    #[doc(hidden)] _config:    PhantomData<MC>
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

impl<'a, G, D, MC> PCD<'a> for SimpleMarlinPCD<G, D, MC>
    where
        G: AffineCurve,
        D: Digest + 'a,
        MC: MarlinConfig,
{
    type PCDAccumulator = DLogAccumulator<G>;
    type PCDVerifierKey = SimpleMarlinPCDVerifierKey<'a, G, D>;

    fn succinct_verify<R: RngCore>(
        &self,
        vk: &Self::PCDVerifierKey,
        rng: &mut R
    ) -> Result<Self::PCDAccumulator, Error>
    {
        let succinct_time = start_timer!(|| "Marlin succinct verifier");

        // Verify sumchecks
        let ahp_result = Marlin::<G::ScalarField, InnerProductArgPC<G, D>, D, MC>::verify_ahp(
            &vk.0,
            self.usr_ins.as_slice(),
            &self.proof,
        );

        if ahp_result.is_err() {
            println!("AHP verification failed");
            return Err(Error::FailedSuccinctCheck)
        }

        let (query_set, evaluations, labeled_comms, mut fs_rng) = ahp_result.unwrap();

        fs_rng.absorb(&self.proof.evaluations);
        let opening_challenge: G::ScalarField = u128::rand(&mut fs_rng).into();
        let opening_challenges = |pow| opening_challenge.pow(&[pow]);

        // Succinct verify DLOG proof
        let succinct_result = InnerProductArgPC::<G, D>::succinct_batch_check_individual_opening_challenges(
            &vk.1,
            &labeled_comms,
            &query_set,
            &evaluations,
            &self.proof.pc_proof.proof,
            &opening_challenges,
            rng
        );

        if succinct_result.is_err() {
            println!("Succinct verification failed: {:?}", succinct_result.err());
            return Err(Error::FailedSuccinctCheck)
        }

        let (xi_s, g_final) = succinct_result.unwrap();
        let acc = DLogAccumulator::<G> {
            g_final: Commitment::<G> {  comm: g_final, shifted_comm: None   },
            xi_s,
        };

        end_timer!(succinct_time);
        Ok(acc)
    }
}