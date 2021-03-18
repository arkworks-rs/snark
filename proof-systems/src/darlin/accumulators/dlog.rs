use algebra::{Field, AffineCurve, ProjectiveCurve, ToBytes, to_bytes, UniformRand};
use algebra_utils::polynomial::DensePolynomial as Polynomial;
use poly_commit::{ipa_pc::{
    InnerProductArgPC,
    Commitment,
    VerifierKey, CommitterKey,
    SuccinctCheckPolynomial,
}, LabeledCommitment, Error};
use crate::darlin::accumulators::{
    Accumulator, AccumulationProof,
};
use rayon::prelude::*;
use rand::RngCore;
use digest::Digest;

/// This implements the public aggregator for the IPA/DLOG commitment scheme.
#[derive(Clone)]
pub struct DLogAccumulator<G: AffineCurve> {
    /// Final committer key after the DLOG reduction.
    pub(crate) g_final:     Commitment<G>,

    /// Challenges of the DLOG reduction.
    pub(crate) xi_s:        SuccinctCheckPolynomial<G::ScalarField>,
}

impl<G: AffineCurve> DLogAccumulator<G> {

    /// This implementation reflects the procedure for size-optimized DLOG accumulation proofs
    /// in which we need to recompute the xi_s. If successfull, returns the new accumulator
    /// (the GFinal is taken from proof itself).
    pub fn succinct_verify_accumulate<D: Digest>(
        vk:                    &VerifierKey<G>,
        previous_accumulators: Vec<Self>,
        proof:                 &AccumulationProof<G>,
    ) -> Result<Option<Self>, Error>
    {
        let succinct_time = start_timer!(|| "Succinct verify accumulate");

        let poly_time = start_timer!(|| "Compute Bullet Polys evaluations");

        // Sample a new challenge z
        let z = InnerProductArgPC::<G, D>::compute_random_oracle_challenge(
            &to_bytes![vk.hash.clone(), previous_accumulators.as_slice()].unwrap(),
        );

        let comms_values = previous_accumulators
            .into_par_iter()
            .enumerate()
            .map(|(i, acc)| {
                let final_comm_key = acc.g_final.comm.clone();
                let xi_s = acc.xi_s;

                // Create a LabeledCommitment out of the g_final
                let labeled_comm = {
                    let comm = Commitment {
                        comm: final_comm_key,
                        shifted_comm: None
                    };

                    LabeledCommitment::new(
                        format!("check_poly_{}", i),
                        comm,
                        None,
                    )
                };

                // Compute the evaluation of the Bullet polynomial at z starting from the xi_s
                let eval = xi_s.evaluate(z);

                (labeled_comm, eval)
            }).collect::<Vec<_>>();

        // Save the evaluations into a separate vec
        let values = comms_values.iter().map(|(_, val)| val.clone()).collect::<Vec<_>>();

        // Save comms into a separate vector
        let comms = comms_values.into_iter().map(|(comm, _)| comm).collect::<Vec<_>>();

        end_timer!(poly_time);

        let check_time = start_timer!(|| "Succinct check IPA proof");

        // Sample new opening challenge
        let opening_challenge = InnerProductArgPC::<G, D>::compute_random_oracle_challenge(
            &values.iter().flat_map(|val| to_bytes!(val).unwrap()).collect::<Vec<_>>()
        );
        let opening_challenges = |pow| opening_challenge.pow(&[pow]);

        // Succinct check
        let xi_s = InnerProductArgPC::<G, D>::succinct_check(
            vk, comms.iter(), z, values, &proof.pc_proof, &opening_challenges
        )?;

        end_timer!(check_time);
        end_timer!(succinct_time);

        if xi_s.is_some() {
            Ok(Some(Self {
                g_final: Commitment::<G>{ comm: vec![proof.pc_proof.final_comm_key.clone()], shifted_comm: None },
                xi_s: xi_s.unwrap(),
            }))
        } else {
            Ok(None)
        }
    }
}

impl<G: AffineCurve> Default for DLogAccumulator<G> {
    fn default() -> Self {
        Self {
            g_final: Commitment::<G>::default(),
            xi_s: SuccinctCheckPolynomial(vec![])
        }
    }
}

impl<G: AffineCurve> ToBytes for DLogAccumulator<G> {
    fn write<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        self.g_final.write(&mut writer)?;
        self.xi_s.0.write(&mut writer)
    }
}

impl<'a, G: AffineCurve> Accumulator<'a> for DLogAccumulator<G> {
    type AccumulatorProverKey = CommitterKey<G>;
    type AccumulatorVerifierKey = VerifierKey<G>;
    type AccumulationProof = AccumulationProof<G>;

    fn check_accumulators<R: RngCore, D: Digest>(
        vk: &Self::AccumulatorVerifierKey,
        accumulators: &[Self],
        rng: &mut R
    ) -> Result<bool, Error>
    {
        let check_time = start_timer!(|| "Check accumulators");

        let final_comm_keys = accumulators.iter().flat_map(|acc| acc.g_final.comm.clone()).collect::<Vec<_>>();
        let xi_s_vec = accumulators.iter().map(|acc| acc.xi_s.clone()).collect::<Vec<_>>();

        let batching_time = start_timer!(|| "Combine check polynomials and final comm keys");

        // Sample batching challenge
        let random_scalar = G::ScalarField::rand(rng);
        let mut batching_chal = G::ScalarField::one();

        // Collect the powers of the batching challenge in a vector
        let mut batching_chal_pows = vec![G::ScalarField::zero(); xi_s_vec.len()];
        for i in 0..batching_chal_pows.len() {
            batching_chal_pows[i] = batching_chal;
            batching_chal *= &random_scalar;
        }

        // Compute the combined_check_poly
        let combined_check_poly = batching_chal_pows
            .par_iter()
            .zip(xi_s_vec)
            .map(|(&chal, xi_s)| {
                Polynomial::from_coefficients_vec(xi_s.compute_scaled_coeffs(-chal))
            }).reduce(|| Polynomial::zero(), |acc, scaled_poly| &acc + &scaled_poly);
        end_timer!(batching_time);

        // DLOG hard part.
        // The equation to check would be:
        // lambda_1 * gfin_1 + ... + lambda_n * gfin_n - combined_h_1 * g_vk_1 - ... - combined_h_m * g_vk_m = 0
        // Where combined_h_i = lambda_1 * h_1_i + ... + lambda_n * h_n_i
        // We do final verification and the batching of the GFin in a single MSM
        let hard_time = start_timer!(|| "Batch verify hard parts");
        let final_val = InnerProductArgPC::<G, D>::cm_commit(
            // The vk might be oversized, but the VariableBaseMSM function, will "trim"
            // the bases in order to be as big as the scalars vector, so no need to explicitly
            // trim the vk here.
            &[final_comm_keys.as_slice(), vk.comm_key.as_slice()].concat(),
            &[batching_chal_pows.as_slice(), combined_check_poly.coeffs.as_slice()].concat(),
            None,
            None,
        );
        end_timer!(hard_time);

        if !ProjectiveCurve::is_zero(&final_val) {
            end_timer!(check_time);
            return Ok(false);
        }
        end_timer!(check_time);
        Ok(true)
    }

    fn accumulate<D: Digest>(
        ck: &Self::AccumulatorProverKey,
        accumulators: Vec<Self>,
    ) -> Result<(Self, Self::AccumulationProof), Error>
    {
        let accumulate_time = start_timer!(|| "Accumulate");

        // Sample a new challenge point z
        let z = InnerProductArgPC::<G, D>::compute_random_oracle_challenge(
            &to_bytes![ck.hash.clone(), accumulators.as_slice()].unwrap(),
        );

        // Collect GFinals from the accumulators
        let g_fins = accumulators.iter().map(|acc| {
            Commitment::<G> {
                comm: acc.g_final.comm.clone(),
                shifted_comm: None
            }
        }).collect::<Vec<_>>();

        // Collect xi_s from the accumulators
        let xi_s = accumulators.into_iter().map(|acc| {
            acc.xi_s
        }).collect::<Vec<_>>();

        let poly_time = start_timer!(|| "Open Bullet Polys");

        let opening_proof = InnerProductArgPC::<G, D>::open_check_polys(
            &ck,
            xi_s.iter(),
            g_fins.iter(),
            z
        )?;

        end_timer!(poly_time);

        let accumulator = DLogAccumulator::<G>::default();

        let mut accumulation_proof = AccumulationProof::<G>::default();
        accumulation_proof.pc_proof = opening_proof;

        end_timer!(accumulate_time);

        Ok((accumulator, accumulation_proof))
    }

    fn verify_accumulate<R: RngCore, D: Digest>(
        &self,
        vk: &Self::AccumulatorVerifierKey,
        previous_accumulators: Vec<Self>,
        proof: &Self::AccumulationProof,
        rng: &mut R
    ) -> Result<bool, Error>
    {
        let check_acc_time = start_timer!(|| "Verify Accumulation");

        // Succinct part
        let new_acc = Self::succinct_verify_accumulate::<D>(vk, previous_accumulators, proof)?;
        if new_acc.is_none() {
            end_timer!(check_acc_time);
            return Ok(false)
        }

        // Hard part
        let hard_time = start_timer!(|| "DLOG hard part");
        let result = Self::check_accumulators::<R, D>(vk, &vec![new_acc.unwrap()], rng)?;
        end_timer!(hard_time);

        end_timer!(check_acc_time);

        Ok(result)
    }
}

pub struct DualDLogAccumulator<G1: AffineCurve, G2: AffineCurve>(
    pub(crate) Vec<DLogAccumulator<G1>>,
    pub(crate) Vec<DLogAccumulator<G2>>,
);

impl<G1: AffineCurve, G2: AffineCurve> ToBytes for DualDLogAccumulator<G1, G2> {
    fn write<W: std::io::Write>(&self, mut writer: W) -> std::io::Result<()> {
        self.0.write(&mut writer)?;
        self.1.write(&mut writer)
    }
}

impl<'a, G1, G2> Accumulator<'a> for DualDLogAccumulator<G1, G2>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField>,
{
    type AccumulatorProverKey = (&'a CommitterKey<G1>, &'a CommitterKey<G2>);
    type AccumulatorVerifierKey = (&'a VerifierKey<G1>, &'a VerifierKey<G2>);
    type AccumulationProof = (AccumulationProof<G1>, AccumulationProof<G2>);

    fn check_accumulators<R: RngCore, D: Digest>(
        vk: &Self::AccumulatorVerifierKey,
        accumulators: &[Self],
        rng: &mut R
    ) -> Result<bool, Error>
    {
        let g1_accumulators = accumulators.iter().flat_map(|acc| { acc.0.clone() }).collect::<Vec<_>>();
        if !DLogAccumulator::<G1>::check_accumulators::<R, D>(&vk.0, g1_accumulators.as_slice(), rng)? {
            return Ok(false)
        }

        let g2_accumulators = accumulators.iter().flat_map(|acc| { acc.1.clone() }).collect::<Vec<_>>();
        if !DLogAccumulator::<G2>::check_accumulators::<R, D>(&vk.1, g2_accumulators.as_slice(), rng)? {
            return Ok(false)
        }

        Ok(true)
    }

    fn accumulate<D: Digest>(
        ck: &Self::AccumulatorProverKey,
        accumulators: Vec<Self>,
    ) -> Result<(Self, Self::AccumulationProof), Error>
    {
        let g1_accumulators = accumulators.iter().flat_map(|acc| { acc.0.clone() }).collect::<Vec<_>>();
        let (_, g1_acc_proof) = DLogAccumulator::<G1>::accumulate::<D>(&ck.0, g1_accumulators)?;

        let g2_accumulators = accumulators.into_iter().flat_map(|acc| { acc.1 }).collect::<Vec<_>>();
        let (_, g2_acc_proof) = DLogAccumulator::<G2>::accumulate::<D>(&ck.1, g2_accumulators)?;

        let accumulator = DualDLogAccumulator::<G1, G2>(vec![DLogAccumulator::<G1>::default()], vec![DLogAccumulator::<G2>::default()]);
        let accumulation_proof = (g1_acc_proof, g2_acc_proof);

        Ok((accumulator, accumulation_proof))
    }

    fn verify_accumulate<R: RngCore, D: Digest>(
        &self,
        vk: &Self::AccumulatorVerifierKey,
        previous_accumulators: Vec<Self>,
        proof: &Self::AccumulationProof,
        rng: &mut R
    ) -> Result<bool, Error>
    {
        let g1_accumulators = previous_accumulators.iter().flat_map(|acc| { acc.0.clone() }).collect();
        if !(&DLogAccumulator::<G1>::default()).verify_accumulate::<R, D>(&vk.0, g1_accumulators, &proof.0, rng)? {
            return Ok(false)
        }

        let g2_accumulators = previous_accumulators.into_iter().flat_map(|acc| { acc.1 }).collect();
        if !(&DLogAccumulator::<G2>::default()).verify_accumulate::<R, D>(&vk.1, g2_accumulators, &proof.1, rng)? {
            return Ok(false)
        }

        Ok(true)
    }
}