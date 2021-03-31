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

        // Succinct verify the commitments of the Bullet polys opened at the new challenge,
        // and get the new chals
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

    /// Batch verification of DLog accumulators: combine Bullet polys and the corresponding GFins
    /// and perform a single MSM.
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

    /// Our implementation is size optimized: the accumulation proof is all that it's needed since
    /// that the g_fin of the new accumulator can be derived from the AccumulationProof (which is
    /// a DLOG opening proof) and the xi_s can be recomputed via succinct verification of the previous
    /// accumulators.
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

        // Compute multi poly single point opening proof of the Bullet polys
        // at the GFin(s)
        let opening_proof = InnerProductArgPC::<G, D>::open_check_polys(
            &ck,
            xi_s.iter(),
            g_fins.iter(),
            z
        )?;

        end_timer!(poly_time);

        // Even if our implementation is size optimized, the API requires us to
        // return an accumulator too: so we return a dummy one instead (to be
        // discarded by the caller).
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

        // Succinct part: compute the aggregated accumulator by recomputing the xi_s from
        // the previous_accumulator and the g_fin from the accumulation proof
        let new_acc = Self::succinct_verify_accumulate::<D>(vk, previous_accumulators, proof)?;
        if new_acc.is_none() {
            end_timer!(check_acc_time);
            return Ok(false)
        }

        // Verify the aggregated accumulator
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

#[cfg(test)]
mod test {
    use super::*;
    use poly_commit::{QuerySet, Evaluations, LabeledPolynomial, ipa_pc::{
        BatchProof, UniversalParams,
    }, PolynomialCommitment};

    use rand::{distributions::Distribution, thread_rng};
    use std::marker::PhantomData;
    use digest::Digest;
    use blake2::Blake2s;


    #[derive(Copy, Clone, Default)]
    struct TestInfo {
        max_degree: Option<usize>,
        supported_degree: Option<usize>,
        num_polynomials: usize,
        enforce_degree_bounds: bool,
        max_num_queries: usize,
        segmented: bool
    }

    #[derive(Derivative)]
    #[derivative(Clone(bound = ""))]
    struct VerifierData<'a, G: AffineCurve> {
        vk:                      VerifierKey<G>,
        comms:                   Vec<LabeledCommitment<Commitment<G>>>,
        query_set:               QuerySet<'a, G::ScalarField>,
        values:                  Evaluations<'a, G::ScalarField>,
        proof:                   BatchProof<G>,
        opening_challenge:       G::ScalarField,
        polynomials:             Vec<LabeledPolynomial<G::ScalarField>>,
        num_polynomials:         usize,
        num_points_in_query_set: usize,
        _m:                      PhantomData<&'a G::ScalarField>, // To avoid compilation issue 'a
    }

    fn get_data_for_verifier<'a, G, D>(
        info: TestInfo,
        pp: Option<UniversalParams<G>>
    ) -> Result<VerifierData<'a, G>, Error>
        where
            G: AffineCurve,
            D: Digest
    {
        let TestInfo {
            max_degree,
            supported_degree,
            num_polynomials,
            enforce_degree_bounds,
            max_num_queries,
            segmented,
            ..
        } = info;

        let rng = &mut thread_rng();
        let max_degree =
            max_degree.unwrap_or(rand::distributions::Uniform::from(2..=64).sample(rng));
        let pp = if pp.is_some() { pp.unwrap() } else { InnerProductArgPC::<G, D>::setup(max_degree)? };

        let supported_degree = match supported_degree {
            Some(0) => 0,
            Some(d) => d,
            None => rand::distributions::Uniform::from(1..=max_degree).sample(rng)
        };
        assert!(
            max_degree >= supported_degree,
            "max_degree < supported_degree"
        );
        let mut polynomials = Vec::new();
        let mut degree_bounds = if enforce_degree_bounds {
            Some(Vec::new())
        } else {
            None
        };

        let seg_mul = rand::distributions::Uniform::from(5..=15).sample(rng);
        let mut labels = Vec::new();
        println!("Sampled supported degree");

        // Generate polynomials
        let num_points_in_query_set =
            rand::distributions::Uniform::from(1..=max_num_queries).sample(rng);
        for i in 0..num_polynomials {
            let label = format!("Test{}", i);
            labels.push(label.clone());
            let degree;
            if segmented {
                degree = if supported_degree > 0 {
                    rand::distributions::Uniform::from(1..=supported_degree).sample(rng)
                } else {
                    0
                } * seg_mul;
            } else {
                degree = if supported_degree > 0 {
                    rand::distributions::Uniform::from(1..=supported_degree).sample(rng)
                } else {
                    0
                }
            }
            let poly = Polynomial::rand(degree, rng);

            let degree_bound = if let Some(degree_bounds) = &mut degree_bounds {
                let degree_bound;
                if segmented {
                    degree_bound = degree;
                } else {
                    let range = rand::distributions::Uniform::from(degree..=supported_degree);
                    degree_bound = range.sample(rng);
                }
                degree_bounds.push(degree_bound);
                Some(degree_bound)
            } else {
                None
            };

            let hiding_bound = if num_points_in_query_set >= degree {
                Some(degree)
            } else {
                Some(num_points_in_query_set)
            };
            println!("Hiding bound: {:?}", hiding_bound);

            polynomials.push(LabeledPolynomial::new(
                label,
                poly,
                degree_bound,
                hiding_bound,
            ))
        }
        let supported_hiding_bound = polynomials
            .iter()
            .map(|p| p.hiding_bound().unwrap_or(0))
            .max()
            .unwrap_or(0);
        println!("supported degree: {:?}", supported_degree);
        println!("supported hiding bound: {:?}", supported_hiding_bound);
        println!("num_points_in_query_set: {:?}", num_points_in_query_set);
        let (ck, vk) = InnerProductArgPC::<G, D>::trim(
            &pp,
            supported_degree,
        )?;
        println!("Trimmed");

        let (comms, rands) = InnerProductArgPC::<G, D>::commit(&ck, &polynomials, Some(rng))?;

        // Construct query set
        let mut query_set = QuerySet::new();
        let mut values = Evaluations::new();
        // let mut point = F::one();
        for _ in 0..num_points_in_query_set {
            let point = G::ScalarField::rand(rng);
            for (i, label) in labels.iter().enumerate() {
                query_set.insert((label.clone(), (format!("{}", i), point)));
                let value = polynomials[i].evaluate(point);
                values.insert((label.clone(), point), value);
            }
        }
        println!("Generated query set");

        let opening_challenge = G::ScalarField::rand(rng);
        let proof = InnerProductArgPC::<G, D>::batch_open(
            &ck,
            &polynomials,
            &comms,
            &query_set,
            opening_challenge,
            &rands,
            Some(rng),
        )?;

        Ok(VerifierData {
            vk,
            comms,
            query_set,
            values,
            proof,
            opening_challenge,
            polynomials,
            num_polynomials,
            num_points_in_query_set,
            _m: PhantomData,
        })
    }

    fn accumulation_test<G, D>() -> Result<(), Error>
        where
            G: AffineCurve,
            D: Digest,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..=64).sample(rng);

        let info = TestInfo {
            max_degree: Some(max_degree),
            supported_degree: None,
            num_polynomials: 10,
            enforce_degree_bounds: true,
            max_num_queries: 5,
            segmented: true,
            ..Default::default()
        };

        let pp = InnerProductArgPC::<G, D>::setup(max_degree)?;

        for num_proofs in 1..10 {
            // Generate all proofs and the data needed by the verifier to verify them
            let verifier_data_vec = vec![get_data_for_verifier::<G, D>(info, Some(pp.clone())).unwrap(); num_proofs];

            let vk = &verifier_data_vec[0].vk;

            let mut comms = Vec::new();
            let mut query_sets = Vec::new();
            let mut evals = Vec::new();
            let mut proofs = Vec::new();
            let mut opening_challenges = Vec::new();

            verifier_data_vec.iter().for_each(|verifier_data| {
                assert_eq!(&verifier_data.vk, vk); // Vk should be equal for all proofs
                comms.push(verifier_data.comms.as_slice());
                query_sets.push(&verifier_data.query_set);
                evals.push(&verifier_data.values);
                proofs.push(&verifier_data.proof);
                opening_challenges.push(verifier_data.opening_challenge.clone());
            });

            // Prover side
            let (xi_s_vec, g_fins) = InnerProductArgPC::<G, D>::succinct_batch_check(
                vk,
                comms.clone(),
                query_sets.clone(),
                evals.clone(),
                proofs.clone(),
                opening_challenges.clone(),
            )?;

            let accumulators = xi_s_vec
                .into_iter()
                .zip(g_fins)
                .map(|(xi_s, g_final)| {
                    DLogAccumulator::<G> { g_final: Commitment::<G> {comm: vec![g_final], shifted_comm: None},  xi_s }
                }).collect::<Vec<_>>();

            let (_, proof) = DLogAccumulator::<G>::accumulate::<D>(
                vk,
                accumulators.clone(),
            )?;

            // Verifier side
            let dummy = DLogAccumulator::<G>::default();
            assert!(
                dummy.verify_accumulate::<_, D>(
                    vk,
                    // Actually the verifier should recompute the accumulators with the succinct verification
                    accumulators,
                    &proof,
                    rng
                )?
            );
        }
        Ok(())
    }

    fn batch_verification_test<G, D>() -> Result<(), Error>
        where
            G: AffineCurve,
            D: Digest,
    {
        let rng = &mut thread_rng();
        let max_degree = rand::distributions::Uniform::from(2..=64).sample(rng);

        let info = TestInfo {
            max_degree: Some(max_degree),
            supported_degree: None,
            num_polynomials: 10,
            enforce_degree_bounds: true,
            max_num_queries: 5,
            segmented: true,
            ..Default::default()
        };

        let pp = InnerProductArgPC::<G, D>::setup(max_degree)?;

        for num_proofs in 1..10 {
            // Generate all proofs and the data needed by the verifier to verify them
            let verifier_data_vec = vec![get_data_for_verifier::<G, D>(info, Some(pp.clone())).unwrap(); num_proofs];

            let vk = &verifier_data_vec[0].vk;

            let mut comms = Vec::new();
            let mut query_sets = Vec::new();
            let mut evals = Vec::new();
            let mut proofs = Vec::new();
            let mut opening_challenges = Vec::new();

            verifier_data_vec.iter().for_each(|verifier_data| {
                assert_eq!(&verifier_data.vk, vk); // Vk should be equal for all proofs
                comms.push(verifier_data.comms.as_slice());
                query_sets.push(&verifier_data.query_set);
                evals.push(&verifier_data.values);
                proofs.push(&verifier_data.proof);
                opening_challenges.push(verifier_data.opening_challenge.clone());
            });

            let (xi_s_vec, g_fins) = InnerProductArgPC::<G, D>::succinct_batch_check(
                vk,
                comms.clone(),
                query_sets.clone(),
                evals.clone(),
                proofs.clone(),
                opening_challenges.clone(),
            )?;

            let accumulators = xi_s_vec
                .into_iter()
                .zip(g_fins)
                .map(|(xi_s, g_final)| {
                    DLogAccumulator::<G> { g_final: Commitment::<G> {comm: vec![g_final], shifted_comm: None},  xi_s }
                }).collect::<Vec<_>>();

            assert!(
                DLogAccumulator::<G>::check_accumulators::<_, D>(
                    vk,
                    &accumulators,
                    rng
                )?
            );
        }
        Ok(())
    }

    use algebra::curves::tweedle::{
        dee::Affine as TweedleDee,
        dum::Affine as TweedleDum,
    };

    #[test]
    fn test_tweedle_accumulate_verify() {
        accumulation_test::<TweedleDee, Blake2s>().unwrap();
        accumulation_test::<TweedleDum, Blake2s>().unwrap();
    }

    #[test]
    fn test_tweedle_batch_verify() {
        batch_verification_test::<TweedleDee, Blake2s>().unwrap();
        batch_verification_test::<TweedleDum, Blake2s>().unwrap();
    }
}