use algebra::{AffineCurve, ProjectiveCurve, ToConstraintField, UniformRand};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use crate::darlin::{
    pcd::{
        PCDParameters,
        final_darlin::{FinalDarlinDeferredData, FinalDarlinPCD}
    },
    accumulators::dlog::DLogAccumulator
};
use marlin::{
    Marlin, ProverKey as MarlinProverKey, VerifierKey as MarlinVerifierKey,
};
use poly_commit::ipa_pc::{InnerProductArgPC, Commitment, CommitterKey, SuccinctCheckPolynomial, UniversalParams};
use rand::RngCore;
use digest::Digest;
use std::ops::MulAssign;

#[derive(Clone)]
struct Circuit<G1: AffineCurve, G2: AffineCurve> {
    a: Option<G1::ScalarField>,
    b: Option<G1::ScalarField>,
    num_constraints: usize,
    num_variables: usize,

    // Deferred elements (sys ins)
    deferred: FinalDarlinDeferredData<G1, G2>
}

impl<G1: AffineCurve, G2: AffineCurve> ConstraintSynthesizer<G1::ScalarField> for Circuit<G1, G2>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    fn generate_constraints<CS: ConstraintSystem<G1::ScalarField>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError>
    {
        // Alloc deferred data as public input and do nothing with them, it's
        // solely for testing purposes
        let deferred_as_native_fes = self.deferred.to_field_elements().unwrap();
        for (i, fe) in deferred_as_native_fes.into_iter().enumerate() {
            cs.alloc_input(
                || format!("Alloc input deferred elem {}", i),
                || Ok(fe)
            )?;
        }

        let a = cs.alloc(|| "a", || self.a.ok_or(SynthesisError::AssignmentMissing))?;
        let b = cs.alloc(|| "b", || self.b.ok_or(SynthesisError::AssignmentMissing))?;
        let c = cs.alloc_input(
            || "c",
            || {
                let mut a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
                let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;

                a.mul_assign(&b);
                Ok(a)
            },
        )?;
        let d = cs.alloc_input(
            || "d",
            || {
                let mut a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
                let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;

                a.mul_assign(&b);
                a.mul_assign(&b);
                Ok(a)
            },
        )?;

        for i in 0..(self.num_variables - 3) {
            let _ = cs.alloc(
                || format!("var {}", i),
                || self.a.ok_or(SynthesisError::AssignmentMissing),
            )?;
        }

        for i in 0..(self.num_constraints - 1){
            cs.enforce(
                || format!("constraint {}", i),
                |lc| lc + a,
                |lc| lc + b,
                |lc| lc + c,
            );
        }
        cs.enforce(
            || format!("constraint {}", self.num_constraints - 1),
            |lc| lc + c,
            |lc| lc + b,
            |lc| lc + d,
        );

        Ok(())
    }
}

pub(crate) fn generate_test_pcd<G1: AffineCurve, G2:AffineCurve, D: Digest, R: RngCore>(
    pc_ck_g1: &CommitterKey<G1>,
    deferred: FinalDarlinDeferredData<G1, G2>,
    marlin_pk: &MarlinProverKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
    num_constraints: usize,
    zk: bool,
    rng: &mut R,
) -> FinalDarlinPCD<G1, G2, D>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let a = G1::ScalarField::rand(rng);
    let b = G1::ScalarField::rand(rng);
    let mut c = a;
    c.mul_assign(&b);
    let mut d = c;
    d.mul_assign(&b);

    let circ = Circuit::<G1, G2> {
        a: Some(a),
        b: Some(b),
        num_constraints,
        num_variables: num_constraints,
        deferred: deferred.clone(),
    };

    let proof = Marlin::<G1::ScalarField, InnerProductArgPC<G1, D>, D>::prove(
        marlin_pk,
        pc_ck_g1,
        circ,
        zk,
        if zk { Some(rng) } else { None }
    ).unwrap();

    FinalDarlinPCD::<G1, G2, D> {
        marlin_proof: proof,
        deferred,
        usr_ins: vec![c, d]
    }
}

pub(crate) fn generate_test_data<G1: AffineCurve, G2: AffineCurve, D: Digest, R: RngCore>(
    num_constraints: usize,
    segment_size: usize,
    params_g1: &UniversalParams<G1>,
    params_g2: &UniversalParams<G2>,
    num_proofs: usize,
    rng: &mut R,
) -> (
    Vec<FinalDarlinPCD<G1, G2, D>>,
    Vec<MarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>>
)
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    // Trim committer key and verifier key
    let config = PCDParameters { segment_size };
    let (committer_key_g1, _) = config.universal_setup::<_, D>(params_g1).unwrap();
    let (committer_key_g2, _) = config.universal_setup::<_, D>(params_g2).unwrap();

    // Generate random (but valid) deferred data

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

    // Collect accumulators in deferred struct
    let deferred = FinalDarlinDeferredData::<G1, G2> {
        previous_acc: acc_g2,
        pre_previous_acc: acc_g1
    };

    // Generate Marlin prover and verifier key
    let circ = Circuit {
        a: None,
        b: None,
        num_constraints,
        num_variables: num_constraints,
        deferred: deferred.clone()
    };

    let (index_pk, index_vk) = config.circuit_specific_setup(circ.clone(), &committer_key_g1).unwrap();

    // Generate Final Darlin PCDs
    let final_darlin_pcd = generate_test_pcd::<G1, G2, D, R>(
        &committer_key_g1,
        deferred,
        &index_pk,
        num_constraints,
        false,
        rng,
    );

    (vec![final_darlin_pcd; num_proofs], vec![index_vk; num_proofs])
}