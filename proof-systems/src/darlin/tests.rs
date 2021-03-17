use algebra::{Field, fields::tweedle::Fr, curves::tweedle::{
    dee::Affine as DeeAffine, dum::Affine as DumAffine
}, AffineCurve, UniformRand};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use marlin::{
    Marlin, ProverKey as MarlinProverKey,
};
use poly_commit::ipa_pc::{InnerProductArgPC, CommitterKey};
use blake2::Blake2s;
use crate::darlin::pcd::{PCDParameters, simple_marlin::SimpleMarlinPCD};
use digest::Digest;
use rand::thread_rng;
use poly_commit::PolynomialCommitment;
use crate::darlin::{accumulate_proofs, verify_aggregated_proofs};
use std::ops::MulAssign;

type TestIPAPCDee = InnerProductArgPC<DeeAffine, Blake2s>;
type TestIPAPCDum = InnerProductArgPC<DumAffine, Blake2s>;

#[derive(Copy, Clone)]
struct Circuit<F: Field> {
    a: Option<F>,
    b: Option<F>,
    num_constraints: usize,
    num_variables: usize,
}

impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for Circuit<ConstraintF> {
    fn generate_constraints<CS: ConstraintSystem<ConstraintF>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
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

fn test_circuit<G: AffineCurve, D: Digest>(
    pc_ck: &CommitterKey<G>,
    marlin_pk: &MarlinProverKey<G::ScalarField, InnerProductArgPC<G, D>>,
    num_constraints: usize,
    zk: bool,
) -> SimpleMarlinPCD<G, D>
{
    let rng = &mut thread_rng();

    let a = G::ScalarField::rand(rng);
    let b = G::ScalarField::rand(rng);
    let mut c = a;
    c.mul_assign(&b);
    let mut d = c;
    d.mul_assign(&b);

    let circ = Circuit {
        a: Some(a),
        b: Some(b),
        num_constraints,
        num_variables: num_constraints,
    };

    let proof = Marlin::<G::ScalarField, InnerProductArgPC<G, D>, D>::prove(
        marlin_pk,
        pc_ck,
        circ,
        zk,
        &mut if zk { Some(thread_rng()) } else { None }
    ).unwrap();

    SimpleMarlinPCD::<G, D> {
        proof,
        usr_ins: vec![c, d]
    }
}


#[test]
fn test_accumulate_verify_simple_marlin_fixed_segment_size() {

    let rng = &mut thread_rng();

    // Set params
    let num_constraints = 100;
    let segment_size = num_constraints;

    // Generate committer key and verifier key
    let config = PCDParameters { segment_size };
    let (committer_key_g1, verifier_key_g1) = config.universal_setup::<_, Blake2s>(
        &TestIPAPCDee::setup(num_constraints, rng).unwrap()
    ).unwrap();
    let (committer_key_g2, verifier_key_g2) = config.universal_setup::<_, Blake2s>(
        &TestIPAPCDum::setup(num_constraints, rng).unwrap()
    ).unwrap();

    // Generate Marlin prover and verifier key
    let circ = Circuit {
        a: None,
        b: None,
        num_constraints,
        num_variables: num_constraints,
    };

    let (index_pk, index_vk) = config.circuit_specific_setup(circ.clone(), &committer_key_g1).unwrap();

    // Generate Marlin PCDs
    let samples = 100usize;
    let simple_marlin_pcd = test_circuit::<DeeAffine, Blake2s>(
        &committer_key_g1,
        &index_pk,
        num_constraints,
        false
    );
    let mut simple_marlin_pcds = vec![simple_marlin_pcd; samples];
    let simple_marlin_vks = vec![index_vk; samples];
    
    // Accumulate PCDs
    let (proof_g1, _) = accumulate_proofs::<DeeAffine, DumAffine, Blake2s>(
        &[],
        &[],
        simple_marlin_pcds.as_slice(),
        simple_marlin_vks.as_slice(),
        &committer_key_g1,
        &committer_key_g2
    ).unwrap();

    let proof_g1 = proof_g1.unwrap();

    // Verify accumulation
    assert!(verify_aggregated_proofs::<DeeAffine, DumAffine, Blake2s, _>(
        &[],
        &[],
        simple_marlin_pcds.as_slice(),
        simple_marlin_vks.as_slice(),
        Some(&proof_g1),
        None,
        &verifier_key_g1,
        &verifier_key_g2,
        rng
    ).unwrap());

    // Pass wrong accumulation proof and check verification fails
    // Change one element in proof_g1
    let mut wrong_proof_g1 = proof_g1.clone();
    wrong_proof_g1.pc_proof.c = Fr::rand(rng);
    assert!(!verify_aggregated_proofs::<DeeAffine, DumAffine, Blake2s, _>(
        &[],
        &[],
        simple_marlin_pcds.as_slice(),
        simple_marlin_vks.as_slice(),
        Some(&wrong_proof_g1),
        None,
        &verifier_key_g1,
        &verifier_key_g2,
        &mut thread_rng()
    ).unwrap());

    // Pass wrong public inputs for one marlin PCD and check verification of accumulation fail
    let wrong_usr_ins = vec![Fr::rand(rng); simple_marlin_pcds[0].usr_ins.len()];
    simple_marlin_pcds[0].usr_ins = wrong_usr_ins;

    assert!(verify_aggregated_proofs::<DeeAffine, DumAffine, Blake2s, _>(
        &[],
        &[],
        simple_marlin_pcds.as_slice(),
        simple_marlin_vks.as_slice(),
        Some(&proof_g1),
        None,
        &verifier_key_g1,
        &verifier_key_g2,
        &mut thread_rng()
    ).is_err());
}