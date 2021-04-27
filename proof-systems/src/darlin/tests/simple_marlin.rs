use algebra::{Field, AffineCurve, UniformRand};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use poly_commit::ipa_pc::{InnerProductArgPC, CommitterKey, UniversalParams};
use marlin::{
    Marlin, ProverKey as MarlinProverKey, VerifierKey as MarlinVerifierKey,
};
use crate::darlin::pcd::{
    PCDParameters, simple_marlin::{SimpleMarlinPCD, MarlinProof}
};
use rand::{ Rng, RngCore };
use digest::Digest;
use std::ops::MulAssign;

#[derive(Copy, Clone)]
pub(crate) struct Circuit<F: Field> {
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

        for i in 0..(self.num_variables - 5) {
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

#[allow(dead_code)]
pub fn generate_test_pcd<'a, G: AffineCurve, D: Digest + 'a, R: RngCore>(
    pc_ck: &CommitterKey<G>,
    marlin_pk: &MarlinProverKey<G::ScalarField, InnerProductArgPC<G, D>>,
    num_constraints: usize,
    zk: bool,
    rng: &mut R,
) -> SimpleMarlinPCD<'a, G, D>
{
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
        num_variables: num_constraints/2,
    };

    let proof = Marlin::<G::ScalarField, InnerProductArgPC<G, D>, D>::prove(
        marlin_pk,
        pc_ck,
        circ,
        zk,
        if zk { Some(rng) } else { None }
    ).unwrap();

    SimpleMarlinPCD::<'a, G, D>::new(MarlinProof::<G, D>(proof), vec![c, d])
}

#[allow(dead_code)]
pub fn generate_test_data<'a, G: AffineCurve, D: Digest + 'a, R: RngCore>(
    num_constraints: usize,
    segment_size: usize,
    params: &UniversalParams<G>,
    num_proofs: usize,
    rng: &mut R,
) -> (
    Vec<SimpleMarlinPCD<'a, G, D>>,
    Vec<MarlinVerifierKey<G::ScalarField, InnerProductArgPC<G, D>>>
)
{
    // Trim committer key and verifier key
    let config = PCDParameters { segment_size };
    let (committer_key, _) = config.universal_setup::<_, D>(params).unwrap();

    // Generate Marlin prover and verifier key
    let circ = Circuit {
        a: None,
        b: None,
        num_constraints,
        num_variables: num_constraints/2,
    };

    let (index_pk, index_vk) = Marlin::<G::ScalarField, InnerProductArgPC<G, D>, D>::index(
        &committer_key, circ.clone()
    ).unwrap();

    // Generate Marlin PCDs
    let simple_marlin_pcd = generate_test_pcd::<G, D, R>(
        &committer_key,
        &index_pk,
        num_constraints,
        rng.gen(),
        rng,
    );

    (vec![simple_marlin_pcd; num_proofs], vec![index_vk; num_proofs])
}