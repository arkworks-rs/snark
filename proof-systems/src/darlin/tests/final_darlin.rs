use algebra::{AffineCurve, ToConstraintField, UniformRand};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use crate::darlin::pcd::{
        PCDParameters,
        final_darlin::{FinalDarlinDeferredData, FinalDarlinPCD}
};
use marlin::{
    Marlin, ProverKey as MarlinProverKey, VerifierKey as MarlinVerifierKey,
};
use poly_commit::ipa_pc::{InnerProductArgPC, CommitterKey, UniversalParams};
use rand::RngCore;
use digest::Digest;
use std::ops::MulAssign;
use r1cs_std::{
    alloc::AllocGadget,
    fields::fp::FpGadget,
    eq::EqGadget,
};

#[derive(Clone)]
pub struct Circuit<G1: AffineCurve, G2: AffineCurve> {
    pub a: Option<G1::ScalarField>,
    pub b: Option<G1::ScalarField>,
    pub num_constraints: usize,
    pub num_variables: usize,

    // Deferred elements (sys ins)
    pub deferred: FinalDarlinDeferredData<G1, G2>
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
        let deferred_as_native_fes = self.deferred.to_field_elements().unwrap();
        let deferred_len = deferred_as_native_fes.len();

        // Alloc deferred data as public input
        let mut deferred_input_gs = Vec::new();
        for (i, fe) in deferred_as_native_fes.iter().enumerate() {
            let ins_g = FpGadget::<G1::ScalarField>::alloc_input(
                cs.ns(|| format!("Alloc input deferred elem {}", i)),
                || Ok(fe)
            )?;
            deferred_input_gs.push(ins_g);
        }

        // Alloc deferred data as witness
        let mut deferred_gs = Vec::new();
        for (i, fe) in deferred_as_native_fes.into_iter().enumerate() {
            let witness_g = FpGadget::<G1::ScalarField>::alloc(
                cs.ns(|| format!("Alloc deferred elem {}", i)),
                || Ok(fe)
            )?;
            deferred_gs.push(witness_g);
        }

        // As a simple way to use the public inputs and being able to test cases where sys data
        // (i.e. the deferred accumulators) are wrong, enforce equality between witnesses and
        // public inputs
        let mut test_constraints = cs.num_constraints();
        for (i, (deferred_w, deferred_ins)) in deferred_input_gs.into_iter().zip(deferred_gs).enumerate() {
            deferred_w.enforce_equal(
                cs.ns(|| format!("enforce deferred equal {}", i)),
                &deferred_ins
            )?;
        }
        test_constraints = cs.num_constraints() - test_constraints;

        // The following is equal to the SimpleMarlin circuit

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

        for i in 0..(self.num_variables - 5 - (2 * deferred_len)) {
            let _ = cs.alloc(
                || format!("var {}", i),
                || self.a.ok_or(SynthesisError::AssignmentMissing),
            )?;
        }

        for i in 0..(self.num_constraints - 1 - test_constraints){
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
pub fn generate_test_pcd<G1: AffineCurve, G2:AffineCurve, D: Digest, R: RngCore>(
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
        num_variables: num_constraints/2,
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

#[allow(dead_code)]
pub fn generate_test_data<G1: AffineCurve, G2: AffineCurve, D: Digest, R: RngCore>(
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
    let deferred = FinalDarlinDeferredData::<G1, G2>::generate_random::<R, D>(
        rng,
        &committer_key_g1,
        &committer_key_g2,
    );

    // Generate Marlin prover and verifier key
    let circ = Circuit {
        a: None,
        b: None,
        num_constraints,
        num_variables: num_constraints/2,
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