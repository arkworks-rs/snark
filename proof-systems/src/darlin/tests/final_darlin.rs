use algebra::{AffineCurve, ToConstraintField, UniformRand};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use crate::darlin::{
    pcd::{
        PCD, PCDParameters, PCDCircuit,
        final_darlin::FinalDarlinPCD,
        error::PCDError,
    },
    accumulators::ItemAccumulator,
    data_structures::FinalDarlinDeferredData,
    FinalDarlinProverKey, FinalDarlinVerifierKey, FinalDarlin,
};
use poly_commit::{
    ipa_pc::{InnerProductArgPC, CommitterKey, UniversalParams},
    Error as PCError
};
use rand::{ Rng, RngCore };
use digest::Digest;
use r1cs_std::{
    alloc::AllocGadget,
    fields::fp::FpGadget,
    eq::EqGadget,
};
use std::ops::MulAssign;

// Dummy Acc used for testing
pub struct TestAcc {}

impl ItemAccumulator for TestAcc {
    type AccumulatorProverKey = ();
    type AccumulatorVerifierKey = ();
    type AccumulationProof = ();
    type Item = ();

    fn check_items<R: RngCore>(
        _vk: &Self::AccumulatorVerifierKey,
        _accumulators: &[Self::Item],
        _rng: &mut R
    ) -> Result<bool, PCError>
    {
        Ok(true)
    }

    fn accumulate_items(
        _ck: &Self::AccumulatorProverKey,
        _accumulators: Vec<Self::Item>
    ) -> Result<(Self::Item, Self::AccumulationProof), PCError>
    {
        Ok(((), ()))
    }

    fn verify_accumulated_items<R: RngCore>(
        _current_accumulator: &Self::Item,
        _vk: &Self::AccumulatorVerifierKey,
        _previous_accumulators: Vec<Self::Item>,
        _proof: &Self::AccumulationProof,
        _rng: &mut R
    ) -> Result<bool, PCError>
    {
        Ok(true)
    }
}

// Test PCDVk
pub struct TestPCDVk {}

impl AsRef<()> for TestPCDVk
{
    fn as_ref(&self) -> &() { &() }
}

/// A PCD with dummy implementation, it is used just to carry the FinalDarlinDeferredData
/// for testing purposes
pub struct TestPrevPCD<G1: AffineCurve, G2: AffineCurve>(FinalDarlinDeferredData<G1, G2>);

impl<G1, G2> PCD for TestPrevPCD<G1, G2>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    type PCDAccumulator = TestAcc;
    type PCDVerifierKey = TestPCDVk;

    fn succinct_verify(&self, _vk: &Self::PCDVerifierKey) -> Result<<Self::PCDAccumulator as ItemAccumulator>::Item, PCDError> {
        Ok(())
    }
}

pub struct CircuitInfo<G1: AffineCurve, G2: AffineCurve> {
    pub num_constraints: usize,
    pub num_variables:   usize,
    pub dummy_deferred:  FinalDarlinDeferredData<G1, G2>,
}

#[derive(Clone, Default)]
pub struct TestCircuit<G1: AffineCurve, G2: AffineCurve> {
    pub a: Option<G1::ScalarField>,
    pub b: Option<G1::ScalarField>,
    pub c: Option<G1::ScalarField>,
    pub d: Option<G1::ScalarField>,
    pub num_constraints: usize,
    pub num_variables: usize,

    // Deferred elements (sys ins)
    pub deferred: FinalDarlinDeferredData<G1, G2>,
}

impl<G1, G2> ConstraintSynthesizer<G1::ScalarField> for TestCircuit<G1, G2>
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
            || self.c.ok_or(SynthesisError::AssignmentMissing)
        )?;
        let d = cs.alloc_input(
            || "d",
            || self.d.ok_or(SynthesisError::AssignmentMissing)
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

impl<G1, G2> PCDCircuit<G1> for CircuitInfo<G1, G2>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    type IncrementalData =  (G1::ScalarField, G1::ScalarField);
    type SystemInputs    =  FinalDarlinDeferredData<G1, G2>;
    type PreviousPCD     =  TestPrevPCD<G1, G2>;
    type Circuit         =  TestCircuit<G1, G2>;

    fn init(&self) -> Self::Circuit {
        TestCircuit::<G1, G2> {
            a: None,
            b: None,
            c: None,
            d: None,
            num_constraints: self.num_constraints,
            num_variables: self.num_variables,
            deferred: self.dummy_deferred.clone()
        }
    }

    fn init_state(
        &self,
        previous_proofs_data: Vec<Self::PreviousPCD>,
        _previous_proofs_vks: Vec<<Self::PreviousPCD as PCD>::PCDVerifierKey>,
        incremental_data:     Self::IncrementalData
    ) -> Self::Circuit
    {
        assert_eq!(previous_proofs_data.len(), 1);

        let a = incremental_data.0;
        let b = incremental_data.1;

        let mut c = a;
        c.mul_assign(&b);
        let mut d = c;
        d.mul_assign(&b);

        TestCircuit::<G1, G2>{
            a: Some(a),
            b: Some(b),
            c: Some(c),
            d: Some(d),
            num_constraints: self.num_constraints,
            num_variables: self.num_constraints/2,
            deferred: previous_proofs_data[0].0.clone(),
        }
    }

    fn get_sys_ins(circuit: &Self::Circuit) -> Result<Self::SystemInputs, PCDError> {
        Ok(circuit.deferred.clone())
    }

    fn get_usr_ins(circuit: &Self::Circuit) -> Result<Vec<G1::ScalarField>, PCDError> {
        let c = circuit.c.ok_or(PCDError::MissingUserInputs("c".to_owned()))?;
        let d = circuit.d.ok_or(PCDError::MissingUserInputs("d".to_owned()))?;
        Ok(vec![c, d])
    }
}

#[allow(dead_code)]
pub fn generate_test_pcd<'a, G1: AffineCurve, G2:AffineCurve, D: Digest + 'a, R: RngCore>(
    pc_ck_g1: &CommitterKey<G1>,
    final_darlin_pk: &FinalDarlinProverKey<G1::ScalarField, InnerProductArgPC<G1, D>>,
    info: CircuitInfo<G1, G2>,
    zk: bool,
    rng: &mut R,
) -> FinalDarlinPCD<'a, G1, G2, D>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    let prev_pcds = vec![TestPrevPCD::<G1, G2>(info.dummy_deferred.clone())];

    let a = G1::ScalarField::rand(rng);
    let b = G1::ScalarField::rand(rng);

    FinalDarlin::<G1, G2, D>::prove(
        final_darlin_pk,
        pc_ck_g1,
        &info,
        prev_pcds,
        vec![],
        (a, b),
        zk,
        if zk { Some(rng) } else { None }
    ).unwrap()
}

#[allow(dead_code)]
pub fn generate_test_data<'a, G1: AffineCurve, G2: AffineCurve, D: Digest + 'a, R: RngCore>(
    num_constraints: usize,
    segment_size: usize,
    params_g1: &UniversalParams<G1>,
    params_g2: &UniversalParams<G2>,
    num_proofs: usize,
    rng: &mut R,
) -> (
    Vec<FinalDarlinPCD<'a, G1, G2, D>>,
    Vec<FinalDarlinVerifierKey<G1::ScalarField, InnerProductArgPC<G1, D>>>
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
    let dummy_deferred = FinalDarlinDeferredData::<G1, G2>::generate_random::<R, D>(
        rng,
        &committer_key_g1,
        &committer_key_g2,
    );

    let info = CircuitInfo::<G1, G2> {
        num_constraints,
        num_variables: num_constraints/2,
        dummy_deferred,
    };

    let (index_pk, index_vk) = FinalDarlin::<G1, G2, D>::index(
        &committer_key_g1,
        &info
    ).unwrap();

    // Generate Final Darlin PCDs
    let final_darlin_pcd = generate_test_pcd::<G1, G2, D, R>(
        &committer_key_g1,
        &index_pk,
        info,
        rng.gen(),
        rng,
    );

    (vec![final_darlin_pcd; num_proofs], vec![index_vk; num_proofs])
}