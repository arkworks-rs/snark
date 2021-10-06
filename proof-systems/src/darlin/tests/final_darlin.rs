//! A test circuit which, besides processing additional data according to
//! a simple quadratic relation, allocates a given instance of `FinalDarlinDeferredData`,
//! and wires it to the outside via system inputs.
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

/// For testing purposes, TestPrevPCD already serves correct sys_ins and usr_ins
/// to the our test PCDCircuit.
pub struct TestPrevPCD<G1: AffineCurve, G2: AffineCurve> {
    sys_ins: FinalDarlinDeferredData<G1, G2>,
    usr_ins: (G1::ScalarField, G1::ScalarField),
}

impl<G1, G2> PCD for TestPrevPCD<G1, G2>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    type PCDAccumulator = TestAcc;
    type PCDVerifierKey = TestPCDVk;

    // there is nothing to succinctly verify
    fn succinct_verify(&self, _vk: &Self::PCDVerifierKey) -> Result<<Self::PCDAccumulator as ItemAccumulator>::Item, PCDError> {
        Ok(())
    }
}

/// The parameters for our test circuit
#[derive(Clone)]
pub struct CircuitInfo<G1: AffineCurve, G2: AffineCurve> {
    pub num_constraints: usize,
    pub num_variables:   usize,
    /// just used to deduce the number of field elements to allocate on the
    /// circuit for simplicity. Would've been the same passing a parameter
    /// like "number_of_deferred_field_element_to_allocate"
    pub dummy_deferred:  FinalDarlinDeferredData<G1, G2>,
}


/// This test circuit simply allocates `deferred`, i.e. a valid instance of FinalDarlinDeferredData,
/// and wires it to the outside via system inputs.
/// The user inputs are the field elements c, d, used along the user inputs (c_prev, d_prev)
/// of the previous proof in order to satisfy:
///     (c,d) = (a * b - c_prev, a * b^2 - d_prev),
/// using the prepared a,b. To produce the given `num_constraints`, the same two constraints
///     a * b = c
///     c * b = d
/// are repeated accordingly. To acchieve the give `num_variables`, a corresponding number of
/// dummy witness variables are allocated.
#[derive(Clone, Default)]
pub struct TestCircuit<G1: AffineCurve, G2: AffineCurve> {

    /// Incremental data (to be allocated as witnesses)
    pub a: Option<G1::ScalarField>,
    pub b: Option<G1::ScalarField>,

    /// Previous user inputs (to be allocated as witnesses)
    pub c_prev: Option<G1::ScalarField>,
    pub d_prev: Option<G1::ScalarField>,

    /// Actual user inputs (to be allocated as public inputs)
    pub c: Option<G1::ScalarField>,
    pub d: Option<G1::ScalarField>,

    // System inputs (i.e previous accumulators, to be allocated as public inputs)
    pub deferred: FinalDarlinDeferredData<G1, G2>,

    /// Setup data
    pub num_constraints: usize,
    pub num_variables: usize,
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
        // convert the FinalDarlinDeferred efficiently to circuit inputs
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

        // Enforce the system inputs to the circuit to be equal to the allocated `deferred`.
        // This is a simple way to allow test cases where sys data (i.e. the deferred 
        // accumulators) are wrong.
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
        let c_prev = cs.alloc(|| "c_prev", || self.c_prev.ok_or(SynthesisError::AssignmentMissing))?;
        let d_prev = cs.alloc(|| "d_prev", || self.d_prev.ok_or(SynthesisError::AssignmentMissing))?;
        let c = cs.alloc_input(
            || "c",
            || self.c.ok_or(SynthesisError::AssignmentMissing)
        )?;
        let d = cs.alloc_input(
            || "d",
            || self.d.ok_or(SynthesisError::AssignmentMissing)
        )?;

        // TODO: This calculation is wrong, as enforce_equal allocates new variables.
        //       However, fixing this may cause unit tests to crash since num_constraints
        //       and num_variables are generated at random and an underflow may happen.
        //       Fix both.
        for i in 0..(self.num_variables - 7 - (2 * deferred_len)) {
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
                |lc| lc + c_prev + c,
            );
        }
        cs.enforce(
            || format!("constraint {}", self.num_constraints - 1),
            |lc| lc + c,
            |lc| lc + b,
            |lc| lc + d_prev + d,
        );

        Ok(())
    }
}

impl<G1, G2> PCDCircuit<G1> for TestCircuit<G1, G2>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
{
    type SetupData       =  CircuitInfo<G1, G2>;
    type AdditionalData =  (G1::ScalarField, G1::ScalarField);
    type SystemInputs    =  FinalDarlinDeferredData<G1, G2>;
    type PreviousPCD     =  TestPrevPCD<G1, G2>;

    fn init(config: Self::SetupData) -> Self {
        Self {
            a: None,
            b: None,
            c_prev: None,
            d_prev: None,
            c: None,
            d: None,
            num_constraints: config.num_constraints,
            num_variables: config.num_variables,
            deferred: config.dummy_deferred.clone()
        }
    }

    fn init_state(
        config:               Self::SetupData,
        previous_proofs_data: Vec<Self::PreviousPCD>,
        _previous_proofs_vks: Vec<<Self::PreviousPCD as PCD>::PCDVerifierKey>,
        additional_data:     Self::AdditionalData
    ) -> Self
    {
        assert_eq!(previous_proofs_data.len(), 1);

        let a = additional_data.0;
        let b = additional_data.1;
        let c_prev = previous_proofs_data[0].usr_ins.0;
        let d_prev = previous_proofs_data[0].usr_ins.1;

        let c = (a * &b) - &c_prev;
        let d = (c * &b) - &d_prev;

        Self {
            a: Some(a),
            b: Some(b),
            c_prev: Some(c_prev),
            d_prev: Some(d_prev),
            c: Some(c),
            d: Some(d),
            num_constraints: config.num_constraints,
            num_variables: config.num_variables,
            deferred: previous_proofs_data[0].sys_ins.clone(),
        }
    }

    fn get_sys_ins(&self) -> Result<&Self::SystemInputs, PCDError> {
        Ok(&self.deferred)
    }

    fn get_usr_ins(&self) -> Result<Vec<G1::ScalarField>, PCDError> {
        let c = self.c.ok_or(PCDError::MissingUserInputs("c".to_owned()))?;
        let d = self.d.ok_or(PCDError::MissingUserInputs("d".to_owned()))?;
        Ok(vec![c, d])
    }
}

/// Generates a FinalDarlinPCD from TestCircuit1, given an instance of 
/// FinalDarlinDeferred as previous PCD (via CircuitInfo). 
/// The additional data a,b is sampled randomly.
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

    let prev_pcd = TestPrevPCD::<G1, G2> {
        // as we have already generated a dummy deferred for CircuitInfo, let's
        // just re-use it
        sys_ins: info.dummy_deferred.clone(),
        usr_ins: (G1::ScalarField::rand(rng), G1::ScalarField::rand(rng))
    };

    // our additional data witnesses
    let a = G1::ScalarField::rand(rng);
    let b = G1::ScalarField::rand(rng);

    FinalDarlin::<G1, G2, D>::prove::<TestCircuit<G1, G2>>(
        final_darlin_pk,
        pc_ck_g1,
        info,
        vec![prev_pcd],
        vec![],
        (a, b),
        zk,
        if zk { Some(rng) } else { None }
    ).unwrap()
}

/// Generates `num_proofs` random instances of FinalDarlinPCDs for TestCircuit1 at given 
/// `num_constraints`, using `segment_size` for the dlog commitment scheme.
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
        num_variables: num_constraints,
        dummy_deferred,
    };

    let (index_pk, index_vk) = FinalDarlin::<G1, G2, D>::index::<TestCircuit<G1, G2>>(
        &committer_key_g1,
        info.clone()
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