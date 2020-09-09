use r1cs_core::{Namespace, SynthesisError};

use crate::{
    commitment::blake2s,
    commitment::CommitmentGadget,
    prf::blake2s::constraints::{evaluate_blake2s, OutputVar},
    Vec,
};
use algebra_core::{Field, PrimeField};
use r1cs_std::prelude::*;

use core::borrow::Borrow;

#[derive(Clone)]
pub struct ParametersVar;

#[derive(Clone)]
pub struct RandomnessVar<F: Field>(pub Vec<UInt8<F>>);

pub struct CommGadget;

impl<F: PrimeField> CommitmentGadget<blake2s::Commitment, F> for CommGadget {
    type OutputVar = OutputVar<F>;
    type ParametersVar = ParametersVar;
    type RandomnessVar = RandomnessVar<F>;

    #[tracing::instrument(target = "r1cs", skip(input, r))]
    fn commit(
        _: &Self::ParametersVar,
        input: &[UInt8<F>],
        r: &Self::RandomnessVar,
    ) -> Result<Self::OutputVar, SynthesisError> {
        let mut input_bits = Vec::with_capacity(512);
        for byte in input.iter().chain(r.0.iter()) {
            input_bits.extend_from_slice(&byte.to_bits_le()?);
        }
        let mut result = Vec::new();
        for int in evaluate_blake2s(&input_bits)?.into_iter() {
            let chunk = int.to_bytes()?;
            result.extend_from_slice(&chunk);
        }
        Ok(OutputVar(result))
    }
}

impl<ConstraintF: Field> AllocVar<(), ConstraintF> for ParametersVar {
    #[tracing::instrument(target = "r1cs", skip(_cs, _f))]
    fn new_variable<T: Borrow<()>>(
        _cs: impl Into<Namespace<ConstraintF>>,
        _f: impl FnOnce() -> Result<T, SynthesisError>,
        _mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        Ok(ParametersVar)
    }
}

impl<ConstraintF: PrimeField> AllocVar<[u8; 32], ConstraintF> for RandomnessVar<ConstraintF> {
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T: Borrow<[u8; 32]>>(
        cs: impl Into<Namespace<ConstraintF>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let bytes = f().map(|b| *b.borrow()).unwrap_or([0u8; 32]);
        match mode {
            AllocationMode::Constant => Ok(Self(UInt8::constant_vec(&bytes))),
            AllocationMode::Input => UInt8::new_input_vec(cs, &bytes).map(Self),
            AllocationMode::Witness => UInt8::new_witness_vec(cs, &bytes).map(Self),
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{
        commitment::blake2s::{
            constraints::{CommGadget, RandomnessVar},
            Commitment,
        },
        commitment::{CommitmentGadget, CommitmentScheme},
    };
    use algebra::{ed_on_bls12_381::Fq as Fr, test_rng};
    use r1cs_core::ConstraintSystem;
    use r1cs_std::prelude::*;
    use rand::Rng;

    #[test]
    fn commitment_gadget_test() {
        let cs = ConstraintSystem::<Fr>::new_ref();

        let input = [1u8; 32];

        let rng = &mut test_rng();

        type TestCOMM = Commitment;
        type TestCOMMGadget = CommGadget;

        let mut randomness = [0u8; 32];
        rng.fill(&mut randomness);

        let parameters = ();
        let primitive_result = Commitment::commit(&parameters, &input, &randomness).unwrap();

        let mut input_var = vec![];
        for byte in &input {
            input_var.push(UInt8::new_witness(cs.clone(), || Ok(*byte)).unwrap());
        }

        let mut randomness_var = vec![];
        for r_byte in randomness.iter() {
            randomness_var.push(UInt8::new_witness(cs.clone(), || Ok(r_byte)).unwrap());
        }
        let randomness_var = RandomnessVar(randomness_var);

        let parameters_var =
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Fr>>::ParametersVar::new_witness(
                r1cs_core::ns!(cs, "gadget_parameters"),
                || Ok(&parameters),
            )
            .unwrap();
        let result_var = <TestCOMMGadget as CommitmentGadget<TestCOMM, Fr>>::commit(
            &parameters_var,
            &input_var,
            &randomness_var,
        )
        .unwrap();

        for i in 0..32 {
            assert_eq!(primitive_result[i], result_var.0[i].value().unwrap());
        }
        assert!(cs.is_satisfied().unwrap());
    }
}
