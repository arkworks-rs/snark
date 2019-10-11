use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{
    commitment::blake2s::Blake2sCommitment,
    prf::blake2s::constraints::{blake2s_gadget, Blake2sOutputGadget},
    CommitmentGadget,
};
use algebra::{PrimeField, Field};
use r1cs_std::prelude::*;

use std::borrow::Borrow;

#[derive(Clone)]
pub struct Blake2sParametersGadget;

#[derive(Clone)]
pub struct Blake2sRandomnessGadget(pub Vec<UInt8>);

pub struct Blake2sCommitmentGadget;

impl<ConstraintF: PrimeField> CommitmentGadget<Blake2sCommitment, ConstraintF> for Blake2sCommitmentGadget {
    type OutputGadget = Blake2sOutputGadget;
    type ParametersGadget = Blake2sParametersGadget;
    type RandomnessGadget = Blake2sRandomnessGadget;

    fn check_commitment_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        _: &Self::ParametersGadget,
        input: &[UInt8],
        r: &Self::RandomnessGadget,
    ) -> Result<Self::OutputGadget, SynthesisError> {
        let mut input_bits = Vec::with_capacity(512);
        for byte in input.into_iter().chain(r.0.iter()) {
            input_bits.extend_from_slice(&byte.into_bits_le());
        }
        let mut result = Vec::new();
        for (i, int) in blake2s_gadget(cs.ns(|| "Blake2s Eval"), &input_bits)?
            .into_iter()
            .enumerate()
        {
            let chunk = int.to_bytes(&mut cs.ns(|| format!("Result ToBytes {}", i)))?;
            result.extend_from_slice(&chunk);
        }
        Ok(Blake2sOutputGadget(result))
    }
}

impl<ConstraintF: Field> AllocGadget<(), ConstraintF> for Blake2sParametersGadget {
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(_: CS, _: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<()>,
    {
        Ok(Blake2sParametersGadget)
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(_: CS, _: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<()>,
    {
        Ok(Blake2sParametersGadget)
    }
}

impl<ConstraintF: PrimeField> AllocGadget<[u8; 32], ConstraintF> for Blake2sRandomnessGadget {
    #[inline]
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, value_gen: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[u8; 32]>,
    {
        let zeros = [0u8; 32];
        let value = match value_gen() {
            Ok(val) => *(val.borrow()),
            Err(_) => zeros,
        };
        let bytes = <UInt8>::alloc_vec(cs, &value)?;

        Ok(Blake2sRandomnessGadget(bytes))
    }

    #[inline]
    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[u8; 32]>,
    {
        let zeros = [0u8; 32];
        let value = match value_gen() {
            Ok(val) => *(val.borrow()),
            Err(_) => zeros,
        };
        let bytes = <UInt8>::alloc_input_vec(cs, &value)?;

        Ok(Blake2sRandomnessGadget(bytes))
    }
}

#[cfg(test)]
mod test {
    use algebra::fields::bls12_381::Fr;
    use rand::{thread_rng, Rng};

    use crate::*;
    use crate::{
        commitment::blake2s::Blake2sCommitment,
        commitment::blake2s::constraints::{Blake2sCommitmentGadget, Blake2sRandomnessGadget},
    };
    use r1cs_core::ConstraintSystem;
    use r1cs_std::prelude::*;
    use r1cs_std::test_constraint_system::TestConstraintSystem;

    #[test]
    fn commitment_gadget_test() {
        let mut cs = TestConstraintSystem::<Fr>::new();

        let input = [1u8; 32];

        let rng = &mut thread_rng();

        type TestCOMM = Blake2sCommitment;
        type TestCOMMGadget = Blake2sCommitmentGadget;

        let mut randomness = [0u8; 32];
        rng.fill(&mut randomness);

        let parameters = ();
        let primitive_result = Blake2sCommitment::commit(&parameters, &input, &randomness).unwrap();

        let mut input_bytes = vec![];
        for (byte_i, input_byte) in input.into_iter().enumerate() {
            let cs = cs.ns(|| format!("input_byte_gadget_{}", byte_i));
            input_bytes.push(UInt8::alloc(cs, || Ok(*input_byte)).unwrap());
        }

        let mut randomness_bytes = vec![];
        for (byte_i, random_byte) in randomness.into_iter().enumerate() {
            let cs = cs.ns(|| format!("randomness_byte_gadget_{}", byte_i));
            randomness_bytes.push(UInt8::alloc(cs, || Ok(*random_byte)).unwrap());
        }
        let randomness_bytes = Blake2sRandomnessGadget(randomness_bytes);

        let gadget_parameters =
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Fr>>::ParametersGadget::alloc(
                &mut cs.ns(|| "gadget_parameters"),
                || Ok(&parameters),
            )
            .unwrap();
        let gadget_result =
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Fr>>::check_commitment_gadget(
                &mut cs.ns(|| "gadget_evaluation"),
                &gadget_parameters,
                &input_bytes,
                &randomness_bytes,
            )
            .unwrap();

        for i in 0..32 {
            assert_eq!(primitive_result[i], gadget_result.0[i].get_value().unwrap());
        }
        assert!(cs.is_satisfied());
    }
}
