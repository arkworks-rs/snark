use crate::crypto_primitives::commitment::blake2s::Blake2sCommitment;
use snark::{ConstraintSystem, SynthesisError};

use crate::gadgets::{
    prf::blake2s::{blake2s_gadget, Blake2sOutputGadget},
    CommitmentGadget,
};
use algebra::PairingEngine;
use snark_gadgets::{
    uint8::UInt8,
    utils::{AllocGadget, ToBytesGadget},
};
use std::borrow::Borrow;

#[derive(Clone)]
pub struct Blake2sParametersGadget;

#[derive(Clone)]
pub struct Blake2sRandomnessGadget(pub Vec<UInt8>);

pub struct Blake2sCommitmentGadget;

impl<E: PairingEngine> CommitmentGadget<Blake2sCommitment, E> for Blake2sCommitmentGadget {
    type OutputGadget = Blake2sOutputGadget;
    type ParametersGadget = Blake2sParametersGadget;
    type RandomnessGadget = Blake2sRandomnessGadget;

    fn check_commitment_gadget<CS: ConstraintSystem<E>>(
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

impl<E: PairingEngine> AllocGadget<(), E> for Blake2sParametersGadget {
    fn alloc<F, T, CS: ConstraintSystem<E>>(_: CS, _: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<()>,
    {
        Ok(Blake2sParametersGadget)
    }

    fn alloc_input<F, T, CS: ConstraintSystem<E>>(_: CS, _: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<()>,
    {
        Ok(Blake2sParametersGadget)
    }
}

impl<E: PairingEngine> AllocGadget<[u8; 32], E> for Blake2sRandomnessGadget {
    #[inline]
    fn alloc<F, T, CS: ConstraintSystem<E>>(cs: CS, value_gen: F) -> Result<Self, SynthesisError>
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
    fn alloc_input<F, T, CS: ConstraintSystem<E>>(
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
    use algebra::curves::bls12_381::Bls12_381;
    use rand::{thread_rng, Rng};

    use crate::{
        crypto_primitives::commitment::{blake2s::Blake2sCommitment, CommitmentScheme},
        gadgets::commitment::{
            blake2s::{Blake2sCommitmentGadget, Blake2sRandomnessGadget},
            CommitmentGadget,
        },
    };
    use snark::ConstraintSystem;
    use snark_gadgets::{
        test_constraint_system::TestConstraintSystem, uint8::UInt8, utils::AllocGadget,
    };

    #[test]
    fn commitment_gadget_test() {
        let mut cs = TestConstraintSystem::<Bls12_381>::new();

        let input = [1u8; 32];

        let rng = &mut thread_rng();

        type TestCOMM = Blake2sCommitment;
        type TestCOMMGadget = Blake2sCommitmentGadget;

        let mut randomness = [0u8; 32];
        rng.fill_bytes(&mut randomness);

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
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Bls12_381>>::ParametersGadget::alloc(
                &mut cs.ns(|| "gadget_parameters"),
                || Ok(&parameters),
            )
            .unwrap();
        let gadget_result =
            <TestCOMMGadget as CommitmentGadget<TestCOMM, Bls12_381>>::check_commitment_gadget(
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
