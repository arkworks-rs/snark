use algebra::{fields::FpParameters, BigInteger, Field, PrimeField};
use algebra_core::{PairingEngine, ToConstraintField};
use crypto_primitives::nizk::{
    constraints::NIZKVerifierGadget,
    groth16::{
        constraints::{Groth16VerifierGadget, ProofGadget, VerifyingKeyGadget},
        Groth16,
    },
};
use groth16::{Parameters, Proof};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use r1cs_std::{
    alloc::AllocGadget, bits::ToBitsGadget, boolean::Boolean, fields::fp::FpGadget,
    pairing::PairingGadget as PG, uint8::UInt8,
};
use std::marker::PhantomData;

pub trait CurvePair {
    type PairingEngineTick: PairingEngine;
    type PairingEngineTock: PairingEngine;
    type PairingGadgetTick: PG<
        Self::PairingEngineTick,
        <Self::PairingEngineTock as PairingEngine>::Fr,
    >;
    type PairingGadgetTock: PG<
        Self::PairingEngineTock,
        <Self::PairingEngineTick as PairingEngine>::Fr,
    >;

    const TICK_CURVE: &'static str;
    const TOCK_CURVE: &'static str;
}

// Verifying InnerCircuit in MiddleCircuit
type InnerProofSystem<C> = Groth16<
    <C as CurvePair>::PairingEngineTick,
    InnerCircuit<<<C as CurvePair>::PairingEngineTick as PairingEngine>::Fr>,
    <<C as CurvePair>::PairingEngineTick as PairingEngine>::Fr,
>;
type InnerVerifierGadget<C> = Groth16VerifierGadget<
    <C as CurvePair>::PairingEngineTick,
    <<C as CurvePair>::PairingEngineTock as PairingEngine>::Fr,
    <C as CurvePair>::PairingGadgetTick,
>;
type InnerProofGadget<C> = ProofGadget<
    <C as CurvePair>::PairingEngineTick,
    <<C as CurvePair>::PairingEngineTock as PairingEngine>::Fr,
    <C as CurvePair>::PairingGadgetTick,
>;
type InnerVkGadget<C> = VerifyingKeyGadget<
    <C as CurvePair>::PairingEngineTick,
    <<C as CurvePair>::PairingEngineTock as PairingEngine>::Fr,
    <C as CurvePair>::PairingGadgetTick,
>;

// Verifying MiddleCircuit in OuterCircuit
type MiddleProofSystem<C> = Groth16<
    <C as CurvePair>::PairingEngineTock,
    MiddleCircuit<C>,
    <<C as CurvePair>::PairingEngineTock as PairingEngine>::Fr,
>;
type MiddleVerifierGadget<C> = Groth16VerifierGadget<
    <C as CurvePair>::PairingEngineTock,
    <<C as CurvePair>::PairingEngineTick as PairingEngine>::Fr,
    <C as CurvePair>::PairingGadgetTock,
>;
type MiddleProofGadget<C> = ProofGadget<
    <C as CurvePair>::PairingEngineTock,
    <<C as CurvePair>::PairingEngineTick as PairingEngine>::Fr,
    <C as CurvePair>::PairingGadgetTock,
>;
type MiddleVkGadget<C> = VerifyingKeyGadget<
    <C as CurvePair>::PairingEngineTock,
    <<C as CurvePair>::PairingEngineTick as PairingEngine>::Fr,
    <C as CurvePair>::PairingGadgetTock,
>;

pub struct InnerCircuit<F: Field> {
    num_constraints: usize,
    inputs: Vec<F>,
}

impl<F: Field> InnerCircuit<F> {
    pub fn new(num_constraints: usize, inputs: Vec<F>) -> Self {
        Self {
            num_constraints,
            inputs,
        }
    }
}

impl<F: Field> ConstraintSynthesizer<F> for InnerCircuit<F> {
    fn generate_constraints<CS: ConstraintSystem<F>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
        assert!(self.inputs.len() >= 2);
        assert!(self.num_constraints >= self.inputs.len());

        let mut variables: Vec<_> = Vec::with_capacity(self.inputs.len());
        for (i, input) in self.inputs.into_iter().enumerate() {
            let input_var = cs.alloc_input(|| format!("Input {}", i), || Ok(input))?;
            variables.push((input, input_var));
        }

        for i in 0..self.num_constraints {
            let new_entry = {
                let (input_1_val, input_1_var) = variables[i];
                let (input_2_val, input_2_var) = variables[i + 1];
                let result_val = input_1_val * input_2_val;
                let result_var = cs.alloc(|| format!("Result {}", i), || Ok(result_val))?;
                cs.enforce(
                    || format!("Enforce constraint {}", i),
                    |lc| lc + input_1_var,
                    |lc| lc + input_2_var,
                    |lc| lc + result_var,
                );
                (result_val, result_var)
            };
            variables.push(new_entry);
        }
        Ok(())
    }
}

pub struct MiddleCircuit<C: CurvePair> {
    inputs: Vec<<C::PairingEngineTick as PairingEngine>::Fr>,
    params: Parameters<C::PairingEngineTick>,
    proof: Proof<C::PairingEngineTick>,
    _curve_pair: PhantomData<C>,
}

impl<C: CurvePair> MiddleCircuit<C> {
    pub fn new(
        inputs: Vec<<C::PairingEngineTick as PairingEngine>::Fr>,
        params: Parameters<C::PairingEngineTick>,
        proof: Proof<C::PairingEngineTick>,
    ) -> Self {
        Self {
            inputs,
            params,
            proof,
            _curve_pair: PhantomData,
        }
    }

    pub fn inputs(
        inputs: &[<C::PairingEngineTick as PairingEngine>::Fr],
    ) -> Vec<<C::PairingEngineTock as PairingEngine>::Fr> {
        let input_bytes = inputs
            .iter()
            .flat_map(|input| {
                input
                    .into_repr()
                    .as_ref()
                    .iter()
                    .flat_map(|l| l.to_le_bytes().to_vec())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        input_bytes[..].to_field_elements().unwrap()
    }
}

impl<C: CurvePair> ConstraintSynthesizer<<C::PairingEngineTock as PairingEngine>::Fr>
    for MiddleCircuit<C>
{
    fn generate_constraints<CS: ConstraintSystem<<C::PairingEngineTock as PairingEngine>::Fr>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
        let params = self.params;
        let proof = self.proof;
        let inputs = self.inputs;
        let input_gadgets;

        {
            let mut cs = cs.ns(|| "Allocate Input");
            // Chain all input values in one large byte array.
            let input_bytes = inputs
                .into_iter()
                .flat_map(|input| {
                    input
                        .into_repr()
                        .as_ref()
                        .iter()
                        .flat_map(|l| l.to_le_bytes().to_vec())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();

            // Allocate this byte array as input packed into field elements.
            let input_bytes = UInt8::alloc_input_vec(cs.ns(|| "Input"), &input_bytes[..])?;
            // 40 byte
            let element_size =
                <<<C::PairingEngineTick as PairingEngine>::Fr as PrimeField>::Params as FpParameters>::BigInt::NUM_LIMBS * 8;
            input_gadgets = input_bytes
                .chunks(element_size)
                .map(|chunk| {
                    chunk
                        .iter()
                        .flat_map(|byte| byte.into_bits_le())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
        }
        println!("|---- Num inputs for sub-SNARK: {}", input_gadgets.len());
        let num_constraints = cs.num_constraints();
        println!(
            "|---- Num constraints to prepare inputs: {}",
            num_constraints
        );

        let vk_gadget = InnerVkGadget::<C>::alloc(cs.ns(|| "Vk"), || Ok(&params.vk))?;
        let proof_gadget = InnerProofGadget::<C>::alloc(cs.ns(|| "Proof"), || Ok(proof.clone()))?;
        <InnerVerifierGadget<C> as NIZKVerifierGadget<
            InnerProofSystem<C>,
            <C::PairingEngineTock as PairingEngine>::Fr,
        >>::check_verify(
            cs.ns(|| "Verify"),
            &vk_gadget,
            input_gadgets.iter(),
            &proof_gadget,
        )?;
        println!(
            "|---- Num constraints for sub-SNARK verification: {}",
            cs.num_constraints() - num_constraints
        );
        Ok(())
    }
}

pub struct OuterCircuit<C: CurvePair> {
    inputs: Vec<<C::PairingEngineTick as PairingEngine>::Fr>,
    params: Parameters<C::PairingEngineTock>,
    proof: Proof<C::PairingEngineTock>,
    _curve_pair: PhantomData<C>,
}

impl<C: CurvePair> OuterCircuit<C> {
    pub fn new(
        inputs: Vec<<C::PairingEngineTick as PairingEngine>::Fr>,
        params: Parameters<C::PairingEngineTock>,
        proof: Proof<C::PairingEngineTock>,
    ) -> Self {
        Self {
            inputs,
            params,
            proof,
            _curve_pair: PhantomData,
        }
    }
}

impl<C: CurvePair> ConstraintSynthesizer<<C::PairingEngineTick as PairingEngine>::Fr>
    for OuterCircuit<C>
{
    fn generate_constraints<CS: ConstraintSystem<<C::PairingEngineTick as PairingEngine>::Fr>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
        let params = self.params;
        let proof = self.proof;
        let inputs = self.inputs;
        let mut input_gadgets = Vec::new();

        {
            let bigint_size =
                <<<C::PairingEngineTick as PairingEngine>::Fr as PrimeField>::Params as FpParameters>::BigInt::NUM_LIMBS * 64;
            let mut input_bits = Vec::new();
            let mut cs = cs.ns(|| "Allocate Input");
            for (i, input) in inputs.into_iter().enumerate() {
                let input_gadget =
                    FpGadget::alloc_input(cs.ns(|| format!("Input {}", i)), || Ok(input))?;
                let mut fp_bits = input_gadget.to_bits(cs.ns(|| format!("To bits {}", i)))?;

                // FpGadget::to_bits outputs a big-endian binary representation of
                // fe_gadget's value, so we have to reverse it to get the little-endian
                // form.
                fp_bits.reverse();

                // Use 320 bits per element.
                for _ in fp_bits.len()..bigint_size {
                    fp_bits.push(Boolean::constant(false));
                }
                input_bits.extend_from_slice(&fp_bits);
            }

            // Pack input bits into field elements of the underlying circuit.
            let max_size = 8
                * (<<<C::PairingEngineTock as PairingEngine>::Fr as PrimeField>::Params as FpParameters>::CAPACITY / 8)
                    as usize;
            let bigint_size =
                <<<C::PairingEngineTock as PairingEngine>::Fr as PrimeField>::Params as FpParameters>::BigInt::NUM_LIMBS * 64;
            for chunk in input_bits.chunks(max_size) {
                let mut chunk = chunk.to_vec();
                let len = chunk.len();
                for _ in len..bigint_size {
                    chunk.push(Boolean::constant(false));
                }
                input_gadgets.push(chunk);
            }
        }
        println!("|---- Num inputs for sub-SNARK: {}", input_gadgets.len());
        let num_constraints = cs.num_constraints();
        println!(
            "|---- Num constraints to prepare inputs: {}",
            num_constraints
        );

        let vk_gadget = MiddleVkGadget::<C>::alloc(cs.ns(|| "Vk"), || Ok(&params.vk))?;
        let proof_gadget = MiddleProofGadget::<C>::alloc(cs.ns(|| "Proof"), || Ok(proof.clone()))?;
        <MiddleVerifierGadget<C> as NIZKVerifierGadget<
            MiddleProofSystem<C>,
            <C::PairingEngineTick as PairingEngine>::Fr,
        >>::check_verify(
            cs.ns(|| "Verify"),
            &vk_gadget,
            input_gadgets.iter(),
            &proof_gadget,
        )?;
        println!(
            "|---- Num constraints for sub-SNARK verification: {}",
            cs.num_constraints() - num_constraints
        );
        Ok(())
    }
}
