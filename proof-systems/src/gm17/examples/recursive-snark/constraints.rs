// This example uses Groth17 over Coda's MNT cycle to wrap a base circuit (the "inner circuit") of
// specified number of inputs and constraints twice. See the description of the Groth16 example for
// details.
use algebra::{fields::FpParameters, Field, PairingEngine, PrimeField, ToBits, ToConstraintField};

use proof_systems::gm17::{Parameters, Proof};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use r1cs_crypto::nizk::{
    gm17::{Gm17, Gm17VerifierGadget, ProofGadget, VerifyingKeyGadget},
    NIZKVerifierGadget,
};
use r1cs_std::{
    alloc::AllocGadget, bits::ToBitsGadget, boolean::Boolean, fields::fp::FpGadget,
    pairing::PairingGadget as PG,
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
type InnerProofSystem<C> = Gm17<
    <C as CurvePair>::PairingEngineTick,
    InnerCircuit<<<C as CurvePair>::PairingEngineTick as PairingEngine>::Fr>,
    <<C as CurvePair>::PairingEngineTick as PairingEngine>::Fr,
>;
type InnerVerifierGadget<C> = Gm17VerifierGadget<
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
type MiddleProofSystem<C> = Gm17<
    <C as CurvePair>::PairingEngineTock,
    MiddleCircuit<C>,
    <<C as CurvePair>::PairingEngineTock as PairingEngine>::Fr,
>;
type MiddleVerifierGadget<C> = Gm17VerifierGadget<
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
                let result_val = input_1_val * &input_2_val;
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
        let input_bits = inputs
            .iter()
            .flat_map(|input| input.write_bits())
            .collect::<Vec<_>>();

        input_bits[..].to_field_elements().unwrap()
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
            // Chain all input values in one large bit array.
            let input_bits = inputs
                .into_iter()
                .flat_map(|input| input.write_bits())
                .collect::<Vec<_>>();

            // Allocate this bit array as input packed into field elements.
            let input_bits = Boolean::alloc_input_vec(cs.ns(|| "Input"), &input_bits[..])?;

            let element_size =
                <<C::PairingEngineTick as PairingEngine>::Fr as PrimeField>::Params::MODULUS_BITS;
            input_gadgets = input_bits
                .chunks(element_size as usize)
                .map(|chunk| {
                    let mut chunk = chunk.to_vec();
                    chunk.reverse();
                    chunk
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
            let mut input_bits = Vec::new();
            let mut cs = cs.ns(|| "Allocate Input");
            for (i, input) in inputs.into_iter().enumerate() {
                let input_gadget =
                    FpGadget::alloc_input(cs.ns(|| format!("Input {}", i)), || Ok(input))?;
                let fp_bits = input_gadget.to_bits(cs.ns(|| format!("To bits {}", i)))?;
                input_bits.extend_from_slice(&fp_bits);
            }

            // Pack input bits into field elements of the underlying circuit.
            let max_size = <<<C::PairingEngineTock as PairingEngine>::Fr as PrimeField>::Params as FpParameters>::CAPACITY as usize;
            let bigint_size = (<<C::PairingEngineTock as PairingEngine>::Fr as PrimeField>::Params::MODULUS_BITS +
                <<C::PairingEngineTock as PairingEngine>::Fr as PrimeField>::Params::REPR_SHAVE_BITS) as usize;
            for chunk in input_bits.chunks(max_size) {
                let len = chunk.len();
                let mut padding = vec![Boolean::constant(false); bigint_size - len];
                padding.extend_from_slice(chunk);
                padding.reverse();
                input_gadgets.push(padding);
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
