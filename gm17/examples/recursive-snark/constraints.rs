use algebra::{
    fields::{FftParameters, FpParameters},
    BigInteger, Field, PrimeField,
};
use algebra_core::{PairingEngine, ToConstraintField};
use core::ops::MulAssign;
use crypto_primitives::nizk::{
    constraints::NIZKVerifierGadget,
    gm17::{
        constraints::{Gm17VerifierGadget, ProofVar, VerifyingKeyVar},
        Gm17,
    },
};
use gm17::{Parameters, Proof};
use r1cs_core::{lc, ConstraintSynthesizer, ConstraintSystemRef, SynthesisError};
use r1cs_std::{fields::fp::FpVar, pairing::PairingVar as PG, prelude::*};
use std::marker::PhantomData;

pub trait CurvePair
where
    <Self::TickGroup as PairingEngine>::G1Projective:
        MulAssign<<Self::TockGroup as PairingEngine>::Fq>,
    <Self::TickGroup as PairingEngine>::G2Projective:
        MulAssign<<Self::TockGroup as PairingEngine>::Fq>,
    <Self::TickGroup as PairingEngine>::G1Affine:
        ToConstraintField<<<Self::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
    <Self::TickGroup as PairingEngine>::G2Affine:
        ToConstraintField<<<Self::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
{
    type TickGroup: PairingEngine<
        Fq = <Self::TockGroup as PairingEngine>::Fr,
        Fr = <Self::TockGroup as PairingEngine>::Fq,
    >;
    type TockGroup: PairingEngine;

    const TICK_CURVE: &'static str;
    const TOCK_CURVE: &'static str;
}

// Verifying InnerCircuit in MiddleCircuit
type InnerProofSystem<C> = Gm17<
    <C as CurvePair>::TickGroup,
    InnerCircuit<<<C as CurvePair>::TickGroup as PairingEngine>::Fr>,
    <<C as CurvePair>::TickGroup as PairingEngine>::Fr,
>;

type InnerVerifierGadget<C, PV> = Gm17VerifierGadget<<C as CurvePair>::TickGroup, PV>;
type InnerProofVar<C, PV> = ProofVar<<C as CurvePair>::TickGroup, PV>;
type InnerVkVar<C, PV> = VerifyingKeyVar<<C as CurvePair>::TickGroup, PV>;

// Verifying MiddleCircuit in OuterCircuit
type MiddleProofSystem<C, PV> = Gm17<
    <C as CurvePair>::TockGroup,
    MiddleCircuit<C, PV>,
    <<C as CurvePair>::TockGroup as PairingEngine>::Fr,
>;
type MiddleVerifierGadget<C, PV> = Gm17VerifierGadget<<C as CurvePair>::TockGroup, PV>;
type MiddleProofVar<C, PV> = ProofVar<<C as CurvePair>::TockGroup, PV>;
type MiddleVkVar<C, PV> = VerifyingKeyVar<<C as CurvePair>::TockGroup, PV>;

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
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        assert!(self.inputs.len() >= 2);
        assert!(self.num_constraints >= self.inputs.len());

        let mut variables: Vec<_> = Vec::with_capacity(self.inputs.len());
        for input in self.inputs.into_iter() {
            let input_var = cs.new_input_variable(|| Ok(input))?;
            variables.push((input, input_var));
        }

        for i in 0..self.num_constraints {
            let new_entry = {
                let (input_1_val, input_1_var) = variables[i];
                let (input_2_val, input_2_var) = variables[i + 1];
                let result_val = input_1_val * input_2_val;
                let result_var = cs.new_witness_variable(|| Ok(result_val))?;
                cs.enforce_constraint(
                    lc!() + input_1_var,
                    lc!() + input_2_var,
                    lc!() + result_var,
                )?;
                (result_val, result_var)
            };
            variables.push(new_entry);
        }
        Ok(())
    }
}

pub struct MiddleCircuit<C: CurvePair, TickPairing: PG<C::TickGroup>>
where
    <C::TickGroup as PairingEngine>::G1Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G2Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G1Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
    <C::TickGroup as PairingEngine>::G2Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
{
    inputs: Vec<<C::TickGroup as PairingEngine>::Fr>,
    params: Parameters<C::TickGroup>,
    proof: Proof<C::TickGroup>,
    _curve_pair: PhantomData<C>,
    _tick_pairing: PhantomData<TickPairing>,
}

impl<C: CurvePair, TickPairing: PG<C::TickGroup>> MiddleCircuit<C, TickPairing>
where
    <C::TickGroup as PairingEngine>::G1Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G2Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G1Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
    <C::TickGroup as PairingEngine>::G2Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
{
    pub fn new(
        inputs: Vec<<C::TickGroup as PairingEngine>::Fr>,
        params: Parameters<C::TickGroup>,
        proof: Proof<C::TickGroup>,
    ) -> Self {
        Self {
            inputs,
            params,
            proof,
            _curve_pair: PhantomData,
            _tick_pairing: PhantomData,
        }
    }

    pub fn inputs(
        inputs: &[<C::TickGroup as PairingEngine>::Fr],
    ) -> Vec<<C::TockGroup as PairingEngine>::Fr> {
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

impl<C, TickPairing> ConstraintSynthesizer<<C::TockGroup as PairingEngine>::Fr>
    for MiddleCircuit<C, TickPairing>
where
    C: CurvePair,
    TickPairing: PG<C::TickGroup>,
    <C::TickGroup as PairingEngine>::G1Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G2Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G1Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
    <C::TickGroup as PairingEngine>::G2Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
{
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<<C::TockGroup as PairingEngine>::Fr>,
    ) -> Result<(), SynthesisError> {
        let params = self.params;
        let proof = self.proof;
        let inputs = self.inputs;
        let input_gadgets;

        {
            let ns = r1cs_core::ns!(cs, "Allocate Input");
            let cs = ns.cs();
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
            let input_bytes = UInt8::new_input_vec(r1cs_core::ns!(cs, "Input"), &input_bytes[..])?;
            // 40 byte
            let element_size =
                <<<C::TickGroup as PairingEngine>::Fr as PrimeField>::Params as FftParameters>::BigInt::NUM_LIMBS * 8;
            input_gadgets = input_bytes
                .chunks(element_size)
                .map(|c| {
                    c.iter()
                        .flat_map(|b| b.to_bits_le().unwrap())
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

        let vk_var =
            InnerVkVar::<C, TickPairing>::new_witness(r1cs_core::ns!(cs, "Vk"), || Ok(&params.vk))?;
        let proof_var =
            InnerProofVar::<C, TickPairing>::new_witness(r1cs_core::ns!(cs, "Proof"), || {
                Ok(proof.clone())
            })?;
        <InnerVerifierGadget<C, TickPairing> as NIZKVerifierGadget<
            InnerProofSystem<C>,
            <C::TockGroup as PairingEngine>::Fr,
        >>::verify(&vk_var, input_gadgets.iter(), &proof_var)?
        .enforce_equal(&Boolean::TRUE)?;
        println!(
            "|---- Num constraints for sub-SNARK verification: {}",
            cs.num_constraints() - num_constraints
        );
        Ok(())
    }
}

pub struct OuterCircuit<C: CurvePair, TockPairing: PG<C::TockGroup>, TickPairing: PG<C::TickGroup>>
where
    <C::TickGroup as PairingEngine>::G1Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G2Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G1Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
    <C::TickGroup as PairingEngine>::G2Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
{
    inputs: Vec<<C::TickGroup as PairingEngine>::Fr>,
    params: Parameters<C::TockGroup>,
    proof: Proof<C::TockGroup>,
    _curve_pair: PhantomData<C>,
    _tock_pairing: PhantomData<TockPairing>,
    _tick_pairing: PhantomData<TickPairing>,
}

impl<C: CurvePair, TockPairing: PG<C::TockGroup>, TickPairing: PG<C::TickGroup>>
    OuterCircuit<C, TockPairing, TickPairing>
where
    <C::TickGroup as PairingEngine>::G1Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G2Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G1Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
    <C::TickGroup as PairingEngine>::G2Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
{
    pub fn new(
        inputs: Vec<<C::TickGroup as PairingEngine>::Fr>,
        params: Parameters<C::TockGroup>,
        proof: Proof<C::TockGroup>,
    ) -> Self {
        Self {
            inputs,
            params,
            proof,
            _curve_pair: PhantomData,
            _tock_pairing: PhantomData,
            _tick_pairing: PhantomData,
        }
    }
}

impl<C: CurvePair, TockPairing: PG<C::TockGroup>, TickPairing: PG<C::TickGroup>>
    ConstraintSynthesizer<<C::TickGroup as PairingEngine>::Fr>
    for OuterCircuit<C, TockPairing, TickPairing>
where
    <C::TickGroup as PairingEngine>::G1Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G2Projective: MulAssign<<C::TockGroup as PairingEngine>::Fq>,
    <C::TickGroup as PairingEngine>::G1Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
    <C::TickGroup as PairingEngine>::G2Affine:
        ToConstraintField<<<C::TockGroup as PairingEngine>::Fr as Field>::BasePrimeField>,
{
    fn generate_constraints(
        self,
        cs: ConstraintSystemRef<<C::TickGroup as PairingEngine>::Fr>,
    ) -> Result<(), SynthesisError> {
        let params = self.params;
        let proof = self.proof;
        let inputs = self.inputs;
        let mut input_gadgets = Vec::new();

        {
            let bigint_size =
                <<C::TickGroup as PairingEngine>::Fr as PrimeField>::BigInt::NUM_LIMBS * 64;
            let mut input_bits = Vec::new();
            let ns = r1cs_core::ns!(cs, "Allocate Input");
            let cs = ns.cs();

            for input in inputs.into_iter() {
                let input_gadget = FpVar::new_input(r1cs_core::ns!(cs, "Input"), || Ok(input))?;
                let mut fp_bits = input_gadget.to_bits_le()?;

                // Use 320 bits per element.
                for _ in fp_bits.len()..bigint_size {
                    fp_bits.push(Boolean::constant(false));
                }
                input_bits.extend_from_slice(&fp_bits);
            }

            // Pack input bits into field elements of the underlying circuit.
            let max_size = 8
                * (<<<C::TockGroup as PairingEngine>::Fr as PrimeField>::Params as FpParameters>::CAPACITY / 8)
                    as usize;
            let bigint_size =
                <<<C::TockGroup as PairingEngine>::Fr as PrimeField>::Params as FftParameters>::BigInt::NUM_LIMBS * 64;
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

        let vk_var =
            MiddleVkVar::<C, TockPairing>::new_witness(
                r1cs_core::ns!(cs, "Vk"),
                || Ok(&params.vk),
            )?;
        let proof_var =
            MiddleProofVar::<C, TockPairing>::new_witness(r1cs_core::ns!(cs, "Proof"), || {
                Ok(proof.clone())
            })?;
        <MiddleVerifierGadget<C, TockPairing> as NIZKVerifierGadget<
            MiddleProofSystem<C, TickPairing>,
            <C::TickGroup as PairingEngine>::Fr,
        >>::verify(&vk_var, input_gadgets.iter(), &proof_var)?
        .enforce_equal(&Boolean::TRUE)?;
        println!(
            "|---- Num constraints for sub-SNARK verification: {}",
            cs.num_constraints() - num_constraints
        );
        Ok(())
    }
}
