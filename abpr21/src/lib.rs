//! An implementation of the [`BPR20`] zkSNARK.
//!
//! [`BPR20`]: https://eprint.iacr.org/2020/1306
#![cfg_attr(not(feature = "std"), no_std)]
#![warn(
    unused,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
    missing_docs
)]
#![allow(clippy::many_single_char_names, clippy::op_ref)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate ark_std;

#[cfg(feature = "r1cs")]
#[macro_use]
extern crate derivative;

/// Reduce an R1CS instance to a *Quadratic Arithmetic Program* instance.
pub(crate) mod r1cs_to_qap;

/// Data structures used by the prover, verifier, and generator.
pub mod data_structures;

/// Generate public parameters for the BPR20 zkSNARK construction.
pub mod generator;

/// Create proofs for the BPR20 zkSNARK construction.
pub mod prover;

/// Verify proofs for the BPR20 zkSNARK construction.
pub mod verifier;

/// Constraints for the BPR20 verifier.
#[cfg(feature = "r1cs")]
pub mod constraints;

#[cfg(test)]
mod test;

pub use self::data_structures::*;
pub use self::{generator::*, prover::*, verifier::*};

use ark_crypto_primitives::snark::*;
use ark_ec::PairingEngine;
use ark_relations::r1cs::{ConstraintSynthesizer, SynthesisError};
use ark_std::rand::RngCore;
use ark_std::{marker::PhantomData, vec::Vec};

/// The SNARK of [[BPR20]](https://eprint.iacr.org/2020/1306.pdf).
pub struct BPR20<E: PairingEngine> {
    e_phantom: PhantomData<E>,
}

impl<E: PairingEngine> SNARK<E::Fr> for BPR20<E> {
    type ProvingKey = ProvingKey<E>;
    type VerifyingKey = VerifyingKey<E>;
    type Proof = Proof<E>;
    type ProcessedVerifyingKey = PreparedVerifyingKey<E>;
    type Error = SynthesisError;

    fn circuit_specific_setup<C: ConstraintSynthesizer<E::Fr>, R: RngCore>(
        circuit: C,
        rng: &mut R,
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), Self::Error> {
        let pk = generate_random_parameters::<E, C, R>(circuit, rng)?;
        let vk = pk.vk.clone();

        Ok((pk, vk))
    }

    fn prove<C: ConstraintSynthesizer<E::Fr>, R: RngCore>(
        pk: &Self::ProvingKey,
        circuit: C,
        rng: &mut R,
    ) -> Result<Self::Proof, Self::Error> {
        create_random_proof::<E, _, _>(circuit, pk, rng)
    }

    fn process_vk(
        circuit_vk: &Self::VerifyingKey,
    ) -> Result<Self::ProcessedVerifyingKey, Self::Error> {
        Ok(prepare_verifying_key(circuit_vk))
    }

    fn verify_with_processed_vk(
        circuit_pvk: &Self::ProcessedVerifyingKey,
        x: &[E::Fr],
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error> {
        Ok(verify_proof(&circuit_pvk, proof, &x)?)
    }
}



impl<E: PairingEngine> CircuitSpecificSetupSNARK<E::Fr> for BPR20<E> {}
