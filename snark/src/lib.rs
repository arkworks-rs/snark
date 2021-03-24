//! This crate contains traits that define the basic behaviour of SNARKs.

#![cfg_attr(not(feature = "std"), no_std)]
#![warn(
    unused,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
    missing_docs
)]
#![forbid(unsafe_code)]

use ark_ff::{PrimeField, ToBytes};
use ark_relations::r1cs::ConstraintSynthesizer;
use ark_std::rand::{CryptoRng, RngCore};
use core::fmt::Debug;

/// The basic functionality for a SNARK.
pub trait SNARK<F: PrimeField> {
    /// The information required by the prover to produce a proof for a specific
    /// circuit *C*.
    type ProvingKey: Clone;

    /// The information required by the verifier to check a proof for a specific
    /// circuit *C*.
    type VerifyingKey: Clone + ToBytes;

    /// The proof output by the prover.
    type Proof: Clone;

    /// This contains the verification key, but preprocessed to enable faster
    /// verification.
    type ProcessedVerifyingKey: Clone;

    /// Errors encountered during setup, proving, or verification.
    type Error: 'static + ark_std::error::Error;

    /// Takes in a description of a computation (specified in R1CS constraints),
    /// and samples proving and verification keys for that circuit.
    fn circuit_specific_setup<C: ConstraintSynthesizer<F>, R: RngCore + CryptoRng>(
        circuit: C,
        rng: &mut R,
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), Self::Error>;

    /// Generates a proof of satisfaction of the arithmetic circuit C (specified
    /// as R1CS constraints).
    fn prove<C: ConstraintSynthesizer<F>, R: RngCore + CryptoRng>(
        circuit_pk: &Self::ProvingKey,
        circuit: C,
        rng: &mut R,
    ) -> Result<Self::Proof, Self::Error>;

    /// Checks that `proof` is a valid proof of the satisfaction of circuit
    /// encoded in `circuit_vk`, with respect to the public input `public_input`,
    /// specified as R1CS constraints.
    fn verify(
        circuit_vk: &Self::VerifyingKey,
        public_input: &[F],
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error> {
        let pvk = Self::process_vk(circuit_vk)?;
        Self::verify_with_processed_vk(&pvk, public_input, proof)
    }

    /// Preprocesses `circuit_vk` to enable faster verification.
    fn process_vk(
        circuit_vk: &Self::VerifyingKey,
    ) -> Result<Self::ProcessedVerifyingKey, Self::Error>;

    /// Checks that `proof` is a valid proof of the satisfaction of circuit
    /// encoded in `circuit_pvk`, with respect to the public input `public_input`,
    /// specified as R1CS constraints.
    fn verify_with_processed_vk(
        circuit_pvk: &Self::ProcessedVerifyingKey,
        public_input: &[F],
        proof: &Self::Proof,
    ) -> Result<bool, Self::Error>;
}

/// A SNARK with (only) circuit-specific setup.
pub trait CircuitSpecificSetupSNARK<F: PrimeField>: SNARK<F> {
    /// The setup algorithm for circuit-specific SNARKs. By default, this
    /// just invokes `<Self as SNARK<F>>::circuit_specific_setup(...)`.
    fn setup<C: ConstraintSynthesizer<F>, R: RngCore + CryptoRng>(
        circuit: C,
        rng: &mut R,
    ) -> Result<(Self::ProvingKey, Self::VerifyingKey), Self::Error> {
        <Self as SNARK<F>>::circuit_specific_setup(circuit, rng)
    }
}

/// A helper type for universal-setup SNARKs, which must infer their computation
/// size bounds.
pub enum UniversalSetupIndexError<Bound, E> {
    /// The provided universal public parameters were insufficient to encode
    /// the given circuit.
    NeedLargerBound(Bound),
    /// Other errors occurred during indexing.
    Other(E),
}

/// A SNARK with universal setup. That is, a SNARK where the trusted setup is
/// circuit-independent.
pub trait UniversalSetupSNARK<F: PrimeField>: SNARK<F> {
    /// Specifies how to bound the size of public parameters required to
    /// generate the index proving and verification keys for a given
    /// circuit.
    type ComputationBound: Clone + Default + Debug;
    /// Specifies the type of universal public parameters.
    type PublicParameters: Clone + Debug;

    /// Specifies how to bound the size of public parameters required to
    /// generate the index proving and verification keys for a given
    /// circuit.
    fn universal_setup<R: RngCore + CryptoRng>(
        compute_bound: &Self::ComputationBound,
        rng: &mut R,
    ) -> Result<Self::PublicParameters, Self::Error>;

    /// Indexes the public parameters according to the circuit `circuit`, and
    /// outputs circuit-specific proving and verification keys.
    fn index<C: ConstraintSynthesizer<F>, R: RngCore + CryptoRng>(
        pp: &Self::PublicParameters,
        circuit: C,
        rng: &mut R,
    ) -> Result<
        (Self::ProvingKey, Self::VerifyingKey),
        UniversalSetupIndexError<Self::ComputationBound, Self::Error>,
    >;
}
