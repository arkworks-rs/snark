use crate::Parameters;
use algebra_core::PairingEngine;
use ff_fft::GeneralEvaluationDomain;
use r1cs_core::{ConstraintSynthesizer, SynthesisError};
use rand::Rng;

pub mod generic;

/// Generates a random common reference string for
/// a circuit.
#[inline]
pub fn generate_random_parameters<E, C, R>(
    circuit: C,
    rng: &mut R,
) -> Result<Parameters<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    R: Rng,
{
    self::generic::generate_random_parameters::<E, C, GeneralEvaluationDomain<E::Fr>, R>(
        circuit, rng,
    )
}

/// Create parameters for a circuit, given some toxic waste.
#[inline]
pub fn generate_parameters<E, C, R>(
    circuit: C,
    alpha: E::Fr,
    beta: E::Fr,
    gamma: E::Fr,
    delta: E::Fr,
    rng: &mut R,
) -> Result<Parameters<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    R: Rng,
{
    self::generic::generate_parameters::<E, C, GeneralEvaluationDomain<E::Fr>, R>(
        circuit, alpha, beta, gamma, delta, rng,
    )
}
