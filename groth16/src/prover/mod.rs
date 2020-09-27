use crate::{Parameters, Proof};
use algebra_core::PairingEngine;
use ff_fft::GeneralEvaluationDomain;
use r1cs_core::{ConstraintSynthesizer, SynthesisError};
use rand::Rng;

pub mod generic;

#[inline]
pub fn create_random_proof<E, C, R>(
    circuit: C,
    params: &Parameters<E>,
    rng: &mut R,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    R: Rng,
{
    self::generic::create_random_proof::<E, C, GeneralEvaluationDomain<E::Fr>, R>(
        circuit, params, rng,
    )
}

#[inline]
pub fn create_proof_no_zk<E, C>(
    circuit: C,
    params: &Parameters<E>,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
{
    self::generic::create_proof_no_zk::<E, C, GeneralEvaluationDomain<E::Fr>>(circuit, params)
}

#[inline]
pub fn create_proof<E, C>(
    circuit: C,
    params: &Parameters<E>,
    r: E::Fr,
    s: E::Fr,
) -> Result<Proof<E>, SynthesisError>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
{
    self::generic::create_proof::<E, C, GeneralEvaluationDomain<E::Fr>>(circuit, params, r, s)
}
