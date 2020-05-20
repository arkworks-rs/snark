use crate::Error;
use algebra_core::PairingEngine;
use gm17::{
    create_random_proof, generate_random_parameters, prepare_verifying_key, verify_proof,
    Parameters, PreparedVerifyingKey, Proof, VerifyingKey,
};
use r1cs_core::ConstraintSynthesizer;
use rand::Rng;

use algebra_core::ToConstraintField;
use core::marker::PhantomData;

use super::NIZK;

#[cfg(feature = "r1cs")]
pub mod constraints;

/// Note: V should serialize its contents to `Vec<E::Fr>` in the same order as
/// during the constraint generation.
pub struct Gm17<
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    V: ToConstraintField<E::Fr> + ?Sized,
> {
    #[doc(hidden)]
    _engine: PhantomData<E>,
    #[doc(hidden)]
    _circuit: PhantomData<C>,
    #[doc(hidden)]
    _verifier_input: PhantomData<V>,
}

impl<E: PairingEngine, C: ConstraintSynthesizer<E::Fr>, V: ToConstraintField<E::Fr> + ?Sized> NIZK
    for Gm17<E, C, V>
{
    type Circuit = C;
    type AssignedCircuit = C;
    type VerifierInput = V;
    type ProvingParameters = Parameters<E>;
    type VerificationParameters = VerifyingKey<E>;
    type PreparedVerificationParameters = PreparedVerifyingKey<E>;
    type Proof = Proof<E>;

    fn setup<R: Rng>(
        circuit: Self::Circuit,
        rng: &mut R,
    ) -> Result<
        (
            Self::ProvingParameters,
            Self::PreparedVerificationParameters,
        ),
        Error,
    > {
        let nizk_time = start_timer!(|| "{Groth-Maller 2017}::Setup");
        let pp = generate_random_parameters::<E, Self::Circuit, R>(circuit, rng)?;
        let vk = prepare_verifying_key(&pp.vk);
        end_timer!(nizk_time);
        Ok((pp, vk))
    }

    fn prove<R: Rng>(
        pp: &Self::ProvingParameters,
        input_and_witness: Self::AssignedCircuit,
        rng: &mut R,
    ) -> Result<Self::Proof, Error> {
        let proof_time = start_timer!(|| "{Groth-Maller 2017}::Prove");
        let result = create_random_proof::<E, _, _>(input_and_witness, pp, rng)?;
        end_timer!(proof_time);
        Ok(result)
    }

    fn verify(
        vk: &Self::PreparedVerificationParameters,
        input: &Self::VerifierInput,
        proof: &Self::Proof,
    ) -> Result<bool, Error> {
        let verify_time = start_timer!(|| "{Groth-Maller 2017}::Verify");
        let conversion_time = start_timer!(|| "Convert input to E::Fr");
        let input = input.to_field_elements()?;
        end_timer!(conversion_time);
        let verification = start_timer!(|| format!("Verify proof w/ input len: {}", input.len()));
        let result = verify_proof(&vk, proof, &input)?;
        end_timer!(verification);
        end_timer!(verify_time);
        Ok(result)
    }
}
