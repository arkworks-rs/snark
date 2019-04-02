use algebra::PairingEngine;
use snark::{ConstraintSystem, SynthesisError};

use crate::crypto_primitives::nizk::NIZK;
use snark_gadgets::utils::{AllocGadget, ToBitsGadget, ToBytesGadget};

pub mod gm17;

pub trait NIZKVerifierGadget<N: NIZK, E: PairingEngine> {
    type VerificationKeyGadget: AllocGadget<N::VerificationParameters, E> + ToBytesGadget<E>;

    type ProofGadget: AllocGadget<N::Proof, E>;

    fn check_verify<'a, CS, I, T>(
        cs: CS,
        verification_key: &Self::VerificationKeyGadget,
        input: I,
        proof: &Self::ProofGadget,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<E>,
        I: Iterator<Item = &'a T>,
        T: 'a + ToBitsGadget<E> + ?Sized;
}
