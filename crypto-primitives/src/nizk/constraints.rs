use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use crate::nizk::NIZK;

pub trait NIZKVerifierGadget<N: NIZK, ConstraintF: Field> {
    type VerificationKeyGadget: AllocGadget<N::VerificationParameters, ConstraintF> + ToBytesGadget<ConstraintF>;

    type ProofGadget: AllocGadget<N::Proof, ConstraintF>;

    fn check_verify<'a, CS, I, T>(
        cs: CS,
        verification_key: &Self::VerificationKeyGadget,
        input: I,
        proof: &Self::ProofGadget,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<ConstraintF>,
        I: Iterator<Item = &'a T>,
        T: 'a + ToBitsGadget<ConstraintF> + ?Sized;
}
