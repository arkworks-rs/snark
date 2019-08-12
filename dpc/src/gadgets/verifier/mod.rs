use algebra::Field;
use snark::{ConstraintSystem, SynthesisError};

use crate::crypto_primitives::nizk::NIZK;
use snark_gadgets::utils::{AllocGadget, ToBitsGadget, ToBytesGadget};

pub mod gm17;

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
