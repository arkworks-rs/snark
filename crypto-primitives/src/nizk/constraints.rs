use algebra_core::Field;
use core::borrow::Borrow;
use r1cs_core::{Namespace, SynthesisError};
use r1cs_std::prelude::*;

use crate::nizk::NIZK;

pub trait NIZKVerifierGadget<N: NIZK, ConstraintF: Field> {
    type PreparedVerificationKeyVar;
    type VerificationKeyVar: AllocVar<N::VerificationParameters, ConstraintF>
        + ToBytesGadget<ConstraintF>;
    type ProofVar: AllocVar<N::Proof, ConstraintF>;

    /// Optionally allocates `N::Proof` in `cs` without performing
    /// subgroup checks.
    ///
    /// The default implementation doesn't omit these checks.
    fn new_proof_unchecked<T: Borrow<N::Proof>>(
        cs: impl Into<Namespace<ConstraintF>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self::ProofVar, SynthesisError> {
        Self::ProofVar::new_variable(cs, f, mode)
    }

    /// Optionally allocates `N::VerificationParameters` in `cs`
    /// without performing subgroup checks.
    ///
    /// The default implementation doesn't omit these checks.
    fn new_verification_key_unchecked<T: Borrow<N::VerificationParameters>>(
        cs: impl Into<Namespace<ConstraintF>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self::VerificationKeyVar, SynthesisError> {
        Self::VerificationKeyVar::new_variable(cs, f, mode)
    }

    fn verify<'a, T: 'a + ToBitsGadget<ConstraintF> + ?Sized>(
        verification_key: &Self::VerificationKeyVar,
        input: impl IntoIterator<Item = &'a T>,
        proof: &Self::ProofVar,
    ) -> Result<(), SynthesisError> {
        Self::conditional_verify(verification_key, input, proof, &Boolean::constant(true))
    }

    fn conditional_verify<'a, T: 'a + ToBitsGadget<ConstraintF> + ?Sized>(
        verification_key: &Self::VerificationKeyVar,
        input: impl IntoIterator<Item = &'a T>,
        proof: &Self::ProofVar,
        condition: &Boolean<ConstraintF>,
    ) -> Result<(), SynthesisError>;

    fn conditional_verify_prepared<'a, T: 'a + ToBitsGadget<ConstraintF> + ?Sized>(
        prepared_verification_key: &Self::PreparedVerificationKeyVar,
        input: impl IntoIterator<Item = &'a T>,
        proof: &Self::ProofVar,
        condition: &Boolean<ConstraintF>,
    ) -> Result<(), SynthesisError>;
}
