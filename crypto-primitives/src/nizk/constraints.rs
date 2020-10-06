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
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
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
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
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
    ) -> Result<Boolean<ConstraintF>, SynthesisError>;

    fn verify_prepared<'a, T: 'a + ToBitsGadget<ConstraintF> + ?Sized>(
        prepared_verification_key: &Self::PreparedVerificationKeyVar,
        input: impl IntoIterator<Item = &'a T>,
        proof: &Self::ProofVar,
    ) -> Result<Boolean<ConstraintF>, SynthesisError>;
}
