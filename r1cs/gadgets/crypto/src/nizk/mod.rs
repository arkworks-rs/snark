use algebra::{bytes::ToBytes, Field};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

pub mod gm17;
pub mod groth16;

pub trait NIZK {
    type Circuit;
    type AssignedCircuit;
    type VerifierInput: ?Sized;
    type ProvingParameters: Clone;
    type VerificationParameters: Clone + Default;
    type PreparedVerificationParameters: Clone + Default + From<Self::VerificationParameters>;
    type Proof: ToBytes + Clone + Default;
}

pub trait NIZKVerifierGadget<N: NIZK, ConstraintF: Field> {
    type VerificationKeyGadget: AllocGadget<N::VerificationParameters, ConstraintF>
        + ToBytesGadget<ConstraintF>;

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
