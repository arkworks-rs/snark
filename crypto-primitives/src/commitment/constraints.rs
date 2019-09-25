use crate::CommitmentScheme;
use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;
use std::fmt::Debug;

pub trait CommitmentGadget<C: CommitmentScheme, ConstraintF: Field> {
    type OutputGadget: EqGadget<ConstraintF>
        + ToBytesGadget<ConstraintF>
        + AllocGadget<C::Output, ConstraintF>
        + Clone
        + Sized
        + Debug;
    type ParametersGadget: AllocGadget<C::Parameters, ConstraintF> + Clone;
    type RandomnessGadget: AllocGadget<C::Randomness, ConstraintF> + Clone;

    fn check_commitment_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        input: &[UInt8],
        r: &Self::RandomnessGadget,
    ) -> Result<Self::OutputGadget, SynthesisError>;
}
