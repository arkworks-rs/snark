use algebra::Field;
use primitives::CommitmentScheme;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;
use std::fmt::Debug;

pub mod blake2s;
pub mod injective_map;
pub mod pedersen;

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
