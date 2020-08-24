use crate::commitment::CommitmentScheme;
use algebra_core::Field;
use core::fmt::Debug;
use r1cs_core::SynthesisError;
use r1cs_std::prelude::*;

pub trait CommitmentGadget<C: CommitmentScheme, ConstraintF: Field> {
    type OutputVar: EqGadget<ConstraintF>
        + ToBytesGadget<ConstraintF>
        + AllocVar<C::Output, ConstraintF>
        + R1CSVar<ConstraintF>
        + Clone
        + Sized
        + Debug;
    type ParametersVar: AllocVar<C::Parameters, ConstraintF> + Clone;
    type RandomnessVar: AllocVar<C::Randomness, ConstraintF> + Clone;

    fn commit(
        parameters: &Self::ParametersVar,
        input: &[UInt8<ConstraintF>],
        r: &Self::RandomnessVar,
    ) -> Result<Self::OutputVar, SynthesisError>;
}
