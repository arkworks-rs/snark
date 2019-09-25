use algebra::Field;
use std::fmt::Debug;

use crate::crh::FixedLengthCRH;
use r1cs_core::{ConstraintSystem, SynthesisError};

use r1cs_std::prelude::*;

pub trait FixedLengthCRHGadget<H: FixedLengthCRH, ConstraintF: Field>: Sized {
    type OutputGadget: ConditionalEqGadget<ConstraintF>
        + EqGadget<ConstraintF>
        + ToBytesGadget<ConstraintF>
        + CondSelectGadget<ConstraintF>
        + AllocGadget<H::Output, ConstraintF>
        + Debug
        + Clone
        + Sized;
    type ParametersGadget: AllocGadget<H::Parameters, ConstraintF> + Clone;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        input: &[UInt8],
    ) -> Result<Self::OutputGadget, SynthesisError>;
}
