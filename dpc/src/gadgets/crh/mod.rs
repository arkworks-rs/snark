use algebra::Field;
use std::fmt::Debug;

use crate::crypto_primitives::crh::FixedLengthCRH;
use snark::{ConstraintSystem, SynthesisError};

use snark_gadgets::{
    uint8::UInt8,
    utils::{AllocGadget, CondSelectGadget, ConditionalEqGadget, EqGadget, ToBytesGadget},
};

pub mod injective_map;
pub mod pedersen;

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

    fn cost() -> usize;
}
