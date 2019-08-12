use crate::prelude::*;
use r1cs_core::{ConstraintSystem, SynthesisError};
use algebra::Field;


/// If condition is `true`, return `first`; else, select `second`.
pub trait CondSelectGadget<ConstraintF: Field>
where
    Self: Sized,
{
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError>;

    fn cost() -> usize;
}

/// Uses two bits to perform a lookup into a table
pub trait TwoBitLookupGadget<ConstraintF: Field>
where
    Self: Sized,
{
    type TableConstant;
    fn two_bit_lookup<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        bits: &[Boolean],
        constants: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError>;

    fn cost() -> usize;
}


