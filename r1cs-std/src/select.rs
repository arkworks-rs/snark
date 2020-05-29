use crate::prelude::*;
use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};

/// If condition is `true`, return `true_value`; else, select `false_value`.
pub trait CondSelectGadget<ConstraintF: Field>
where
    Self: Sized,
{
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        cond: &Boolean,
        true_value: &Self,
        false_value: &Self,
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

/// Uses three bits to perform a lookup into a table, where the last bit
/// performs negation
pub trait ThreeBitCondNegLookupGadget<ConstraintF: Field>
where
    Self: Sized,
{
    type TableConstant;
    fn three_bit_cond_neg_lookup<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        bits: &[Boolean],
        b0b1: &Boolean,
        constants: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError>;

    fn cost() -> usize;
}
