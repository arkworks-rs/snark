use crate::prelude::*;
use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};

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

    fn two_bit_lookup_lc<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        precomp: &Boolean,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError>;

    fn cost() -> usize;

    fn cost_of_lc() -> usize {
        0
    }
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
