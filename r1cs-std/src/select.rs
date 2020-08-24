use crate::prelude::*;
use algebra::Field;
use r1cs_core::SynthesisError;

/// If condition is `true`, return `true_value`; else, select `false_value`.
pub trait CondSelectGadget<ConstraintF: Field>
where
    Self: Sized,
{
    fn conditionally_select(
        cond: &Boolean<ConstraintF>,
        true_value: &Self,
        false_value: &Self,
    ) -> Result<Self, SynthesisError>;
}

/// Uses two bits to perform a lookup into a table
pub trait TwoBitLookupGadget<ConstraintF: Field>
where
    Self: Sized,
{
    type TableConstant;
    fn two_bit_lookup(
        bits: &[Boolean<ConstraintF>],
        constants: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError>;
}

/// Uses three bits to perform a lookup into a table, where the last bit
/// performs negation
pub trait ThreeBitCondNegLookupGadget<ConstraintF: Field>
where
    Self: Sized,
{
    type TableConstant;
    fn three_bit_cond_neg_lookup(
        bits: &[Boolean<ConstraintF>],
        b0b1: &Boolean<ConstraintF>,
        constants: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError>;
}
