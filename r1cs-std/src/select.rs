use crate::prelude::*;
use algebra::Field;
use r1cs_core::SynthesisError;

/// Generates constraints for selecting between one of two values.
pub trait CondSelectGadget<ConstraintF: Field>
where
    Self: Sized,
{
    /// If `cond == &Boolean::TRUE`, then this returns `true_value`; else, returns `false_value`.
    ///
    /// # Note
    /// `Self::conditionally_select(cond, true_value, false_value)?` can be more succinctly written as
    /// `cond.select(&true_value, &false_value)?`.
    fn conditionally_select(
        cond: &Boolean<ConstraintF>,
        true_value: &Self,
        false_value: &Self,
    ) -> Result<Self, SynthesisError>;
}

/// Performs a lookup in a 4-element table using two bits.
pub trait TwoBitLookupGadget<ConstraintF: Field>
where
    Self: Sized,
{
    /// The type of values being looked up.
    type TableConstant;

    /// Interprets the slice `bits` as a two-bit integer `b = bits[0] + (bits[1] << 1)`,
    /// and then outputs `constants[b]`.
    ///
    /// For example, if `bits == [0, 1]`, and `constants == [0, 1, 2, 3]`, this method
    /// should output a variable corresponding to `2`.
    ///
    /// # Panics
    ///
    /// This method panics if `bits.len() != 2` or `constants.len() != 4`.
    fn two_bit_lookup(
        bits: &[Boolean<ConstraintF>],
        constants: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError>;
}

/// Uses three bits to perform a lookup into a table, where the last bit
/// conditionally negates the looked-up value.
pub trait ThreeBitCondNegLookupGadget<ConstraintF: Field>
where
    Self: Sized,
{
    /// The type of values being looked up.
    type TableConstant;

    /// Interprets the slice `bits` as a two-bit integer `b = bits[0] + (bits[1] << 1)`,
    /// and then outputs `constants[b] * c`, where `c = if bits[2] { -1 } else { 1 };`.
    ///
    /// That is, `bits[2]` conditionally negates the looked-up value.
    ///
    /// For example, if `bits == [1, 0, 1]`, and `constants == [0, 1, 2, 3]`, this method
    /// should output a variable corresponding to `-1`.
    ///
    /// # Panics
    ///
    /// This method panics if `bits.len() != 3` or `constants.len() != 4`.
    fn three_bit_cond_neg_lookup(
        bits: &[Boolean<ConstraintF>],
        b0b1: &Boolean<ConstraintF>,
        constants: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError>;
}
