use crate::prelude::*;
use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};

/// Specifies how to generate constraints that check for equality for two variables of type `Self`.
pub trait EqGadget<ConstraintF: Field>: Eq {
    /// Output a `Boolean` value representing whether `self.value() == other.value()`.
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError>;

    /// Output a `Boolean` value representing whether `self.value() != other.value()`.
    ///
    /// By default, this is defined as `self.is_eq(other)?.not()`.
    fn is_neq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        Ok(self.is_eq(cs, other)?.not())
    }

    /// If `should_enforce == true`, enforce that `self` and `other` are equal; else,
    /// enforce a vacuously true statement.
    ///
    /// A safe default implementation is provided that generates the following constraints:
    /// `self.is_eq(other)?.conditional_enforce_equal(&Boolean::TRUE, should_enforce)`.
    ///
    /// More efficient specialized implementation may be possible; implementors
    /// are encouraged to carefully analyze the efficiency and safety of these.
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.is_eq(cs.ns(|| "is_eq(self, other)"), &other)?
            .conditional_enforce_equal(
                cs.ns(|| "enforce condition"),
                &Boolean::constant(true),
                should_enforce,
            )
    }

    /// Enforce that `self` and `other` are equal.
    ///
    /// A safe default implementation is provided that generates the following constraints:
    /// `self.conditional_enforce_equal(other, &Boolean::TRUE)`.
    ///
    /// More efficient specialized implementation may be possible; implementors
    /// are encouraged to carefully analyze the efficiency and safety of these.
    fn enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        self.conditional_enforce_equal(cs, other, &Boolean::constant(true))
    }

    /// If `should_enforce == true`, enforce that `self` and `other` are *not* equal; else,
    /// enforce a vacuously true statement.
    ///
    /// A safe default implementation is provided that generates the following constraints:
    /// `self.is_neq(other)?.conditional_enforce_equal(&Boolean::TRUE, should_enforce)`.
    ///
    /// More efficient specialized implementation may be possible; implementors
    /// are encouraged to carefully analyze the efficiency and safety of these.
    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.is_neq(cs.ns(|| "is_neq(self, other)"), &other)?
            .conditional_enforce_equal(
                cs.ns(|| "enforce condition"),
                &Boolean::constant(true),
                should_enforce,
            )
    }

    /// Enforce that `self` and `other` are *not* equal.
    ///
    /// A safe default implementation is provided that generates the following constraints:
    /// `self.conditional_enforce_not_equal(other, &Boolean::TRUE)`.
    ///
    /// More efficient specialized implementation may be possible; implementors
    /// are encouraged to carefully analyze the efficiency and safety of these.
    fn enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        self.conditional_enforce_not_equal(cs, other, &Boolean::constant(true))
    }
}

impl<T: EqGadget<ConstraintF>, ConstraintF: Field> EqGadget<ConstraintF> for [T] {
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        assert_eq!(self.len(), other.len());
        assert!(!self.is_empty());
        let mut results = Vec::with_capacity(self.len());
        for (i, (a, b)) in self.iter().zip(other).enumerate() {
            results.push(a.is_eq(cs.ns(|| format!("is_eq_{}", i)), b)?);
        }
        Boolean::kary_and(cs.ns(|| "kary and"), &results)
    }

    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        assert_eq!(self.len(), other.len());
        for (i, (a, b)) in self.iter().zip(other).enumerate() {
            a.conditional_enforce_equal(
                cs.ns(|| format!("conditional_enforce_equal_{}", i)),
                b,
                condition,
            )?;
        }
        Ok(())
    }

    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        assert_eq!(self.len(), other.len());
        let some_are_different = self.is_neq(cs.ns(|| "is_neq"), other)?;
        if some_are_different.get_value().is_some() && should_enforce.get_value().is_some() {
            assert!(some_are_different.get_value().unwrap());
            Ok(())
        } else {
            some_are_different.conditional_enforce_equal(
                cs.ns(|| "conditional_enforce_equal"),
                should_enforce,
                should_enforce,
            )?;
            Ok(())
        }
    }
}
