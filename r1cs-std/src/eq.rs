use crate::{prelude::*, Vec};
use algebra::Field;
use r1cs_core::SynthesisError;

pub trait EqGadget<F: Field> {
    /// Output a `Boolean` value representing whether `self.value() == other.value()`.
    fn is_eq(&self, other: &Self) -> Result<Boolean<F>, SynthesisError>;

    /// Output a `Boolean` value representing whether `self.value() != other.value()`.
    fn is_neq(&self, other: &Self) -> Result<Boolean<F>, SynthesisError> {
        Ok(self.is_eq(other)?.not())
    }

    /// If `should_enforce == true`, enforce that `self` and `other` are equal; else,
    /// enforce a vacuously true statement.
    fn conditional_enforce_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        self.is_eq(&other)?
            .conditional_enforce_equal(&Boolean::constant(true), should_enforce)
    }

    /// Enforce that `self` and `other` are equal.
    fn enforce_equal(&self, other: &Self) -> Result<(), SynthesisError> {
        self.conditional_enforce_equal(other, &Boolean::constant(true))
    }

    /// If `should_enforce == true`, enforce that `self` and `other` are not equal; else,
    /// enforce a vacuously true statement.
    fn conditional_enforce_not_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        self.is_neq(&other)?
            .conditional_enforce_equal(&Boolean::constant(true), should_enforce)
    }

    /// Enforce that `self` and `other` are not equal.
    fn enforce_not_equal(&self, other: &Self) -> Result<(), SynthesisError> {
        self.conditional_enforce_not_equal(other, &Boolean::constant(true))
    }
}

impl<T: EqGadget<F> + R1CSVar<F>, F: Field> EqGadget<F> for [T] {
    fn is_eq(&self, other: &Self) -> Result<Boolean<F>, SynthesisError> {
        assert_eq!(self.len(), other.len());
        assert!(!self.is_empty());
        let mut results = Vec::with_capacity(self.len());
        for (a, b) in self.iter().zip(other) {
            results.push(a.is_eq(b)?);
        }
        Boolean::kary_and(&results)
    }

    fn conditional_enforce_equal(
        &self,
        other: &Self,
        condition: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        assert_eq!(self.len(), other.len());
        for (a, b) in self.iter().zip(other) {
            a.conditional_enforce_equal(b, condition)?;
        }
        Ok(())
    }

    fn conditional_enforce_not_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        assert_eq!(self.len(), other.len());
        let some_are_different = self.is_neq(other)?;
        if let Some(cs) = some_are_different.cs().or(should_enforce.cs()) {
            cs.enforce_constraint(
                some_are_different.lc(),
                should_enforce.lc(),
                should_enforce.lc(),
            )
        } else {
            // `some_are_different` and `should_enforce` are both constants
            assert!(some_are_different.value().unwrap());
            Ok(())
        }
    }
}

pub trait OrEqualsGadget<ConstraintF: Field>
where
    Self: Sized,
{
    /// If `should_enforce == true`, enforce that `self` equals
    /// (a) `first` (if `cond` is `true`)
    /// (b) `second` (if `cond` is `false`)
    fn conditional_enforce_equal_or(
        &self,
        cond: &Boolean<ConstraintF>,
        first: &Self,
        second: &Self,
        should_enforce: &Boolean<ConstraintF>,
    ) -> Result<(), SynthesisError>;

    fn enforce_equal_or(
        &self,
        cond: &Boolean<ConstraintF>,
        first: &Self,
        second: &Self,
    ) -> Result<(), SynthesisError> {
        self.conditional_enforce_equal_or(cond, first, second, &Boolean::Constant(true))
    }
}

impl<ConstraintF, T> OrEqualsGadget<ConstraintF> for T
where
    ConstraintF: Field,
    T: Sized + EqGadget<ConstraintF> + CondSelectGadget<ConstraintF>,
{
    fn conditional_enforce_equal_or(
        &self,
        cond: &Boolean<ConstraintF>,
        first: &Self,
        second: &Self,
        should_enforce: &Boolean<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let match_opt = cond.select(first, second)?;
        self.conditional_enforce_equal(&match_opt, should_enforce)
    }
}
