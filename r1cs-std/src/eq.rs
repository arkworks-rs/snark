use crate::prelude::*;
use r1cs_core::{ConstraintSystem, SynthesisError};
use algebra::Field;

/// If `condition == 1`, then enforces that `self` and `other` are equal;
/// otherwise, it doesn't enforce anything.
pub trait ConditionalEqGadget<ConstraintF: Field>: Eq {
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError>;

    fn cost() -> usize;
}
impl<T: ConditionalEqGadget<ConstraintF>, ConstraintF: Field> ConditionalEqGadget<ConstraintF> for [T] {
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        for (i, (a, b)) in self.iter().zip(other.iter()).enumerate() {
            let mut cs = cs.ns(|| format!("Iteration {}", i));
            a.conditional_enforce_equal(&mut cs, b, condition)?;
        }
        Ok(())
    }

    fn cost() -> usize {
        unimplemented!()
    }
}

pub trait EqGadget<ConstraintF: Field>: Eq
where
    Self: ConditionalEqGadget<ConstraintF>,
{
    fn enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        self.conditional_enforce_equal(cs, other, &Boolean::constant(true))
    }

    fn cost() -> usize {
        <Self as ConditionalEqGadget<ConstraintF>>::cost()
    }
}

impl<T: EqGadget<ConstraintF>, ConstraintF: Field> EqGadget<ConstraintF> for [T] {}

pub trait NEqGadget<ConstraintF: Field>: Eq {
    fn enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError>;

    fn cost() -> usize;
}

pub trait OrEqualsGadget<ConstraintF: Field>
where
    Self: Sized,
{
    fn enforce_equal_or<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        cond: &Boolean,
        var: &Self,
        first: &Self,
        second: &Self,
    ) -> Result<(), SynthesisError>;

    fn cost() -> usize;
}

impl<ConstraintF: Field, T: Sized + ConditionalOrEqualsGadget<ConstraintF>> OrEqualsGadget<ConstraintF> for T {
    fn enforce_equal_or<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        cond: &Boolean,
        var: &Self,
        first: &Self,
        second: &Self,
    ) -> Result<(), SynthesisError> {
        Self::conditional_enforce_equal_or(cs, cond, var, first, second, &Boolean::Constant(true))
    }

    fn cost() -> usize {
        <Self as ConditionalOrEqualsGadget<ConstraintF>>::cost()
    }
}

pub trait ConditionalOrEqualsGadget<ConstraintF: Field>
where
    Self: Sized,
{
    fn conditional_enforce_equal_or<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        cond: &Boolean,
        var: &Self,
        first: &Self,
        second: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError>;

    fn cost() -> usize;
}

impl<ConstraintF: Field, T: Sized + ConditionalEqGadget<ConstraintF> + CondSelectGadget<ConstraintF>>
    ConditionalOrEqualsGadget<ConstraintF> for T
{
    fn conditional_enforce_equal_or<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        var: &Self,
        first: &Self,
        second: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        let match_opt = Self::conditionally_select(
            &mut cs.ns(|| "conditional_select_in_or"),
            cond,
            first,
            second,
        )?;
        var.conditional_enforce_equal(&mut cs.ns(|| "equals_in_or"), &match_opt, should_enforce)
    }

    fn cost() -> usize {
        <Self as ConditionalEqGadget<ConstraintF>>::cost() + <Self as CondSelectGadget<ConstraintF>>::cost()
    }
}


