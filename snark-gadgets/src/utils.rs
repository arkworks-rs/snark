use crate::bits::{boolean::Boolean, uint8::UInt8};
use algebra::Field;
use snark::{ConstraintSystem, SynthesisError};
use std::borrow::Borrow;

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

pub trait ToBitsGadget<ConstraintF: Field> {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Vec<Boolean>, SynthesisError>;

    /// Additionally checks if the produced list of booleans is 'valid'.
    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError>;
}

impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for Boolean {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, _: CS) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(vec![self.clone()])
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(vec![self.clone()])
    }
}

impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for [Boolean] {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, _cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.to_vec())
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.to_vec())
    }
}
impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for Vec<Boolean> {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, _cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.clone())
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.clone())
    }
}

impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for [UInt8] {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, _cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        let mut result = Vec::with_capacity(&self.len() * 8);
        for byte in self {
            result.extend_from_slice(&byte.into_bits_le());
        }
        Ok(result)
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        self.to_bits(cs)
    }
}

pub trait ToBytesGadget<ConstraintF: Field> {
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Vec<UInt8>, SynthesisError>;

    /// Additionally checks if the produced list of booleans is 'valid'.
    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError>;
}

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

pub trait AllocGadget<V, ConstraintF: Field>
where
    Self: Sized,
    V: ?Sized,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>;

    fn alloc_checked<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>,
    {
        Self::alloc(cs, f)
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>;

    fn alloc_input_checked<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>,
    {
        Self::alloc_input(cs, f)
    }
}

impl<I, ConstraintF: Field, A: AllocGadget<I, ConstraintF>> AllocGadget<[I], ConstraintF> for Vec<A> {
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(mut cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[I]>,
    {
        let mut vec = Vec::new();
        for (i, value) in f()?.borrow().iter().enumerate() {
            vec.push(A::alloc(&mut cs.ns(|| format!("value_{}", i)), || {
                Ok(value)
            })?);
        }
        Ok(vec)
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(mut cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[I]>,
    {
        let mut vec = Vec::new();
        for (i, value) in f()?.borrow().iter().enumerate() {
            vec.push(A::alloc_input(
                &mut cs.ns(|| format!("value_{}", i)),
                || Ok(value),
            )?);
        }
        Ok(vec)
    }

    fn alloc_checked<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[I]>,
    {
        let mut vec = Vec::new();
        for (i, value) in f()?.borrow().iter().enumerate() {
            vec.push(A::alloc_checked(
                &mut cs.ns(|| format!("value_{}", i)),
                || Ok(value),
            )?);
        }
        Ok(vec)
    }

    fn alloc_input_checked<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<[I]>,
    {
        let mut vec = Vec::new();
        for (i, value) in f()?.borrow().iter().enumerate() {
            vec.push(A::alloc_input_checked(
                &mut cs.ns(|| format!("value_{}", i)),
                || Ok(value),
            )?);
        }
        Ok(vec)
    }
}
