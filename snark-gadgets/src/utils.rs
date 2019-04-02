use crate::bits::{boolean::Boolean, uint8::UInt8};
use algebra::PairingEngine;
use snark::{ConstraintSystem, SynthesisError};
use std::borrow::Borrow;

/// If `condition == 1`, then enforces that `self` and `other` are equal;
/// otherwise, it doesn't enforce anything.
pub trait ConditionalEqGadget<E: PairingEngine>: Eq {
    fn conditional_enforce_equal<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError>;

    fn cost() -> usize;
}
impl<T: ConditionalEqGadget<E>, E: PairingEngine> ConditionalEqGadget<E> for [T] {
    fn conditional_enforce_equal<CS: ConstraintSystem<E>>(
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

pub trait EqGadget<E: PairingEngine>: Eq
where
    Self: ConditionalEqGadget<E>,
{
    fn enforce_equal<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        self.conditional_enforce_equal(cs, other, &Boolean::constant(true))
    }

    fn cost() -> usize {
        <Self as ConditionalEqGadget<E>>::cost()
    }
}

impl<T: EqGadget<E>, E: PairingEngine> EqGadget<E> for [T] {}

pub trait NEqGadget<E: PairingEngine>: Eq {
    fn enforce_not_equal<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError>;

    fn cost() -> usize;
}

pub trait ToBitsGadget<E: PairingEngine> {
    fn to_bits<CS: ConstraintSystem<E>>(&self, cs: CS) -> Result<Vec<Boolean>, SynthesisError>;

    /// Additionally checks if the produced list of booleans is 'valid'.
    fn to_bits_strict<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError>;
}

impl<E: PairingEngine> ToBitsGadget<E> for Boolean {
    fn to_bits<CS: ConstraintSystem<E>>(&self, _: CS) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(vec![self.clone()])
    }

    fn to_bits_strict<CS: ConstraintSystem<E>>(
        &self,
        _: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(vec![self.clone()])
    }
}

impl<E: PairingEngine> ToBitsGadget<E> for [Boolean] {
    fn to_bits<CS: ConstraintSystem<E>>(&self, _cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.to_vec())
    }

    fn to_bits_strict<CS: ConstraintSystem<E>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.to_vec())
    }
}
impl<E: PairingEngine> ToBitsGadget<E> for Vec<Boolean> {
    fn to_bits<CS: ConstraintSystem<E>>(&self, _cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.clone())
    }

    fn to_bits_strict<CS: ConstraintSystem<E>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.clone())
    }
}

impl<E: PairingEngine> ToBitsGadget<E> for [UInt8] {
    fn to_bits<CS: ConstraintSystem<E>>(&self, _cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        let mut result = Vec::with_capacity(&self.len() * 8);
        for byte in self {
            result.extend_from_slice(&byte.into_bits_le());
        }
        Ok(result)
    }

    fn to_bits_strict<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        self.to_bits(cs)
    }
}

pub trait ToBytesGadget<E: PairingEngine> {
    fn to_bytes<CS: ConstraintSystem<E>>(&self, cs: CS) -> Result<Vec<UInt8>, SynthesisError>;

    /// Additionally checks if the produced list of booleans is 'valid'.
    fn to_bytes_strict<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError>;
}

/// If condition is `true`, return `first`; else, select `second`.
pub trait CondSelectGadget<E: PairingEngine>
where
    Self: Sized,
{
    fn conditionally_select<CS: ConstraintSystem<E>>(
        cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError>;

    fn cost() -> usize;
}

/// Uses two bits to perform a lookup into a table
pub trait TwoBitLookupGadget<E: PairingEngine>
where
    Self: Sized,
{
    type TableConstant;
    fn two_bit_lookup<CS: ConstraintSystem<E>>(
        cs: CS,
        bits: &[Boolean],
        constants: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError>;

    fn cost() -> usize;
}

pub trait OrEqualsGadget<E: PairingEngine>
where
    Self: Sized,
{
    fn enforce_equal_or<CS: ConstraintSystem<E>>(
        cs: CS,
        cond: &Boolean,
        var: &Self,
        first: &Self,
        second: &Self,
    ) -> Result<(), SynthesisError>;

    fn cost() -> usize;
}

impl<E: PairingEngine, T: Sized + ConditionalOrEqualsGadget<E>> OrEqualsGadget<E> for T {
    fn enforce_equal_or<CS: ConstraintSystem<E>>(
        cs: CS,
        cond: &Boolean,
        var: &Self,
        first: &Self,
        second: &Self,
    ) -> Result<(), SynthesisError> {
        Self::conditional_enforce_equal_or(cs, cond, var, first, second, &Boolean::Constant(true))
    }

    fn cost() -> usize {
        <Self as ConditionalOrEqualsGadget<E>>::cost()
    }
}

pub trait ConditionalOrEqualsGadget<E: PairingEngine>
where
    Self: Sized,
{
    fn conditional_enforce_equal_or<CS: ConstraintSystem<E>>(
        cs: CS,
        cond: &Boolean,
        var: &Self,
        first: &Self,
        second: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError>;

    fn cost() -> usize;
}

impl<E: PairingEngine, T: Sized + ConditionalEqGadget<E> + CondSelectGadget<E>>
    ConditionalOrEqualsGadget<E> for T
{
    fn conditional_enforce_equal_or<CS: ConstraintSystem<E>>(
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
        <Self as ConditionalEqGadget<E>>::cost() + <Self as CondSelectGadget<E>>::cost()
    }
}

pub trait AllocGadget<V, E: PairingEngine>
where
    Self: Sized,
    V: ?Sized,
{
    fn alloc<F, T, CS: ConstraintSystem<E>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>;

    fn alloc_checked<F, T, CS: ConstraintSystem<E>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>,
    {
        Self::alloc(cs, f)
    }

    fn alloc_input<F, T, CS: ConstraintSystem<E>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>;

    fn alloc_input_checked<F, T, CS: ConstraintSystem<E>>(
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

impl<I, E: PairingEngine, A: AllocGadget<I, E>> AllocGadget<[I], E> for Vec<A> {
    fn alloc<F, T, CS: ConstraintSystem<E>>(mut cs: CS, f: F) -> Result<Self, SynthesisError>
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

    fn alloc_input<F, T, CS: ConstraintSystem<E>>(mut cs: CS, f: F) -> Result<Self, SynthesisError>
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

    fn alloc_checked<F, T, CS: ConstraintSystem<E>>(
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

    fn alloc_input_checked<F, T, CS: ConstraintSystem<E>>(
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
