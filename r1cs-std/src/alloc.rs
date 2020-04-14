use crate::Vec;
use algebra::Field;
use core::borrow::Borrow;
use r1cs_core::{ConstraintSystem, SynthesisError};

pub trait AllocGadget<V, ConstraintF: Field>
where
    Self: Sized,
    V: ?Sized,
{
    fn alloc_constant<T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<V>;

    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>;

    fn alloc_checked<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<V>,
    {
        Self::alloc(cs, f)
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
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

impl<I, ConstraintF: Field, A: AllocGadget<I, ConstraintF>> AllocGadget<[I], ConstraintF>
    for Vec<A>
{
    #[inline]
    fn alloc_constant<T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<[I]>,
    {
        let mut vec = Vec::new();
        for (i, value) in t.borrow().iter().enumerate() {
            vec.push(A::alloc_constant(cs.ns(|| format!("value_{}", i)), value)?);
        }
        Ok(vec)
    }

    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
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

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
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
