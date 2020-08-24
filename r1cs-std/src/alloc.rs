use crate::Vec;
use algebra::Field;
use core::borrow::Borrow;
use r1cs_core::{Namespace, SynthesisError};

#[derive(Eq, PartialEq, Ord, PartialOrd, Debug, Copy, Clone)]
pub enum AllocationMode {
    Constant = 0,
    Input = 1,
    Witness = 2,
}

impl AllocationMode {
    // Outputs the maximum according to the relation `Constant < Input < Witness`.
    pub fn max(&self, other: Self) -> Self {
        use AllocationMode::*;
        match (self, other) {
            (Constant, _) => other,
            (Input, Constant) => *self,
            (Input, _) => other,
            (Witness, _) => *self,
        }
    }
}

pub trait AllocVar<V, F: Field>
where
    Self: Sized,
    V: ?Sized,
{
    fn new_variable<T: Borrow<V>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError>;

    fn new_constant(
        cs: impl Into<Namespace<F>>,
        t: impl Borrow<V>,
    ) -> Result<Self, SynthesisError> {
        Self::new_variable(cs, || Ok(t), AllocationMode::Constant)
    }

    fn new_input<T: Borrow<V>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
    ) -> Result<Self, SynthesisError> {
        Self::new_variable(cs, f, AllocationMode::Input)
    }

    fn new_witness<T: Borrow<V>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
    ) -> Result<Self, SynthesisError> {
        Self::new_variable(cs, f, AllocationMode::Witness)
    }
}

impl<I, F: Field, A: AllocVar<I, F>> AllocVar<[I], F> for Vec<A> {
    fn new_variable<T: Borrow<[I]>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        let mut vec = Vec::new();
        for value in f()?.borrow().iter() {
            vec.push(A::new_variable(cs.clone(), || Ok(value), mode)?);
        }
        Ok(vec)
    }
}
