use super::ConstraintSystemRef;
use ark_ff::Field;

/// A namespaced `ConstraintSystemRef`.
#[derive(Clone)]
pub struct Namespace<F: Field> {
    inner: ConstraintSystemRef<F>,
    id: Option<tracing::Id>,
}

impl<F: Field> From<ConstraintSystemRef<F>> for Namespace<F> {
    fn from(other: ConstraintSystemRef<F>) -> Self {
        Self {
            inner: other,
            id: None,
        }
    }
}

impl<F: Field> Namespace<F> {
    /// Construct a new `Namespace`.
    pub fn new(inner: ConstraintSystemRef<F>, id: Option<tracing::Id>) -> Self {
        Self { inner, id }
    }

    /// Obtain the inner `ConstraintSystemRef<F>`.
    pub fn cs(&self) -> ConstraintSystemRef<F> {
        self.inner.clone()
    }

    /// Manually leave the namespace.
    pub fn leave_namespace(self) {
        drop(self)
    }
}

impl<F: Field> Drop for Namespace<F> {
    fn drop(&mut self) {
        if let Some(id) = self.id.as_ref() {
            tracing::dispatcher::get_default(|dispatch| dispatch.exit(id))
        }
        let _ = self.inner;
    }
}

/// Creates namespaces for different parts of a circuit when generating constraints. 
/// Here, a namespace is equivalent to having a unique span for each part of the circuit. 
/// For more information on spans, see the [tracing](https://docs.rs/tracing) crate.
///
/// Takes in a reference to a Constraint System and a string slice representing
/// the name of the namespace. The name is used to identify the namespace.
///
/// # Simple Example of using namespaces
///
/// ```rust,ignore
/// use ark_ff::Field;
/// use ark_r1cs_std::prelude::*;
/// use ark_relations::gr1cs::{ConstraintSystemRef, SynthesisError};
///
/// // Define the circuit structure
/// pub struct SimpleAdditionCircuit<F: Field> {
///     pub a: F,  // Public input 1
///     pub b: F,  // Public input 2
///     pub c: F,  // Witness (private input)
/// }
///
/// impl<F: Field> SimpleAdditionCircuit<F> {
///     pub fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
///         // Create a namespace for public inputs
///         let cs = ns!(cs, "public_inputs");
///         let a_var = FpVar::new_input(cs.clone(), || Ok(self.a))?;
///         let b_var = FpVar::new_input(cs.clone(), || Ok(self.b))?;
///
///         // Create another namespace for the witness
///         let cs = ns!(cs, "witness");
///         let c_var = FpVar::new_witness(cs.clone(), || Ok(self.c))?;
///
///         // Create another namespace for the addition constraint
///         let cs = ns!(cs, "addition_constraint");
///         let sum = &a_var + &b_var;
///
///         // Enforce that a + b = c
///         sum.enforce_equal(&c_var)?;
///
///         Ok(())
///     }
/// }
/// ```
#[macro_export]
macro_rules! ns {
    ($cs:expr, $name:expr) => {{
        // Define a span with `gr1cs` as the target and the given name
        let span = $crate::gr1cs::info_span!(target: "gr1cs", $name);
        let id = span.id();
        // Enter the span
        let _enter_guard = span.enter();
        // We want the span and the guard to live forever so we forget them
        core::mem::forget(_enter_guard);
        core::mem::forget(span);
        // Create a new namespace
        $crate::gr1cs::Namespace::new($cs.clone(), id)
    }};
}
