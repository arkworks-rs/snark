//! Core interface for working with Rank-1 Constraint Systems (R1CS).

#![cfg_attr(not(feature = "std"), no_std)]
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, private_in_public, variant_size_differences)]
#![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
#![deny(unused_attributes, unused_imports, unused_mut, missing_docs)]
#![deny(renamed_and_removed_lints, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use, const_err)]
#![deny(unsafe_code)]

#[cfg(not(feature = "std"))]
pub extern crate alloc;

#[cfg(not(feature = "std"))]
pub use alloc::{
    collections::{BTreeMap, BTreeSet},
    format,
    rc::Rc,
    string::{String, ToString},
    vec,
    vec::Vec,
};

#[cfg(feature = "std")]
pub use std::{
    collections::{BTreeMap, BTreeSet},
    format,
    rc::Rc,
    string::{String, ToString},
    vec,
    vec::Vec,
};

mod constraint_system;
mod error;
mod impl_lc;
#[cfg(feature = "std")]
mod trace;

#[cfg(feature = "std")]
pub use crate::trace::{ConstraintLayer, ConstraintTrace, TraceStep, TracingMode};
#[cfg(feature = "std")]
pub use tracing::info_span;

pub use algebra_core::{Field, ToConstraintField};
pub use constraint_system::{
    ConstraintMatrices, ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, Namespace,
    SynthesisMode,
};
pub use error::SynthesisError;

use core::cmp::Ordering;

/// A linear combination of variables according to associated coefficients.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinearCombination<F: Field>(pub Vec<(F, Variable)>);

/// A sparse representation of constraint matrices.
pub type Matrix<F> = Vec<Vec<(F, usize)>>;

/// Represents the different kinds of variables present in a constraint system.
#[derive(Copy, Clone, PartialEq, Debug, Eq)]
pub enum Variable {
    /// Represents the "zero" constant.
    Zero,
    /// Represents of the "one" constant.
    One,
    /// Represents a public instance variable.
    Instance(usize),
    /// Represents a private witness variable.
    Witness(usize),
    /// Represents of a linear combination.
    SymbolicLc(LcIndex),
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
/// An opaque counter for symbolic linear combinations.
pub struct LcIndex(usize);

/// Generate a `LinearCombination` from arithmetic expressions involving
/// `Variable`s.
#[macro_export]
macro_rules! lc {
    () => {
        $crate::LinearCombination::zero()
    };
}

/// Generate a `Namespace` with name `name` from `ConstraintSystem` `cs`.
/// `name` must be a `&'static str`.
#[macro_export]
macro_rules! ns {
    ($cs:expr, $name:expr) => {{
        #[cfg(feature = "std")]
        {
            let span = $crate::info_span!(target: "r1cs", $name);
            let id = span.id();
            let _enter_guard = span.enter();
            core::mem::forget(_enter_guard);
            core::mem::forget(span);
            $crate::Namespace::new($cs.clone(), id)
        }
        #[cfg(not(feature = "std"))]
        {
            $crate::Namespace::from($cs.clone())
        }
    }};
}

impl Variable {
    /// Is `self` the zero variable?
    #[inline]
    pub fn is_zero(&self) -> bool {
        match self {
            Variable::Zero => true,
            _ => false,
        }
    }

    /// Is `self` the one variable?
    #[inline]
    pub fn is_one(&self) -> bool {
        match self {
            Variable::One => true,
            _ => false,
        }
    }

    /// Is `self` an instance variable?
    #[inline]
    pub fn is_instance(&self) -> bool {
        match self {
            Variable::Instance(_) => true,
            _ => false,
        }
    }

    /// Is `self` a witness variable?
    #[inline]
    pub fn is_witness(&self) -> bool {
        match self {
            Variable::Witness(_) => true,
            _ => false,
        }
    }

    /// Is `self` a linear combination?
    #[inline]
    pub fn is_lc(&self) -> bool {
        match self {
            Variable::SymbolicLc(_) => true,
            _ => false,
        }
    }

    /// Get the `LcIndex` in `self` if `self.is_lc()`.
    #[inline]
    pub fn get_lc_index(&self) -> Option<LcIndex> {
        match self {
            Variable::SymbolicLc(index) => Some(*index),
            _ => None,
        }
    }

    /// Returns `Some(usize)` if `!self.is_lc()`, and `None` otherwise.
    #[inline]
    pub fn get_index_unchecked(&self, witness_offset: usize) -> Option<usize> {
        match self {
            // The one variable always has index 0
            Variable::One => Some(0),
            Variable::Instance(i) => Some(*i),
            Variable::Witness(i) => Some(witness_offset + *i),
            _ => None,
        }
    }
}

impl PartialOrd for Variable {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        use crate::Variable::*;
        match (self, other) {
            (Zero, Zero) => Some(Ordering::Equal),
            (One, One) => Some(Ordering::Equal),
            (Zero, _) => Some(Ordering::Less),
            (One, _) => Some(Ordering::Less),
            (_, Zero) => Some(Ordering::Greater),
            (_, One) => Some(Ordering::Greater),

            (Instance(i), Instance(j)) | (Witness(i), Witness(j)) => i.partial_cmp(j),
            (Instance(_), Witness(_)) => Some(Ordering::Less),
            (Witness(_), Instance(_)) => Some(Ordering::Greater),

            (SymbolicLc(i), SymbolicLc(j)) => i.partial_cmp(j),
            (_, SymbolicLc(_)) => Some(Ordering::Less),
            (SymbolicLc(_), _) => Some(Ordering::Greater),
        }
    }
}

impl Ord for Variable {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}
