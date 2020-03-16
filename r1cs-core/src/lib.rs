//! Core interface for working with Rank-1 Constraint Systems (R1CS).

#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, private_in_public, variant_size_differences)]
#![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
#![deny(unused_attributes, unused_imports, unused_mut, missing_docs)]
#![deny(renamed_and_removed_lints, stable_features, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use, const_err)]

#![forbid(unsafe_code)]

mod constraint_system;
mod error;
mod impl_lc;
mod impl_constraint_var;

pub use constraint_system::{ConstraintSystem, ConstraintSynthesizer, Namespace};
pub use error::SynthesisError;
pub use algebra::ToConstraintField;

use algebra::Field;
use smallvec::SmallVec as StackVec;
use std::cmp::Ordering;

type SmallVec<F> = StackVec<[(Variable, F); 16]>;

/// Represents a variable in a constraint system.
#[derive(PartialOrd, Ord, PartialEq, Eq, Copy, Clone, Debug)]
pub struct Variable(Index);

impl Variable {
    /// This constructs a variable with an arbitrary index.
    /// Circuit implementations are not recommended to use this.
    pub fn new_unchecked(idx: Index) -> Variable {
        Variable(idx)
    }

    /// This returns the index underlying the variable.
    /// Circuit implementations are not recommended to use this.
    pub fn get_unchecked(&self) -> Index {
        self.0
    }
}

/// Represents the index of either an input variable or auxiliary variable.
#[derive(Copy, Clone, PartialEq, Debug, Eq)]
pub enum Index {
    /// Index of an input variable.
    Input(usize),
    /// Index of an auxiliary (or private) variable.
    Aux(usize),
}

impl PartialOrd for Index {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Index {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Index::Input(ref idx1), Index::Input(ref idx2))
            | (Index::Aux(ref idx1), Index::Aux(ref idx2)) => idx1.cmp(idx2),
            (Index::Input(_), Index::Aux(_)) => Ordering::Less,
            (Index::Aux(_), Index::Input(_)) => Ordering::Greater,
        }
    }
}


/// This represents a linear combination of some variables, with coefficients
/// in the field `F`.
/// The `(coeff, var)` pairs in a `LinearCombination` are kept sorted according
/// to the index of the variable in its constraint system.
#[derive(Debug, Clone)]
pub struct LinearCombination<F: Field>(pub SmallVec<F>);


/// Either a `Variable` or a `LinearCombination`.
#[derive(Clone, Debug)]
pub enum ConstraintVar<F: Field> {
    /// A wrapper around a `LinearCombination`.
    LC(LinearCombination<F>),
    /// A wrapper around a `Variable`.
    Var(Variable),
}
