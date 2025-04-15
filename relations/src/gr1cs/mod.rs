//! Core interface for working with Generalized Rank-1 Constraint Systems
//! (GR1CS).
mod constraint_system_ref;

/// Interface for specifying strategies for reducing the number of constraints
/// that public input/instance variables are involved in.
pub mod instance_outliner;

mod namespace;

pub mod predicate;

#[macro_use]
mod constraint_system;

#[cfg(test)]
mod tests;

#[cfg(feature = "std")]
pub mod trace;
///////////////////////////////////////////////////////////////////////////////////////

#[cfg(feature = "std")]
pub use crate::gr1cs::trace::{ConstraintLayer, ConstraintTrace, TraceStep, TracingMode};

use ark_std::vec::Vec;
pub use tracing::info_span;

pub use ark_ff::{Field, ToConstraintField};

pub use crate::{
    gr1cs::{
        constraint_system::ConstraintSystem, constraint_system_ref::ConstraintSystemRef,
        predicate::polynomial_constraint::R1CS_PREDICATE_LABEL,
    },
    lc,
    utils::{
        error::SynthesisError,
        linear_combination::{LcIndex, LinearCombination},
        matrix::{mat_vec_mul, transpose, Matrix},
        variable::Variable,
        Result,
    },
};
pub use namespace::Namespace;
///////////////////////////////////////////////////////////////////////////////////////

/// Computations are expressed in terms of generalized rank-1 constraint systems
/// (GR1CS). The `generate_constraints` method is called to generate constraints
/// for both CRS generation and for proving.
// TODO: Think: should we replace this with just a closure?
pub trait ConstraintSynthesizer<F: Field> {
    /// Drives generation of new constraints inside `cs`.
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> crate::gr1cs::Result<()>;
}

///////////////////////////////////////////////////////////////////////////////////////

/// In GR1CS a constraint is a vector of linear combinations associated with a
///  predicate
pub type Constraint = Vec<LcIndex>;

/// Each predicate is associated with a label
pub type Label = ark_std::string::String;
///////////////////////////////////////////////////////////////////////////////////////

/// Defines the mode of operation of a `ConstraintSystem`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum SynthesisMode {
    /// Indicate to the `ConstraintSystem` that it should only generate
    /// constraint matrices and not populate the variable assignments.
    Setup,
    /// Indicate to the `ConstraintSystem` that it populate the variable
    /// assignments. If additionally `construct_matrices == true`, then generate
    /// the matrices as in the `Setup` case.
    Prove {
        /// If `construct_matrices == true`, then generate
        /// the matrices as in the `Setup` case.
        construct_matrices: bool,
        /// If `construct_matrices == true`, then generate
        /// the matrices as in the `Setup` case.
        generate_lc_assignments: bool,
    },
}

///////////////////////////////////////////////////////////////////////////////////////

/// Defines the parameter to optimize for a `ConstraintSystem`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum OptimizationGoal {
    /// Make no attempt to optimize.
    None,
    /// Minimize the number of constraints by inlining the linear combinations.
    Constraints,
    /// Minimize the total weight of the constraints (the number of nonzero
    /// entries across all constraints) by outlining the linear combinations
    /// and creating new witness variables.
    #[deprecated]
    Weight,
}

///////////////////////////////////////////////////////////////////////////////////////
