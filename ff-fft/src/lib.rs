//! This crate implements functions for manipulating polynomials over finite fields,
//! including FFTs.
#![deny(unused_import_braces, trivial_casts, bare_trait_objects, missing_docs)]
#![deny(unused_qualifications, variant_size_differences, stable_features)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(renamed_and_removed_lints, unused_allocation, unused_comparisons)]
#![deny(const_err, unused_must_use, unused_mut, private_in_public)]
#![deny(unreachable_pub, unused_extern_crates, trivial_numeric_casts)]
#![forbid(unsafe_code)]

pub mod domain;
pub mod evaluations;
pub mod polynomial;

pub(crate) mod multicore;

pub use domain::EvaluationDomain;
pub use evaluations::Evaluations;
pub use polynomial::{DensePolynomial, SparsePolynomial, DenseOrSparsePolynomial};

#[cfg(test)]
mod test;
