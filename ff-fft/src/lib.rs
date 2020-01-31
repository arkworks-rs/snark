//! This crate implements functions for manipulating polynomials over finite
//! fields, including FFTs.
// this crate is not yet no_std
#![no_std]
#![deny(unused_import_braces, trivial_casts, bare_trait_objects, missing_docs)]
#![deny(unused_qualifications, variant_size_differences, stable_features)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(renamed_and_removed_lints, unused_allocation, unused_comparisons)]
#![deny(const_err, unused_must_use, unused_mut, private_in_public)]
#![deny(unreachable_pub, unused_extern_crates, trivial_numeric_casts)]
#![forbid(unsafe_code)]

#[cfg(not(feature = "std"))]
extern crate alloc;

#[macro_use]
extern crate std;

#[cfg(not(feature = "std"))]
pub(crate) use alloc::{borrow::Cow, collections::BTreeMap, vec::Vec};

#[cfg(feature = "std")]
pub(crate) use std::{borrow::Cow, collections::BTreeMap, vec::Vec};

pub mod domain;
pub mod evaluations;
pub mod polynomial;

pub use domain::EvaluationDomain;
pub use evaluations::Evaluations;
pub use polynomial::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial};

#[cfg(test)]
mod test;
