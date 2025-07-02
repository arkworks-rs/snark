//! Core interface for working with various relations that are useful in
//! zkSNARKs. At the moment, we only implement APIs for working with Generalized
//! Rank-1 Constraint Systems (R1CS) (See <https://eprint.iacr.org/2024/1245.pdf>).
//!
//! # Compatibility with R1CS
//!
//! In previous versions, this crate only supported R1CS.
//! For ease of migration, a GR1CS instance is equipped with an R1CS predicate by default.
//! Also, there is a separate
//! [`enforce_r1cs_constraint`](`crate::gr1cs::ConstraintSystemRef::enforce_r1cs_constraint`) function that has the same API as the `enforce_constraint` function
//! in previous versions.
//! Hence, by replacing the latter with the former and updating `use` statements,
//! you can migrate your code to the new version.

#![cfg_attr(not(feature = "std"), no_std)]
#![warn(
    unused,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
    missing_docs,
    clippy::pedantic
)]
#![allow(
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::must_use_candidate,
    clippy::inline_always,
    clippy::redundant_closure_for_method_calls
)]
#![deny(unsafe_code)]

#[macro_use]
extern crate ark_std;

/// The Generalized Rank-1 Constraint System (GR1CS) Infrastructure
pub mod gr1cs;

/// The Squared Rank-1 Constraint System (GR1CS) Infrastructure
pub mod sr1cs;

/// Functions and data structures needed for working with GR1CS
pub mod utils;
