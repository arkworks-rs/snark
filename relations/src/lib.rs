//! Core interface for working with various relations that are useful in
//! zkSNARKs. At the moment, we only implement APIs for working with Rank-1
//! Constraint Systems (R1CS).

#![cfg_attr(not(feature = "std"), no_std)]
#![warn(
    unused,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
    missing_docs
)]
#![deny(unsafe_code)]

#[macro_use]
extern crate ark_std;

/// The Generalized Rank-1 Constraint System (GR1CS) Infrastructure
pub mod gr1cs;

/// Utilities functions and data structures needed for working with GR1CS
pub mod utils;
