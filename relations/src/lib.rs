//! Core interface for working with various relations that are useful in
//! zkSNARKs. At the moment, we implement APIs for working with Rank-1
//! Constraint Systems (R1CS), Arithmetic Circuits and Arithmetic Expressions.

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

pub mod arithmetic_circuit;
pub mod expression;
pub mod r1cs;

#[cfg(test)]
pub(crate) mod reader;

/// The path to the test data directory that contains the circom test files.
#[macro_export]
macro_rules! TEST_DATA_PATH {
    () => {
        concat!(env!("CARGO_MANIFEST_DIR"), "/../circom/{}",)
    };
}
