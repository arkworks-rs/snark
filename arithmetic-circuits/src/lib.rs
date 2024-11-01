//! A library for arithmetic circuits

#![cfg_attr(not(any(feature = "std", test)), no_std)]
#![warn(
    unused,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms,
    missing_docs
)]

pub mod arithmetic_circuit;
pub mod expression;

#[cfg(test)]
pub(crate) mod reader;

/// The path to the test data directory that contains the circom test files.
#[macro_export]
macro_rules! TEST_DATA_PATH {
    () => {
        concat!(env!("CARGO_MANIFEST_DIR"), "/../circom/{}",)
    };
}
