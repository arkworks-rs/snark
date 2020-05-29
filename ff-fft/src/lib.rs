//! This crate implements functions for manipulating polynomials over finite
//! fields, including FFTs.
#![cfg_attr(not(feature = "std"), no_std)]
#![deny(unused_import_braces, trivial_casts, bare_trait_objects, missing_docs)]
#![deny(unused_qualifications, variant_size_differences, stable_features)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(renamed_and_removed_lints, unused_allocation, unused_comparisons)]
#![deny(const_err, unused_must_use, unused_mut, private_in_public)]
#![deny(unreachable_pub, unused_extern_crates, trivial_numeric_casts)]
#![forbid(unsafe_code)]

#[cfg(not(feature = "std"))]
#[macro_use]
extern crate alloc;

#[cfg(not(feature = "std"))]
pub(crate) use alloc::{borrow::Cow, collections::BTreeMap, vec::Vec};

#[cfg(feature = "std")]
pub(crate) use std::{borrow::Cow, collections::BTreeMap, vec::Vec};

/// Creates parallel iterator over refs if `parallel` feature is enabled.
#[macro_export]
macro_rules! cfg_iter {
    ($e: expr) => {{
        #[cfg(feature = "parallel")]
        let result = $e.par_iter();

        #[cfg(not(feature = "parallel"))]
        let result = $e.iter();

        result
    }};
}

/// Creates parallel iterator over mut refs if `parallel` feature is enabled.
#[macro_export]
macro_rules! cfg_iter_mut {
    ($e: expr) => {{
        #[cfg(feature = "parallel")]
        let result = $e.par_iter_mut();

        #[cfg(not(feature = "parallel"))]
        let result = $e.iter_mut();

        result
    }};
}

/// Creates parallel iterator if `parallel` feature is enabled.
#[macro_export]
macro_rules! cfg_into_iter {
    ($e: expr) => {{
        #[cfg(feature = "parallel")]
        let result = $e.into_par_iter();

        #[cfg(not(feature = "parallel"))]
        let result = $e.into_iter();

        result
    }};
}

/// Returns an iterator over `chunk_size` elements of the slice at a
/// time.
#[macro_export]
macro_rules! cfg_chunks {
    ($e: expr, $size: expr) => {{
        #[cfg(feature = "parallel")]
        let result = $e.par_chunks($size);

        #[cfg(not(feature = "parallel"))]
        let result = $e.chunks($size);

        result
    }};
}

/// Returns an iterator over `chunk_size` mutable elements of the slice at a
/// time.
#[macro_export]
macro_rules! cfg_chunks_mut {
    ($e: expr, $size: expr) => {{
        #[cfg(feature = "parallel")]
        let result = $e.par_chunks_mut($size);

        #[cfg(not(feature = "parallel"))]
        let result = $e.chunks_mut($size);

        result
    }};
}

pub mod domain;

pub mod evaluations;
pub mod polynomial;

pub use domain::{
    EvaluationDomain, GeneralEvaluationDomain, MixedRadixEvaluationDomain, Radix2EvaluationDomain,
};
pub use evaluations::Evaluations;
pub use polynomial::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial};

#[cfg(test)]
mod test;
