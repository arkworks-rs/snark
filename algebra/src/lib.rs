#![cfg_attr(not(feature = "std"), no_std)]
#![deny(unused_import_braces, trivial_casts, trivial_numeric_casts)]
#![deny(unused_qualifications, variant_size_differences, unused_extern_crates)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(renamed_and_removed_lints, unused_allocation, unused_comparisons)]
#![deny(const_err, unused_must_use, unused_mut, bare_trait_objects)]
#![forbid(unsafe_code)]

#[cfg(all(test, not(feature = "std")))]
#[macro_use]
extern crate std;

/// this crate needs to be public, cause we expose `to_bytes!` macro
/// see similar issue in [`smallvec#198`]
///
/// [`smallvec#198`]: https://github.com/servo/rust-smallvec/pull/198
#[cfg(not(feature = "std"))]
#[allow(unused_imports)]
#[macro_use]
#[doc(hidden)]
pub extern crate alloc;

#[cfg(not(feature = "std"))]
#[allow(unused_imports)]
#[doc(hidden)]
pub use alloc::{boxed::Box, format, string::String, vec, vec::Vec};

#[cfg(feature = "std")]
#[allow(unused_imports)]
#[doc(hidden)]
pub use std::{boxed::Box, format, vec, vec::Vec};

pub use algebra_core::*;

///////////////////////////////////////////////////////////////////////////////
#[cfg(feature = "bls12_377")]
pub mod bls12_377;
#[cfg(feature = "bls12_377")]
pub use bls12_377::Bls12_377;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(
    not(feature = "bls12_377"),
    any(feature = "edwards_bls12", feature = "sw6", feature = "edwards_sw6")
))]
pub(crate) mod bls12_377;

#[cfg(feature = "edwards_bls12")]
pub mod edwards_bls12;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(feature = "bls12_381")]
pub mod bls12_381;
#[cfg(feature = "bls12_381")]
pub use bls12_381::Bls12_381;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(not(feature = "bls12_381"), feature = "jubjub"))]
pub(crate) mod bls12_381;

#[cfg(feature = "jubjub")]
pub mod jubjub;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(not(feature = "mnt6_298"), feature = "mnt4_298"))]
pub(crate) mod mnt6_298;

#[cfg(feature = "mnt4_298")]
pub mod mnt4_298;
#[cfg(feature = "mnt4_298")]
pub use mnt4_298::MNT4_298;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(not(feature = "mnt6_753"), feature = "mnt4_753"))]
pub(crate) mod mnt6_753;

#[cfg(feature = "mnt4_753")]
pub mod mnt4_753;
#[cfg(feature = "mnt4_753")]
pub use mnt4_753::MNT4_753;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(feature = "mnt6_298")]
pub mod mnt6_298;
#[cfg(feature = "mnt6_298")]
pub use mnt6_298::MNT6_298;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(feature = "mnt6_753")]
pub mod mnt6_753;
#[cfg(feature = "mnt6_753")]
pub use mnt6_753::MNT6_753;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(feature = "sw6")]
pub mod sw6;
#[cfg(feature = "sw6")]
pub use sw6::SW6;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(not(feature = "sw6"), feature = "edwards_sw6"))]
pub(crate) mod sw6;

#[cfg(feature = "edwards_sw6")]
pub mod edwards_sw6;
///////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
pub(crate) mod tests;
