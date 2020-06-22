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
#[cfg(feature = "bn254")]
pub mod bn254;
#[cfg(feature = "bn254")]
pub use bn254::Bn254;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(not(feature = "bn254"), feature = "ed_on_bn254"))]
pub(crate) mod bn254;

#[cfg(feature = "ed_on_bn254")]
pub mod ed_on_bn254;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(feature = "bls12_377")]
pub mod bls12_377;
#[cfg(feature = "bls12_377")]
pub use bls12_377::Bls12_377;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(
    not(feature = "bls12_377"),
    any(
        feature = "ed_on_bls12_377",
        feature = "cp6_782",
        feature = "ed_on_cp6_782",
        feature = "ed_on_bw6_761",
        feature = "bw6_761",
    )
))]
pub(crate) mod bls12_377;

#[cfg(feature = "ed_on_bls12_377")]
pub mod ed_on_bls12_377;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(feature = "bls12_381")]
pub mod bls12_381;
#[cfg(feature = "bls12_381")]
pub use bls12_381::Bls12_381;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(not(feature = "bls12_381"), feature = "ed_on_bls12_381"))]
pub(crate) mod bls12_381;

#[cfg(feature = "ed_on_bls12_381")]
pub mod ed_on_bls12_381;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(
    not(feature = "mnt6_298"),
    any(feature = "mnt4_298", feature = "ed_on_mnt4_298")
))]
pub(crate) mod mnt6_298;

#[cfg(feature = "mnt4_298")]
pub mod mnt4_298;
#[cfg(feature = "mnt4_298")]
pub use mnt4_298::MNT4_298;

#[cfg(all(not(feature = "mnt4_298"), feature = "ed_on_mnt4_298"))]
pub(crate) mod mnt4_298;

#[cfg(feature = "ed_on_mnt4_298")]
pub mod ed_on_mnt4_298;

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(
    not(feature = "mnt6_753"),
    any(feature = "mnt4_753", feature = "ed_on_mnt4_753")
))]
pub(crate) mod mnt6_753;

#[cfg(feature = "mnt4_753")]
pub mod mnt4_753;
#[cfg(feature = "mnt4_753")]
pub use mnt4_753::MNT4_753;

#[cfg(all(not(feature = "mnt4_753"), feature = "ed_on_mnt4_753"))]
pub(crate) mod mnt4_753;

#[cfg(feature = "ed_on_mnt4_753")]
pub mod ed_on_mnt4_753;
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
#[cfg(feature = "cp6_782")]
pub mod cp6_782;
#[cfg(feature = "cp6_782")]
pub use cp6_782::CP6_782;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(feature = "bw6_761")]
pub mod bw6_761;
#[cfg(feature = "bw6_761")]
pub use bw6_761::BW6_761;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(
    not(feature = "cp6_782"),
    any(feature = "ed_on_cp6_782", feature = "ed_on_bw6_761")
))]
pub(crate) mod cp6_782;

#[cfg(feature = "ed_on_cp6_782")]
pub mod ed_on_cp6_782;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(not(feature = "ed_on_cp6_782"), feature = "ed_on_bw6_761"))]
pub(crate) mod ed_on_cp6_782;

#[cfg(feature = "ed_on_bw6_761")]
pub mod ed_on_bw6_761;
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#[cfg(all(not(feature = "bw6_761"), feature = "ed_on_cp6_782"))]
pub(crate) mod bw6_761;
///////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
pub(crate) mod tests;
