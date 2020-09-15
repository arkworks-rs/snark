#![cfg_attr(not(feature = "std"), no_std)]
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, variant_size_differences)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(unused_extern_crates, renamed_and_removed_lints, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, const_err, unused_must_use)]
#![deny(unused_mut, unused_unsafe, private_in_public)]
#![cfg_attr(use_asm, feature(llvm_asm))]
#![cfg_attr(not(use_asm), forbid(unsafe_code))]
#![cfg_attr(use_asm, deny(unsafe_code))]

#[cfg(all(test, not(feature = "std")))]
#[macro_use]
extern crate std;

/// This crate needs to be public, because we expose the `to_bytes!` macro.
/// See the similar issue in [`smallvec#198`]
///
/// [`smallvec#198`]: https://github.com/servo/rust-smallvec/pull/198
#[cfg(not(feature = "std"))]
#[macro_use]
#[doc(hidden)]
pub extern crate alloc;

#[cfg(not(feature = "std"))]
#[allow(unused_imports)]
#[doc(hidden)]
pub use alloc::{
    borrow::{Cow, ToOwned},
    boxed::Box,
    collections::btree_map::BTreeMap,
    format,
    string::String,
    vec,
    vec::Vec,
};

#[cfg(feature = "std")]
#[allow(unused_imports)]
#[doc(hidden)]
pub use std::{
    borrow::{Cow, ToOwned},
    boxed::Box,
    collections::btree_map::BTreeMap,
    format,
    string::String,
    vec,
    vec::Vec,
};

#[macro_use]
extern crate derivative;

#[cfg_attr(test, macro_use)]
pub mod bytes;
pub use self::bytes::*;

#[macro_use]
pub mod serialize;
pub use self::serialize::*;

#[macro_use]
pub mod fields;
pub use self::fields::*;

pub mod biginteger;
pub use self::biginteger::*;

pub mod curves;
pub use self::curves::*;

pub mod groups;
pub use self::groups::*;

mod rand;
pub use self::rand::*;

mod error;
pub use self::error::*;

mod to_field_vec;
pub use to_field_vec::ToConstraintField;

pub mod msm;
pub use self::msm::*;

pub use num_traits::{One, Zero};

pub mod prelude {
    pub use crate::biginteger::BigInteger;

    pub use crate::fields::{Field, FpParameters, PrimeField, SquareRootField};

    pub use crate::groups::Group;

    pub use crate::curves::{AffineCurve, PairingEngine, ProjectiveCurve};

    pub use crate::rand::UniformRand;

    pub use num_traits::{One, Zero};

    pub use crate::error::*;
}

#[cfg(not(feature = "std"))]
pub mod io;

#[cfg(feature = "std")]
pub use std::io;

#[cfg(feature = "derive")]
#[allow(unused_imports)]
#[macro_use]
extern crate algebra_core_derive;

#[cfg(not(feature = "std"))]
fn error(_msg: &'static str) -> io::Error {
    io::Error
}

#[cfg(feature = "std")]
fn error(msg: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::Other, msg)
}

/// Returns the base-2 logarithm of `x`.
///
/// ```
/// use algebra_core::log2;
///
/// assert_eq!(log2(16), 4);
/// assert_eq!(log2(17), 5);
/// assert_eq!(log2(1), 0);
/// assert_eq!(log2(0), 0);
/// assert_eq!(log2(usize::MAX), (core::mem::size_of::<usize>() * 8) as u32);
/// assert_eq!(log2(1 << 15), 15);
/// assert_eq!(log2(2usize.pow(18)), 18);
/// ```
pub fn log2(x: usize) -> u32 {
    if x == 0 {
        0
    } else if x.is_power_of_two() {
        1usize.leading_zeros() - x.leading_zeros()
    } else {
        0usize.leading_zeros() - x.leading_zeros()
    }
}
