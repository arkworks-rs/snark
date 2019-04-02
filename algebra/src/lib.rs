#![deny(
    unused_import_braces,
    unused_qualifications,
    trivial_casts,
    trivial_numeric_casts
)]
#![deny(
    unused_qualifications,
    variant_size_differences,
    stable_features,
    unreachable_pub
)]
#![deny(
    non_shorthand_field_patterns,
    unused_attributes,
    unused_imports,
    unused_extern_crates
)]
#![deny(
    renamed_and_removed_lints,
    stable_features,
    unused_allocation,
    unused_comparisons,
    bare_trait_objects
)]
#![deny(
    const_err,
    unused_must_use,
    unused_mut,
    unused_unsafe,
    private_in_public,
    unsafe_code
)]
#![forbid(unsafe_code)]
#![feature(const_fn)]

#[macro_use]
extern crate derivative;

#[cfg_attr(test, macro_use)]
pub mod bytes;
pub use self::bytes::*;

pub mod biginteger;
pub use self::biginteger::*;

pub mod curves;
pub use self::curves::*;

pub mod fields;
pub use self::fields::*;

pub mod groups;
pub use self::groups::*;

#[cfg(feature = "parallel")]
pub mod msm;
#[cfg(feature = "parallel")]
pub use self::msm::*;

#[cfg(feature = "parallel")]
pub mod fft;
#[cfg(feature = "parallel")]
pub use self::fft::*;

pub mod utils;
