#![cfg_attr(not(feature = "std"), no_std)]
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

#[cfg(all(test, not(feature = "std")))]
#[macro_use]
extern crate std;

#[cfg(not(feature = "std"))]
extern crate alloc as ralloc;

#[macro_use]
extern crate algebra;

#[macro_use]
extern crate derivative;

/// used by test_constraint_system
#[cfg(not(feature = "std"))]
macro_rules! println {
    () => {};
    ($($arg: tt)*) => {};
}

#[cfg(not(feature = "std"))]
use ralloc::{collections::BTreeMap, string::String, vec::Vec};

#[cfg(feature = "std")]
use std::{collections::BTreeMap, string::String, vec::Vec};

pub mod test_constraint_system;

pub mod bits;
pub use self::bits::*;

pub mod fields;

pub mod groups;

pub mod pairing;

pub mod alloc;
pub mod eq;
pub mod select;

pub mod prelude {
    pub use crate::{
        alloc::*,
        bits::{boolean::Boolean, uint32::UInt32, uint8::UInt8, ToBitsGadget, ToBytesGadget},
        eq::*,
        fields::FieldGadget,
        groups::GroupGadget,
        pairing::PairingGadget,
        select::*,
    };
}

pub trait Assignment<T> {
    fn get(self) -> Result<T, r1cs_core::SynthesisError>;
}

impl<T> Assignment<T> for Option<T> {
    fn get(self) -> Result<T, r1cs_core::SynthesisError> {
        self.ok_or_else(|| r1cs_core::SynthesisError::AssignmentMissing)
    }
}
