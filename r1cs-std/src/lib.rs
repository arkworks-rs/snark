#![cfg_attr(not(feature = "std"), no_std)]
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, variant_size_differences, unreachable_pub)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(unused_extern_crates, renamed_and_removed_lints, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, const_err, unused_must_use)]
#![deny(unused_mut, unused_unsafe, private_in_public, unsafe_code)]
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

pub mod test_constraint_counter;
pub mod test_constraint_system;

pub mod bits;
pub use self::bits::*;

pub mod fields;

pub mod groups;

mod instantiated;

#[cfg(feature = "bls12_377")]
pub use instantiated::bls12_377;

#[cfg(feature = "ed_on_bls12_377")]
pub use instantiated::ed_on_bls12_377;

#[cfg(feature = "ed_on_cp6_782")]
pub use instantiated::ed_on_cp6_782;

#[cfg(feature = "ed_on_bls12_381")]
pub use instantiated::ed_on_bls12_381;

#[cfg(feature = "mnt4_298")]
pub use instantiated::mnt4_298;

#[cfg(feature = "mnt4_753")]
pub use instantiated::mnt4_753;

#[cfg(feature = "mnt6_298")]
pub use instantiated::mnt6_298;

#[cfg(feature = "mnt6_753")]
pub use instantiated::mnt6_753;

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
        instantiated::*,
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
