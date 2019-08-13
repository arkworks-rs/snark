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

#[macro_use]
extern crate algebra;
#[macro_use]
extern crate derivative;

pub mod test_constraint_system;

pub mod bits;
pub use self::bits::*;

pub mod fields;

pub mod groups;

pub mod pairing;

pub mod eq;
pub mod select;
pub mod alloc;

pub mod prelude {
    pub use crate::eq::*;
    pub use crate::select::*;
    pub use crate::alloc::*;
    pub use crate::fields::FieldGadget;
    pub use crate::groups::GroupGadget;
    pub use crate::pairing::PairingGadget;
    pub use crate::bits::{ToBitsGadget, ToBytesGadget, boolean::Boolean, uint8::UInt8, uint32::UInt32};
}

pub trait Assignment<T> {
    fn get(self) -> Result<T, r1cs_core::SynthesisError>;
}

impl<T> Assignment<T> for Option<T> {
    fn get(self) -> Result<T, r1cs_core::SynthesisError> {
        match self {
            Some(v) => Ok(v),
            None => Err(r1cs_core::SynthesisError::AssignmentMissing),
        }
    }
}
