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
#![allow(
    clippy::upper_case_acronyms,
    clippy::too_many_arguments,
    clippy::type_complexity,
    clippy::try_err,
    clippy::map_collect_result_unit,
    clippy::not_unsafe_ptr_arg_deref,
    clippy::suspicious_op_assign_impl,
    clippy::suspicious_arithmetic_impl,
    clippy::assertions_on_constants
)]

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

pub mod instantiated;
pub use instantiated::*;

pub mod alloc;
pub mod eq;
pub mod select;
pub mod to_field_gadget_vec;

pub mod prelude {
    pub use crate::{
        alloc::*,
        bits::{
            boolean::Boolean, uint32::UInt32, uint8::UInt8, FromBitsGadget, ToBitsGadget,
            ToBytesGadget,
        },
        eq::*,
        fields::{cubic_extension::*, quadratic_extension::*, FieldGadget},
        groups::GroupGadget,
        pairing::PairingGadget,
        select::*,
    };
}

use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};

pub trait Assignment<T> {
    fn get(self) -> Result<T, SynthesisError>;
}

impl<T> Assignment<T> for Option<T> {
    fn get(self) -> Result<T, SynthesisError> {
        match self {
            Some(v) => Ok(v),
            None => Err(SynthesisError::AssignmentMissing),
        }
    }
}

pub trait FromGadget<T, ConstraintF: Field>: Sized {
    fn from<CS: ConstraintSystem<ConstraintF>>(other: T, cs: CS) -> Result<Self, SynthesisError>;
}
