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

#[macro_use]
pub mod macros;

#[cfg(not(feature = "std"))]
use ralloc::vec::Vec;

#[cfg(feature = "std")]
use std::vec::Vec;

use algebra::prelude::Field;

pub mod bits;
pub use self::bits::*;

pub mod fields;

pub mod groups;

mod instantiated;

#[cfg(feature = "bls12_377")]
pub use instantiated::bls12_377;

#[cfg(feature = "ed_on_bn254")]
pub use instantiated::ed_on_bn254;

#[cfg(feature = "ed_on_bls12_377")]
pub use instantiated::ed_on_bls12_377;

#[cfg(feature = "ed_on_mnt4_298")]
pub use instantiated::ed_on_mnt4_298;

#[cfg(feature = "ed_on_mnt4_753")]
pub use instantiated::ed_on_mnt4_753;

#[cfg(feature = "ed_on_cp6_782")]
pub use instantiated::ed_on_cp6_782;

#[cfg(feature = "ed_on_bw6_761")]
pub use instantiated::ed_on_bw6_761;

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
        fields::{FieldOpsBounds, FieldVar},
        groups::{CurveVar, GroupOpsBounds},
        instantiated::*,
        pairing::PairingVar,
        select::*,
        R1CSVar,
    };
}

pub trait R1CSVar<F: Field> {
    type Value: core::fmt::Debug + Eq + Clone;

    /// Returns the underlying `ConstraintSystemRef`.
    fn cs(&self) -> Option<r1cs_core::ConstraintSystemRef<F>>;

    /// Returns `true` if `self` is a circuit-generation-time constant.
    fn is_constant(&self) -> bool {
        self.cs()
            .map_or(true, |cs| cs == r1cs_core::ConstraintSystemRef::None)
    }

    /// Returns the value that is assigned to `self` in the underlying
    /// `ConstraintSystem`.
    fn value(&self) -> Result<Self::Value, r1cs_core::SynthesisError>;
}

impl<F: Field, T: R1CSVar<F>> R1CSVar<F> for [T] {
    type Value = Vec<T::Value>;

    fn cs(&self) -> Option<r1cs_core::ConstraintSystemRef<F>> {
        let mut result = None;
        for var in self {
            result = var.cs().or(result);
        }
        result
    }

    fn value(&self) -> Result<Self::Value, r1cs_core::SynthesisError> {
        let mut result = Vec::new();
        for var in self {
            result.push(var.value()?);
        }
        Ok(result)
    }
}

impl<'a, F: Field, T: 'a + R1CSVar<F>> R1CSVar<F> for &'a T {
    type Value = T::Value;

    fn cs(&self) -> Option<r1cs_core::ConstraintSystemRef<F>> {
        (*self).cs()
    }

    fn value(&self) -> Result<Self::Value, r1cs_core::SynthesisError> {
        (*self).value()
    }
}

pub trait Assignment<T> {
    fn get(self) -> Result<T, r1cs_core::SynthesisError>;
}

impl<T> Assignment<T> for Option<T> {
    fn get(self) -> Result<T, r1cs_core::SynthesisError> {
        self.ok_or(r1cs_core::SynthesisError::AssignmentMissing)
    }
}
