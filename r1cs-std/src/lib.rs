//! This crate implements common "gadgets" that make
//! programming rank-1 constraint systems easier.
#![cfg_attr(not(feature = "std"), no_std)]
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, variant_size_differences, unreachable_pub)]
#![deny(non_shorthand_field_patterns, unused_attributes, unused_imports)]
#![deny(unused_extern_crates, renamed_and_removed_lints, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, const_err, unused_must_use)]
#![deny(unused_mut, unused_unsafe, private_in_public, unsafe_code)]
#![deny(missing_docs)]
#![forbid(unsafe_code)]

#[doc(hidden)]
#[cfg(all(test, not(feature = "std")))]
#[macro_use]
extern crate std;

#[doc(hidden)]
#[cfg(not(feature = "std"))]
extern crate alloc as ralloc;

#[doc(hidden)]
#[macro_use]
extern crate algebra;

#[doc(hidden)]
#[macro_use]
extern crate derivative;

/// Some utility macros for making downstream impls easier.
#[macro_use]
pub mod macros;

#[cfg(not(feature = "std"))]
use ralloc::vec::Vec;

#[cfg(feature = "std")]
use std::vec::Vec;

use algebra::prelude::Field;

/// This module implements gadgets related to bit manipulation, such as
/// `Boolean` and `UInt`s.
pub mod bits;
pub use self::bits::*;

/// This module implements gadgets related to field arithmetic.
pub mod fields;

/// This module implements gadgets related to group arithmetic, and specifically
/// elliptic curve arithmetic.
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

/// This module implements gadgets related to computing pairings in bilinear
/// groups.
pub mod pairing;

/// This module describes a trait for allocating new variables in a constraint
/// system.
pub mod alloc;
/// This module describes a trait for checking equality of variables.
pub mod eq;
/// This module describes traits for conditionally selecting a variable from a
/// list of variables.
pub mod select;

#[allow(missing_docs)]
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

/// This trait describes some core functionality that is common to high-level
/// variables, such as `Boolean`s, `FieldVar`s, `GroupVar`s, etc.
pub trait R1CSVar<F: Field> {
    /// The type of the "native" value that `Self` represents in the constraint
    /// system.
    type Value: core::fmt::Debug + Eq + Clone;

    /// Returns the underlying `ConstraintSystemRef`.
    ///
    /// If `self` is a constant value, then this *must* return
    /// `r1cs_core::ConstraintSystemRef::None`.
    fn cs(&self) -> r1cs_core::ConstraintSystemRef<F>;

    /// Returns `true` if `self` is a circuit-generation-time constant.
    fn is_constant(&self) -> bool {
        self.cs().is_none()
    }

    /// Returns the value that is assigned to `self` in the underlying
    /// `ConstraintSystem`.
    fn value(&self) -> Result<Self::Value, r1cs_core::SynthesisError>;
}

impl<F: Field, T: R1CSVar<F>> R1CSVar<F> for [T] {
    type Value = Vec<T::Value>;

    fn cs(&self) -> r1cs_core::ConstraintSystemRef<F> {
        let mut result = r1cs_core::ConstraintSystemRef::None;
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

    fn cs(&self) -> r1cs_core::ConstraintSystemRef<F> {
        (*self).cs()
    }

    fn value(&self) -> Result<Self::Value, r1cs_core::SynthesisError> {
        (*self).value()
    }
}

/// A utility trait to convert `Self` to `Result<T, SynthesisErrorA`.>
pub trait Assignment<T> {
    /// Converts `self` to `Result`.
    fn get(self) -> Result<T, r1cs_core::SynthesisError>;
}

impl<T> Assignment<T> for Option<T> {
    fn get(self) -> Result<T, r1cs_core::SynthesisError> {
        self.ok_or(r1cs_core::SynthesisError::AssignmentMissing)
    }
}

/// Specifies how to convert a variable of type `Self` to variables of
/// type `FpVar<ConstraintF>`
pub trait ToConstraintFieldGadget<ConstraintF: algebra::PrimeField> {
    /// Converts `self` to `FpVar<ConstraintF>` variables.
    fn to_constraint_field(
        &self,
    ) -> Result<Vec<crate::fields::fp::FpVar<ConstraintF>>, r1cs_core::SynthesisError>;
}
