//! This module implements the R1CS equivalent of `algebra::ed_on_bls12_377`.
//!
//! It implements field variables for `algebra::ed_on_bls12_377::Fq`,
//! and group variables for `algebra::ed_on_bls12_377::GroupProjective`.
//!
//! The field underlying these constraints is `algebra::ed_on_bls12_377::Fq`.
//!
//! # Examples
//!
//! One can perform standard algebraic operations on `FqVar`:
//!
//! ```
//! # fn main() -> Result<(), r1cs_core::SynthesisError> {
//! use algebra::{UniformRand, ed_on_bls12_377::*};
//! use r1cs_core::*;
//! use r1cs_std::prelude::*;
//! use r1cs_std::ed_on_bls12_377::*;
//!
//! let cs = ConstraintSystem::<Fq>::new_ref();
//! // This rng is just for test purposes; do not use it
//! // in real applications.
//! let mut rng = algebra::test_rng();
//!
//! // Generate some random `Fq` elements.
//! let a_native = Fq::rand(&mut rng);
//! let b_native = Fq::rand(&mut rng);
//!
//! // Allocate `a_native` and `b_native` as witness variables in `cs`.
//! let a = FqVar::new_witness(r1cs_core::ns!(cs, "generate_a"), || Ok(a_native))?;
//! let b = FqVar::new_witness(r1cs_core::ns!(cs, "generate_b"), || Ok(b_native))?;
//!
//! // Allocate `a_native` and `b_native` as constants in `cs`. This does not add any
//! // constraints or variables.
//! let a_const = FqVar::new_constant(r1cs_core::ns!(cs, "a_as_constant"), a_native)?;
//! let b_const = FqVar::new_constant(r1cs_core::ns!(cs, "b_as_constant"), b_native)?;
//!
//! let one = FqVar::one();
//! let zero = FqVar::zero();
//!
//! // Sanity check one + one = two
//! let two = &one + &one + &zero;
//! two.enforce_equal(&one.double()?)?;
//!
//! assert!(cs.is_satisfied()?);
//!
//! // Check that the value of &a + &b is correct.
//! assert_eq!((&a + &b).value()?, a_native + &b_native);
//!
//! // Check that the value of &a * &b is correct.
//! assert_eq!((&a * &b).value()?, a_native * &b_native);
//!
//! // Check that operations on variables and constants are equivalent.
//! (&a + &b).enforce_equal(&(&a_const + &b_const))?;
//! assert!(cs.is_satisfied()?);
//! # Ok(())
//! # }
//! ```
//!
//! One can also perform standard algebraic operations on `EdwardsVar`:
//!
//! ```
//! # fn main() -> Result<(), r1cs_core::SynthesisError> {
//! # use algebra::{UniformRand, ed_on_bls12_377::*};
//! # use r1cs_core::*;
//! # use r1cs_std::prelude::*;
//! # use r1cs_std::ed_on_bls12_377::*;
//!
//! # let cs = ConstraintSystem::<Fq>::new_ref();
//! # let mut rng = algebra::test_rng();
//!
//! // Generate some random `Edwards` elements.
//! let a_native = EdwardsProjective::rand(&mut rng);
//! let b_native = EdwardsProjective::rand(&mut rng);
//!
//! // Allocate `a_native` and `b_native` as witness variables in `cs`.
//! let a = EdwardsVar::new_witness(r1cs_core::ns!(cs, "a"), || Ok(a_native))?;
//! let b = EdwardsVar::new_witness(r1cs_core::ns!(cs, "b"), || Ok(b_native))?;
//!
//! // Allocate `a_native` and `b_native` as constants in `cs`. This does not add any
//! // constraints or variables.
//! let a_const = EdwardsVar::new_constant(r1cs_core::ns!(cs, "a_as_constant"), a_native)?;
//! let b_const = EdwardsVar::new_constant(r1cs_core::ns!(cs, "b_as_constant"), b_native)?;
//!
//! // This returns the identity.
//! let zero = EdwardsVar::zero();
//!
//! // Sanity check one + one = two
//! let two_a = &a + &a + &zero;
//! two_a.enforce_equal(&a.double()?)?;
//!
//! assert!(cs.is_satisfied()?);
//!
//! // Check that the value of &a + &b is correct.
//! assert_eq!((&a + &b).value()?, a_native + &b_native);
//!
//! // Check that operations on variables and constants are equivalent.
//! (&a + &b).enforce_equal(&(&a_const + &b_const))?;
//! assert!(cs.is_satisfied()?);
//! # Ok(())
//! # }
//! ```

mod curves;
mod fields;

pub use curves::*;
pub use fields::*;
