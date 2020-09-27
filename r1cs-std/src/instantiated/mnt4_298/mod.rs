//! This module implements the R1CS equivalent of `algebra::mnt4_298`.
//!
//! It implements field variables for `algebra::mnt4_298::{Fq, Fq2, Fq4}`,
//! group variables for `algebra::mnt4_298::{G1, G2}`, and implements constraint
//! generation for computing `MNT4_298::pairing`.
//!
//! The field underlying these constraints is `algebra::mnt4_298::Fq`.
//!
//! # Examples
//!
//! One can perform standard algebraic operations on `FqVar`:
//!
//! ```
//! # fn main() -> Result<(), r1cs_core::SynthesisError> {
//! use algebra::{UniformRand, mnt4_298::*};
//! use r1cs_core::*;
//! use r1cs_std::prelude::*;
//! use r1cs_std::mnt4_298::*;
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
//! One can also perform standard algebraic operations on `G1Var` and `G2Var`:
//!
//! ```
//! # fn main() -> Result<(), r1cs_core::SynthesisError> {
//! # use algebra::{UniformRand, mnt4_298::*};
//! # use r1cs_core::*;
//! # use r1cs_std::prelude::*;
//! # use r1cs_std::mnt4_298::*;
//!
//! # let cs = ConstraintSystem::<Fq>::new_ref();
//! # let mut rng = algebra::test_rng();
//!
//! // Generate some random `G1` elements.
//! let a_native = G1Projective::rand(&mut rng);
//! let b_native = G1Projective::rand(&mut rng);
//!
//! // Allocate `a_native` and `b_native` as witness variables in `cs`.
//! let a = G1Var::new_witness(r1cs_core::ns!(cs, "a"), || Ok(a_native))?;
//! let b = G1Var::new_witness(r1cs_core::ns!(cs, "b"), || Ok(b_native))?;
//!
//! // Allocate `a_native` and `b_native` as constants in `cs`. This does not add any
//! // constraints or variables.
//! let a_const = G1Var::new_constant(r1cs_core::ns!(cs, "a_as_constant"), a_native)?;
//! let b_const = G1Var::new_constant(r1cs_core::ns!(cs, "b_as_constant"), b_native)?;
//!
//! // This returns the identity of `G1`.
//! let zero = G1Var::zero();
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
//!
//! Finally, one can check pairing computations as well:
//!
//! ```
//! # fn main() -> Result<(), r1cs_core::SynthesisError> {
//! # use algebra::{UniformRand, PairingEngine, mnt4_298::*};
//! # use r1cs_core::*;
//! # use r1cs_std::prelude::*;
//! # use r1cs_std::mnt4_298::{self, *};
//!
//! # let cs = ConstraintSystem::<Fq>::new_ref();
//! # let mut rng = algebra::test_rng();
//!
//! // Generate random `G1` and `G2` elements.
//! let a_native = G1Projective::rand(&mut rng);
//! let b_native = G2Projective::rand(&mut rng);
//!
//! // Allocate `a_native` and `b_native` as witness variables in `cs`.
//! let a = G1Var::new_witness(r1cs_core::ns!(cs, "a"), || Ok(a_native))?;
//! let b = G2Var::new_witness(r1cs_core::ns!(cs, "b"), || Ok(b_native))?;
//!
//! // Allocate `a_native` and `b_native` as constants in `cs`. This does not add any
//! // constraints or variables.
//! let a_const = G1Var::new_constant(r1cs_core::ns!(cs, "a_as_constant"), a_native)?;
//! let b_const = G2Var::new_constant(r1cs_core::ns!(cs, "b_as_constant"), b_native)?;
//!
//! let pairing_result_native = MNT4_298::pairing(a_native, b_native);
//!
//! // Prepare `a` and `b` for pairing.
//! let a_prep = mnt4_298::PairingVar::prepare_g1(&a)?;
//! let b_prep = mnt4_298::PairingVar::prepare_g2(&b)?;
//! let pairing_result = mnt4_298::PairingVar::pairing(a_prep, b_prep)?;
//!
//! // Check that the value of &a + &b is correct.
//! assert_eq!(pairing_result.value()?, pairing_result_native);
//!
//! // Check that operations on variables and constants are equivalent.
//! let a_prep_const = mnt4_298::PairingVar::prepare_g1(&a_const)?;
//! let b_prep_const = mnt4_298::PairingVar::prepare_g2(&b_const)?;
//! let pairing_result_const = mnt4_298::PairingVar::pairing(a_prep_const, b_prep_const)?;
//! println!("Done here 3");
//!
//! pairing_result.enforce_equal(&pairing_result_const)?;
//! assert!(cs.is_satisfied()?);
//! # Ok(())
//! # }
//! ```

mod curves;
mod fields;
mod pairing;

pub use curves::*;
pub use fields::*;
pub use pairing::*;
