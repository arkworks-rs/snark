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

use crate::ConstraintVar::*;
use algebra::{Field, PairingEngine};
use snark::{LinearCombination, SynthesisError, Variable};
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub};

pub mod test_constraint_system;
pub mod utils;

pub mod bits;
pub use self::bits::*;

pub mod fields;

pub mod groups;

pub mod pairing;

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

#[derive(Clone, Debug)]
pub enum ConstraintVar<E: PairingEngine> {
    LC(LinearCombination<E>),
    Var(Variable),
}

impl<E: PairingEngine> From<Variable> for ConstraintVar<E> {
    #[inline]
    fn from(var: Variable) -> Self {
        Var(var)
    }
}

impl<E: PairingEngine> From<(E::Fr, Variable)> for ConstraintVar<E> {
    #[inline]
    fn from(coeff_var: (E::Fr, Variable)) -> Self {
        LC(coeff_var.into())
    }
}

impl<E: PairingEngine> From<(E::Fr, LinearCombination<E>)> for ConstraintVar<E> {
    #[inline]
    fn from((coeff, mut lc): (E::Fr, LinearCombination<E>)) -> Self {
        lc *= coeff;
        LC(lc)
    }
}

impl<E: PairingEngine> From<(E::Fr, ConstraintVar<E>)> for ConstraintVar<E> {
    #[inline]
    fn from((coeff, var): (E::Fr, ConstraintVar<E>)) -> Self {
        match var {
            LC(lc) => (coeff, lc).into(),
            Var(var) => (coeff, var).into(),
        }
    }
}

impl<E: PairingEngine> ConstraintVar<E> {
    #[inline]
    pub fn zero() -> Self {
        LC(LinearCombination::zero())
    }

    pub fn negate_in_place(&mut self) {
        match self {
            LC(ref mut lc) => lc.negate_in_place(),
            Var(var) => *self = (-E::Fr::one(), *var).into(),
        }
    }

    pub fn double_in_place(&mut self) {
        match self {
            LC(lc) => lc.double_in_place(),
            Var(var) => *self = (E::Fr::one().double(), *var).into(),
        }
    }
}

impl<E: PairingEngine> Add<LinearCombination<E>> for ConstraintVar<E> {
    type Output = LinearCombination<E>;

    #[inline]
    fn add(self, other_lc: LinearCombination<E>) -> LinearCombination<E> {
        match self {
            LC(lc) => other_lc + lc,
            Var(var) => other_lc + var,
        }
    }
}

impl<E: PairingEngine> Sub<LinearCombination<E>> for ConstraintVar<E> {
    type Output = LinearCombination<E>;

    #[inline]
    fn sub(self, other_lc: LinearCombination<E>) -> LinearCombination<E> {
        let result = match self {
            LC(lc) => other_lc - lc,
            Var(var) => other_lc - var,
        };
        -result
    }
}

impl<E: PairingEngine> Add<LinearCombination<E>> for &ConstraintVar<E> {
    type Output = LinearCombination<E>;

    #[inline]
    fn add(self, other_lc: LinearCombination<E>) -> LinearCombination<E> {
        match self {
            LC(lc) => other_lc + lc,
            Var(var) => other_lc + *var,
        }
    }
}

impl<E: PairingEngine> Sub<LinearCombination<E>> for &ConstraintVar<E> {
    type Output = LinearCombination<E>;

    #[inline]
    fn sub(self, other_lc: LinearCombination<E>) -> LinearCombination<E> {
        let result = match self {
            LC(lc) => other_lc - lc,
            Var(var) => other_lc - *var,
        };
        -result
    }
}

impl<E: PairingEngine> Add<(E::Fr, Variable)> for ConstraintVar<E> {
    type Output = Self;

    #[inline]
    fn add(self, var: (E::Fr, Variable)) -> Self {
        let lc = match self {
            LC(lc) => lc + var,
            Var(var2) => LinearCombination::from(var2) + var,
        };
        LC(lc)
    }
}

impl<E: PairingEngine> AddAssign<(E::Fr, Variable)> for ConstraintVar<E> {
    #[inline]
    fn add_assign(&mut self, var: (E::Fr, Variable)) {
        match self {
            LC(ref mut lc) => *lc += var,
            Var(var2) => *self = LC(LinearCombination::from(*var2) + var),
        };
    }
}

impl<E: PairingEngine> Neg for ConstraintVar<E> {
    type Output = Self;

    #[inline]
    fn neg(mut self) -> Self {
        self.negate_in_place();
        self
    }
}

impl<E: PairingEngine> Mul<E::Fr> for ConstraintVar<E> {
    type Output = Self;

    #[inline]
    fn mul(self, scalar: E::Fr) -> Self {
        match self {
            LC(lc) => LC(lc * scalar),
            Var(var) => (scalar, var).into(),
        }
    }
}

impl<E: PairingEngine> MulAssign<E::Fr> for ConstraintVar<E> {
    #[inline]
    fn mul_assign(&mut self, scalar: E::Fr) {
        match self {
            LC(lc) => *lc *= scalar,
            Var(var) => *self = (scalar, *var).into(),
        }
    }
}

impl<E: PairingEngine> Sub<(E::Fr, Variable)> for ConstraintVar<E> {
    type Output = Self;

    #[inline]
    fn sub(self, (coeff, var): (E::Fr, Variable)) -> Self {
        self + (-coeff, var)
    }
}

impl<E: PairingEngine> Add<Variable> for ConstraintVar<E> {
    type Output = Self;

    fn add(self, other: Variable) -> Self {
        self + (E::Fr::one(), other)
    }
}

impl<E: PairingEngine> Sub<Variable> for ConstraintVar<E> {
    type Output = Self;

    #[inline]
    fn sub(self, other: Variable) -> Self {
        self - (E::Fr::one(), other)
    }
}

impl<'a, E: PairingEngine> Add<&'a Self> for ConstraintVar<E> {
    type Output = Self;

    #[inline]
    fn add(self, other: &'a Self) -> Self {
        let lc = match self {
            LC(lc2) => lc2,
            Var(var) => var.into(),
        };
        let lc2 = match other {
            LC(lc2) => lc + lc2,
            Var(var) => lc + *var,
        };
        LC(lc2)
    }
}

impl<'a, E: PairingEngine> Sub<&'a Self> for ConstraintVar<E> {
    type Output = Self;

    #[inline]
    fn sub(self, other: &'a Self) -> Self {
        let lc = match self {
            LC(lc2) => lc2,
            Var(var) => var.into(),
        };
        let lc2 = match other {
            LC(lc2) => lc - lc2,
            Var(var) => lc - *var,
        };
        LC(lc2)
    }
}

impl<E: PairingEngine> Add<&ConstraintVar<E>> for &ConstraintVar<E> {
    type Output = ConstraintVar<E>;

    #[inline]
    fn add(self, other: &ConstraintVar<E>) -> Self::Output {
        (ConstraintVar::zero() + self) + other
    }
}

impl<E: PairingEngine> Sub<&ConstraintVar<E>> for &ConstraintVar<E> {
    type Output = ConstraintVar<E>;

    #[inline]
    fn sub(self, other: &ConstraintVar<E>) -> Self::Output {
        (ConstraintVar::zero() + self) - other
    }
}

impl<'a, E: PairingEngine> Add<(E::Fr, &'a Self)> for ConstraintVar<E> {
    type Output = Self;

    #[inline]
    fn add(self, (coeff, other): (E::Fr, &'a Self)) -> Self {
        let mut lc = match self {
            LC(lc2) => lc2,
            Var(var) => LinearCombination::zero() + var,
        };

        lc = match other {
            LC(lc2) => lc + (coeff, lc2),
            Var(var) => lc + (coeff, *var),
        };
        LC(lc)
    }
}

impl<'a, E: PairingEngine> Sub<(E::Fr, &'a Self)> for ConstraintVar<E> {
    type Output = Self;

    #[inline]
    fn sub(self, (coeff, other): (E::Fr, &'a Self)) -> Self {
        let mut lc = match self {
            LC(lc2) => lc2,
            Var(var) => LinearCombination::zero() + var,
        };
        lc = match other {
            LC(lc2) => lc - (coeff, lc2),
            Var(var) => lc - (coeff, *var),
        };
        LC(lc)
    }
}
