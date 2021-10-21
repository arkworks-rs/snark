use crate::UniformRand;
use crate::{
    BitIterator, CanonicalDeserialize, CanonicalSerialize, FromBytesChecked, SemanticallyValid,
};
use std::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Add, AddAssign, Neg, Sub, SubAssign},
};

use crate::{
    bytes::{FromBytes, ToBytes},
    fields::PrimeField,
};
use serde::{Deserialize, Serialize};

#[cfg(test)]
pub mod tests;

pub trait Group:
    ToBytes
    + FromBytes
    + FromBytesChecked
    + SemanticallyValid
    + Serialize
    + for<'a> Deserialize<'a>
    + CanonicalSerialize
    + CanonicalDeserialize
    + Copy
    + Clone
    + Debug
    + Display
    + Default
    + Send
    + Sync
    + 'static
    + Eq
    + Hash
    + Neg<Output = Self>
    + UniformRand
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
{
    type ScalarField: PrimeField + Into<<Self::ScalarField as PrimeField>::BigInt>;

    /// Returns the additive identity.
    fn zero() -> Self;

    /// Returns `self == zero`.
    fn is_zero(&self) -> bool;

    /// Returns `self + self`.
    #[must_use]
    fn double(&self) -> Self;

    /// Sets `self := self + self`.
    fn double_in_place(&mut self) -> &mut Self;

    #[must_use]
    fn mul<'a>(&self, other: &'a Self::ScalarField) -> Self {
        let mut copy = *self;
        copy.mul_assign(other);
        copy
    }

    /// WARNING: This implementation doesn't take costant time with respect
    /// to the exponent, and therefore is susceptible to side-channel attacks.
    /// Be sure to use it in applications where timing (or similar) attacks
    /// are not possible.
    /// TODO: Add a side-channel secure variant.
    fn mul_assign<'a>(&mut self, other: &'a Self::ScalarField) {
        let mut res = Self::zero();
        for i in BitIterator::new(other.into_repr()) {
            res.double_in_place();
            if i {
                res += self
            }
        }
        *self = res
    }
}
