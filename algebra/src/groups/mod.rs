use crate::{BitIterator, UniformRand};
use num_traits::Zero;
use std::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Add, AddAssign, Neg, Sub, SubAssign},
};

use crate::{
    bytes::{FromBytes, ToBytes},
    fields::PrimeField,
};

#[cfg(test)]
pub mod tests;

pub trait Group:
    ToBytes
    + FromBytes
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
    + Zero
    + Add<Self, Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
{
    type ScalarField: PrimeField + Into<<Self::ScalarField as PrimeField>::BigInt>;

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
