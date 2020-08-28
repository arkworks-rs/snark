use crate::{
    bytes::{FromBytes, ToBytes},
    fields::BitIteratorBE,
    io::{Read, Result as IoResult, Write},
    CanonicalDeserialize, CanonicalSerialize, ConstantSerializedSize, SerializationError,
    UniformRand, Vec,
};
use core::fmt::{Debug, Display};
use rand::{
    distributions::{Distribution, Standard},
    Rng,
};

#[macro_use]
mod macros;

bigint_impl!(BigInteger64, 1);
bigint_impl!(BigInteger128, 2);
bigint_impl!(BigInteger256, 4);
bigint_impl!(BigInteger320, 5);
bigint_impl!(BigInteger384, 6);
bigint_impl!(BigInteger768, 12);
bigint_impl!(BigInteger832, 13);

impl<T: BigInteger> CanonicalSerialize for T {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.write(writer)?;
        Ok(())
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        Self::SERIALIZED_SIZE
    }
}

impl<T: BigInteger> ConstantSerializedSize for T {
    const SERIALIZED_SIZE: usize = Self::NUM_LIMBS * 8;
    const UNCOMPRESSED_SIZE: usize = Self::SERIALIZED_SIZE;
}

impl<T: BigInteger> CanonicalDeserialize for T {
    #[inline]
    fn deserialize<R: Read>(reader: R) -> Result<Self, SerializationError> {
        let value = Self::read(reader)?;
        Ok(value)
    }
}

#[cfg(test)]
mod tests;

/// This defines a `BigInteger`, a smart wrapper around a
/// sequence of `u64` limbs, least-significant limb first.
pub trait BigInteger:
    ToBytes
    + FromBytes
    + CanonicalSerialize
    + ConstantSerializedSize
    + CanonicalDeserialize
    + Copy
    + Clone
    + Debug
    + Default
    + Display
    + Eq
    + Ord
    + Send
    + Sized
    + Sync
    + 'static
    + UniformRand
    + AsMut<[u64]>
    + AsRef<[u64]>
    + From<u64>
{
    /// Number of limbs.
    const NUM_LIMBS: usize;

    /// Add another representation to this one, returning the carry bit.
    fn add_nocarry(&mut self, other: &Self) -> bool;

    /// Subtract another representation from this one, returning the borrow bit.
    fn sub_noborrow(&mut self, other: &Self) -> bool;

    /// Performs a leftwise bitshift of this number, effectively multiplying
    /// it by 2. Overflow is ignored.
    fn mul2(&mut self);

    /// Performs a leftwise bitshift of this number by some amount.
    fn muln(&mut self, amt: u32);

    /// Performs a rightwise bitshift of this number, effectively dividing
    /// it by 2.
    fn div2(&mut self);

    /// Performs a rightwise bitshift of this number by some amount.
    fn divn(&mut self, amt: u32);

    /// Returns true iff this number is odd.
    fn is_odd(&self) -> bool;

    /// Returns true iff this number is even.
    fn is_even(&self) -> bool;

    /// Returns true iff this number is zero.
    fn is_zero(&self) -> bool;

    /// Compute the number of bits needed to encode this number. Always a
    /// multiple of 64.
    fn num_bits(&self) -> u32;

    /// Compute the `i`-th bit of `self`.
    fn get_bit(&self, i: usize) -> bool;

    /// Returns the big integer representation of a given big endian boolean
    /// array.
    fn from_bits(bits: &[bool]) -> Self;

    /// Returns the bit representation in a big endian boolean array, without
    /// leading zeros.
    fn to_bits(&self) -> Vec<bool>;

    /// Returns a vector for wnaf.
    fn find_wnaf(&self) -> Vec<i64>;

    /// Writes this `BigInteger` as a big endian integer. Always writes
    /// `(num_bits` / 8) bytes.
    fn write_le<W: Write>(&self, writer: &mut W) -> IoResult<()> {
        self.write(writer)
    }

    /// Reads a big endian integer occupying (`num_bits` / 8) bytes into this
    /// representation.
    fn read_le<R: Read>(&mut self, reader: &mut R) -> IoResult<()> {
        *self = Self::read(reader)?;
        Ok(())
    }
}

pub mod arithmetic {
    use crate::Vec;
    pub fn find_wnaf(num: &[u64]) -> Vec<i64> {
        let is_zero = |num: &[u64]| num.iter().all(|x| *x == 0u64);
        let is_odd = |num: &[u64]| num[0] & 1 == 1;
        let sub_noborrow = |num: &mut [u64], z: u64| {
            let mut other = vec![0u64; num.len()];
            other[0] = z;
            let mut borrow = 0;

            for (a, b) in num.iter_mut().zip(other) {
                *a = sbb(*a, b, &mut borrow);
            }
        };
        let add_nocarry = |num: &mut [u64], z: u64| {
            let mut other = vec![0u64; num.len()];
            other[0] = z;
            let mut carry = 0;

            for (a, b) in num.iter_mut().zip(other) {
                *a = adc(*a, b, &mut carry);
            }
        };
        let div2 = |num: &mut [u64]| {
            let mut t = 0;
            for i in num.iter_mut().rev() {
                let t2 = *i << 63;
                *i >>= 1;
                *i |= t;
                t = t2;
            }
        };

        let mut num = num.to_vec();
        let mut res = vec![];

        while !is_zero(&num) {
            let z: i64;
            if is_odd(&num) {
                z = 2 - (num[0] % 4) as i64;
                if z >= 0 {
                    sub_noborrow(&mut num, z as u64)
                } else {
                    add_nocarry(&mut num, (-z) as u64)
                }
            } else {
                z = 0;
            }
            res.push(z);
            div2(&mut num);
        }

        res
    }

    /// Calculate a + b + carry, returning the sum and modifying the
    /// carry value.
    #[inline(always)]
    pub(crate) fn adc(a: u64, b: u64, carry: &mut u64) -> u64 {
        let tmp = u128::from(a) + u128::from(b) + u128::from(*carry);

        *carry = (tmp >> 64) as u64;

        tmp as u64
    }

    /// Calculate a - b - borrow, returning the result and modifying
    /// the borrow value.
    #[inline(always)]
    pub(crate) fn sbb(a: u64, b: u64, borrow: &mut u64) -> u64 {
        let tmp = (1u128 << 64) + u128::from(a) - u128::from(b) - u128::from(*borrow);

        *borrow = if tmp >> 64 == 0 { 1 } else { 0 };

        tmp as u64
    }

    /// Calculate a + (b * c) + carry, returning the least significant digit
    /// and setting carry to the most significant digit.
    #[inline(always)]
    pub(crate) fn mac_with_carry(a: u64, b: u64, c: u64, carry: &mut u64) -> u64 {
        let tmp = (u128::from(a)) + u128::from(b) * u128::from(c) + u128::from(*carry);

        *carry = (tmp >> 64) as u64;

        tmp as u64
    }

    #[inline(always)]
    pub(crate) fn mac(a: u64, b: u64, c: u64, carry: &mut u64) -> u64 {
        let tmp = (u128::from(a)) + u128::from(b) * u128::from(c);

        *carry = (tmp >> 64) as u64;

        tmp as u64
    }

    #[inline(always)]
    pub(crate) fn mac_discard(a: u64, b: u64, c: u64, carry: &mut u64) {
        let tmp = (u128::from(a)) + u128::from(b) * u128::from(c);

        *carry = (tmp >> 64) as u64;
    }
}
