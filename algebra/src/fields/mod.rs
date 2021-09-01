use crate::{
    biginteger::BigInteger, bytes::{FromBytes, ToBytes}, UniformRand, bits::{ToBits, FromBits},
    Error, BitSerializationError, SemanticallyValid, FromBytesChecked,
    serialize:: {
        CanonicalSerialize, CanonicalDeserialize,
        CanonicalSerializeWithFlags, CanonicalDeserializeWithFlags,
        Flags, EmptyFlags
    }
};
use std::{
    fmt::{Debug, Display},
    hash::Hash,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    str::FromStr,
};
use serde::{Serialize, Deserialize};

#[macro_use]
mod macros;

#[cfg(feature = "bls12_377")]
pub mod bls12_377;

#[cfg(feature = "bls12_381")]
pub mod bls12_381;

#[cfg(feature = "bn_382")]
pub mod bn_382;

#[cfg(feature = "edwards_bls12")]
pub mod edwards_bls12;

#[cfg(feature = "edwards_sw6")]
pub mod edwards_sw6;

#[cfg(feature = "jubjub")]
pub mod jubjub;

#[cfg(feature = "mnt4_753")]
pub mod mnt4753;

#[cfg(feature = "mnt6_753")]
pub mod mnt6753;

#[cfg(feature = "mnt6")]
pub mod mnt6;

#[cfg(feature = "sw6")]
pub mod sw6;

#[cfg(feature = "tweedle")]
pub mod tweedle;

#[cfg(test)]
pub mod tests;

#[macro_use]
pub mod arithmetic;

pub mod models;
pub use self::models::*;

#[macro_export]
macro_rules! field_new {
    ($name:ident, $c0:expr) => {
        $name {
            0: $c0, 
            1: std::marker::PhantomData
        }
    };
    ($name:ident, $c0:expr, $c1:expr $(,)?) => {
        $name {
            c0: $c0,
            c1: $c1,
            _parameters: std::marker::PhantomData,
        }
    };
    ($name:ident, $c0:expr, $c1:expr, $c2:expr $(,)?) => {
        $name {
            c0: $c0,
            c1: $c1,
            c2: $c2,
            _parameters: std::marker::PhantomData,
        }
    };
}

pub trait MulShort<Rhs = Self> {

    type Output;

    #[must_use]
    fn mul_short(self, rhs: Rhs) -> Self::Output;
}

pub trait MulShortAssign<Rhs = Self> {

    fn mul_short_assign(&mut self, rhs: Rhs);
}

/// The interface for a generic field.
pub trait Field:
    ToBytes
    + FromBytes
    + FromBytesChecked
    + ToBits
    + FromBits
    + Serialize
    + for <'a> Deserialize<'a>
    + CanonicalSerialize
    + CanonicalSerializeWithFlags
    + CanonicalDeserialize
    + CanonicalDeserializeWithFlags
    + SemanticallyValid
    + Copy
    + Clone
    + Debug
    + Display
    + Default
    + Send
    + Sync
    + 'static
    + Eq
    + Ord
    + Neg<Output = Self>
    + UniformRand
    + Sized
    + Hash
    + From<u128>
    + From<u64>
    + From<u32>
    + From<u16>
    + From<u8>
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + AddAssign<Self>
    + SubAssign<Self>
    + MulAssign<Self>
    + DivAssign<Self>
    + for<'a> Add<&'a Self, Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + for<'a> Div<&'a Self, Output = Self>
    + for<'a> AddAssign<&'a Self>
    + for<'a> SubAssign<&'a Self>
    + for<'a> MulAssign<&'a Self>
    + for<'a> DivAssign<&'a Self>
    + std::iter::Sum<Self>
    + for<'a> std::iter::Sum<&'a Self>
    + std::iter::Product<Self>
    + for<'a> std::iter::Product<&'a Self>
{
    type BasePrimeField: PrimeField;

    /// Returns the zero element of the field, the additive identity.
    fn zero() -> Self;

    /// Returns true if and only if `self == Self::zero()`.
    fn is_zero(&self) -> bool;

    /// Returns the one element of the field, a field generator.
    fn one() -> Self;

    /// Returns true if and only if `self == Self::one()`.
    fn is_one(&self) -> bool;

    /// Returns true iff self is odd
    fn is_odd(&self) -> bool;

    /// Returns the characteristic of the field.
    fn characteristic<'a>() -> &'a [u64] {
        Self::BasePrimeField::characteristic()
    }

    /// Returns `self + self`.
    #[must_use]
    fn double(&self) -> Self;

    /// Doubles `self` in place.
    fn double_in_place(&mut self) -> &mut Self;

    /// Returns `self * self`.
    #[must_use]
    fn square(&self) -> Self;

    /// Squares `self` in place.
    fn square_in_place(&mut self) -> &mut Self;

    /// Computes the multiplicative inverse of `self` if `self` is nonzero.
    #[must_use]
    fn inverse(&self) -> Option<Self>;

    // Sets `self` to `self`'s inverse if it exists. Otherwise it is a no-op.
    fn inverse_in_place(&mut self) -> Option<&mut Self>;

    /// Returns a field element if the set of bytes forms a valid field element,
    /// otherwise returns None. This function is primarily intended for sampling
    /// random field elements from a hash-function or RNG output.
    fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
        Self::from_random_bytes_with_flags::<EmptyFlags>(bytes).map(|f| f.0)
    }

    /// Returns a field element with an extra sign bit used for group parsing if
    /// the set of bytes forms a valid field element, otherwise returns
    /// None. This function is primarily intended for sampling
    /// random field elements from a hash-function or RNG output.
    fn from_random_bytes_with_flags<F: Flags>(bytes: &[u8]) -> Option<(Self, F)>;

    /// Exponentiates this element by a power of the base prime modulus via
    /// the Frobenius automorphism.
    fn frobenius_map(&mut self, power: usize);

    /// Exponentiates this element by a number represented with `u64` limbs,
    /// least significant limb first.
    /// WARNING: This implementation doesn't take costant time with respect
    /// to the exponent, and therefore is susceptible to side-channel attacks.
    /// Be sure to use it in applications where timing (or similar) attacks
    /// are not possible.
    /// TODO: Add a side-channel secure variant.
    fn pow<S: AsRef<[u64]>>(&self, exp: S) -> Self {
        let mut res = Self::one();

        let mut found_one = false;

        for i in BitIterator::new(exp) {
            if !found_one {
                if i {
                    found_one = true;
                } else {
                    continue;
                }
            }

            res.square_in_place();

            if i {
                res *= self;
            }
        }
        res
    }
}

use std::io::{ Read, Result as IoResult };
impl<F: Field> FromBytesChecked for F {
    fn read_checked<R: Read>(reader: R) -> IoResult<Self>
    { Self::read(reader) }
}

/// A trait that defines parameters for a prime field.
pub trait FpParameters: 'static + Send + Sync + Sized {
    type BigInt: BigInteger;

    /// The modulus of the field.
    const MODULUS: Self::BigInt;

    /// The number of bits needed to represent the `Self::MODULUS`.
    const MODULUS_BITS: u32;

    /// The number of bits that must be shaved from the beginning of
    /// the representation when randomly sampling.
    const REPR_SHAVE_BITS: u32;

    /// R = 2^256 % Self::MODULUS
    const R: Self::BigInt;

    /// R2 = R^2 % Self::MODULUS
    const R2: Self::BigInt;

    /// INV = -(MODULUS^{-1} mod MODULUS) mod MODULUS
    const INV: u64;

    /// A multiplicative generator that is also a quadratic nonresidue.
    /// `Self::GENERATOR` is an element having multiplicative order
    /// `Self::MODULUS - 1`.
    /// There also does not exist `x` such that `Self::GENERATOR = x^2 %
    /// Self::MODULUS`
    const GENERATOR: Self::BigInt;

    /// The number of bits that can be reliably stored.
    /// (Should equal `SELF::MODULUS_BITS - 1`)
    const CAPACITY: u32;

    /// 2^s * t = MODULUS - 1 with t odd. This is the two-adicity of the prime.
    const TWO_ADICITY: u32;

    /// 2^s root of unity computed by GENERATOR^t
    const ROOT_OF_UNITY: Self::BigInt;

    /// t for 2^s * t = MODULUS - 1
    const T: Self::BigInt;

    /// (t - 1) / 2
    const T_MINUS_ONE_DIV_TWO: Self::BigInt;

    /// (Self::MODULUS - 1) / 2
    const MODULUS_MINUS_ONE_DIV_TWO: Self::BigInt;

    const SMALL_SUBGROUP_DEFINED: bool = false;

    const SMALL_SUBGROUP_BASE: Option<u64> = None;

    const SMALL_SUBGROUP_POWER: Option<u64> = None;

    // generator^((modulus-1) / (2^s * small_subgroup_base^small_subgroup_power))
    const FULL_ROOT_OF_UNITY: Option<Self::BigInt> = None;
}

/// The interface for a prime field.
pub trait PrimeField: Field<BasePrimeField = Self> + FromStr {
    type Params: FpParameters<BigInt = Self::BigInt>;
    type BigInt: BigInteger;

    /// Returns a prime field element from its underlying representation.
    fn from_repr(repr: <Self::Params as FpParameters>::BigInt) -> Self;

    /// Returns the underlying representation of the prime field element.
    fn into_repr(&self) -> Self::BigInt;

    /// Returns a prime field element from its underlying raw representation.
    fn from_repr_raw(repr: Self::BigInt) -> Self;

    /// Returns the underlying raw representation of the prime field element.
    fn into_repr_raw(&self) -> Self::BigInt;

    /// Returns the multiplicative generator of `char()` - 1 order.
    fn multiplicative_generator() -> Self;

    /// Returns the 2^s root of unity.
    fn root_of_unity() -> Self;

    ///Returns the full root of unity
    fn full_root_of_unity() -> Option<Self>;

    /// Return the a QNR^T
    fn qnr_to_t() -> Self {
        Self::root_of_unity()
    }

    /// Returns the field size in bits.
    fn size_in_bits() -> usize {
        Self::Params::MODULUS_BITS as usize
    }

    /// Returns the trace.
    fn trace() -> Self::BigInt {
        Self::Params::T
    }

    /// Returns the trace minus one divided by two.
    fn trace_minus_one_div_two() -> Self::BigInt {
        Self::Params::T_MINUS_ONE_DIV_TWO
    }

    /// Returns the modulus minus one divided by two.
    fn modulus_minus_one_div_two() -> Self::BigInt {
        Self::Params::MODULUS_MINUS_ONE_DIV_TWO
    }

}

impl<F: PrimeField> ToBits for F {
    #[inline]
    fn write_bits(&self) -> Vec<bool> {
        let num_bits = <Self as PrimeField>::Params::MODULUS_BITS;

        let mut field_char = BitIterator::new(Self::characteristic());
        let mut tmp = Vec::with_capacity(num_bits as usize);
        let mut found_one = false;
        for b in BitIterator::new(self.into_repr()) {
            // Skip leading bits
            found_one |= field_char.next().unwrap();
            if !found_one {
                continue;
            }

            tmp.push(b);
        }

        debug_assert_eq!(tmp.len(), num_bits as usize);

        tmp
    }
}

impl<F: PrimeField> FromBits for F {
    #[inline]
    fn read_bits(bits: Vec<bool>) -> Result<Self, Error> {
        let modulus_bits = <Self as PrimeField>::Params::MODULUS_BITS as usize;

        //NOTE: We allow bits having enough leading bits to zero s.t. the length will be <= F::MODULUS_BITS
        let leading_zeros = leading_zeros(bits.as_slice()) as usize;
        let bits = &bits.as_slice()[leading_zeros..];
        match bits.len() <=  modulus_bits {
            true => {
                let read_bigint = <Self as PrimeField>::BigInt::from_bits(bits);
                match read_bigint < F::Params::MODULUS {
                    true => Ok(Self::from_repr(read_bigint)),
                    false => {
                        let e = Box::new(
                            BitSerializationError::InvalidFieldElement("element is over the field modulus".to_owned())
                        );
                        Err(e)
                    }
                }
            },
            false => {
               let e = Box::new(
                   BitSerializationError::InvalidFieldElement(format!("bit vec length is greater than the modulus bits ({})", modulus_bits))
               );
                Err(e)
            }
        }
    }
}

/// Attempts to convert a boolean vec into a valid field element for field `ToF`.
/// If `from` is not a valid element for field ToF, this function returns None.
pub fn convert<ToF: PrimeField>(from: Vec<bool>) -> Result<ToF, Error> {
    ToF::read_bits(from)
}

#[inline]
pub fn leading_zeros(bits: &[bool]) -> u32 {
    let mut ctr = 0;
    for &b in bits.iter() {
        if !b {
            ctr += 1;
        } else {
            break;
        }
    }
    ctr
}

/// The interface for a field that supports an efficient square-root operation.
pub trait SquareRootField: Field {
    /// Returns the Legendre symbol.
    fn legendre(&self) -> LegendreSymbol;

    /// Returns the square root of self, if it exists.
    #[must_use]
    fn sqrt(&self) -> Option<Self>;

    /// Sets `self` to be the square root of `self`, if it exists.
    fn sqrt_in_place(&mut self) -> Option<&mut Self>;
}

#[derive(Debug, PartialEq)]
pub enum LegendreSymbol {
    Zero = 0,
    QuadraticResidue = 1,
    QuadraticNonResidue = -1,
}

impl LegendreSymbol {
    pub fn is_zero(&self) -> bool {
        *self == LegendreSymbol::Zero
    }

    pub fn is_qnr(&self) -> bool {
        *self == LegendreSymbol::QuadraticNonResidue
    }

    pub fn is_qr(&self) -> bool {
        *self == LegendreSymbol::QuadraticResidue
    }
}

#[derive(Debug)]
pub struct BitIterator<E> {
    t: E,
    n: usize,
}

impl<E: AsRef<[u64]>> BitIterator<E> {
    pub fn new(t: E) -> Self {
        let n = t.as_ref().len() * 64;

        BitIterator { t, n }
    }
}

impl<E: AsRef<[u64]>> Iterator for BitIterator<E> {
    type Item = bool;

    fn next(&mut self) -> Option<bool> {
        if self.n == 0 {
            None
        } else {
            self.n -= 1;
            let part = self.n / 64;
            let bit = self.n - (64 * part);

            Some(self.t.as_ref()[part] & (1 << bit) > 0)
        }
    }
}

use crate::biginteger::{
    BigInteger256, BigInteger320, BigInteger384, BigInteger768, BigInteger832,
};

impl_field_bigint_conv!(Fp256, BigInteger256, Fp256Parameters);
impl_field_bigint_conv!(Fp320, BigInteger320, Fp320Parameters);
impl_field_bigint_conv!(Fp384, BigInteger384, Fp384Parameters);
impl_field_bigint_conv!(Fp768, BigInteger768, Fp768Parameters);
impl_field_bigint_conv!(Fp832, BigInteger832, Fp832Parameters);

pub fn batch_inversion<F: Field>(v: &mut [F]) {
    // Montgomery’s Trick and Fast Implementation of Masked AES
    // Genelle, Prouff and Quisquater
    // Section 3.2

    // First pass: compute [a, ab, abc, ...]
    let mut prod = Vec::with_capacity(v.len());
    let mut tmp = F::one();
    for f in v.iter().filter(|f| !f.is_zero()) {
        tmp.mul_assign(f);
        prod.push(tmp);
    }

    // Invert `tmp`.
    tmp = tmp.inverse().unwrap(); // Guaranteed to be nonzero.

    // Second pass: iterate backwards to compute inverses
    for (f, s) in v.iter_mut()
        // Backwards
        .rev()
        // Ignore normalized elements
        .filter(|f| !f.is_zero())
        // Backwards, skip last element, fill in one for last term.
        .zip(prod.into_iter().rev().skip(1).chain(Some(F::one())))
    {
        // tmp := tmp * g.z; g.z := tmp * s = 1/z
        let newtmp = tmp * *f;
        *f = tmp * &s;
        tmp = newtmp;
    }
}
