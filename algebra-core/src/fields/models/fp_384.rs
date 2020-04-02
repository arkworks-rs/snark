use core::{
    cmp::{Ord, Ordering, PartialOrd},
    fmt::{Display, Formatter, Result as FmtResult},
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    str::FromStr,
};
use num_traits::{One, Zero};
use unroll::unroll_for_loops;

use crate::{
    biginteger::{arithmetic as fa, BigInteger as _BigInteger, BigInteger384 as BigInteger},
    bytes::{FromBytes, ToBytes},
    fields::{Field, FpParameters, LegendreSymbol, PrimeField, SquareRootField},
    io::{Read, Result as IoResult, Write},
};

pub trait Fp384Parameters: FpParameters<BigInt = BigInteger> {}

#[derive(Derivative)]
#[derivative(
    Default(bound = ""),
    Hash(bound = ""),
    Clone(bound = ""),
    Copy(bound = ""),
    Debug(bound = ""),
    PartialEq(bound = ""),
    Eq(bound = "")
)]
pub struct Fp384<P>(
    pub BigInteger,
    #[derivative(Debug = "ignore")]
    #[doc(hidden)]
    pub PhantomData<P>,
);

impl<P> Fp384<P> {
    #[inline]
    pub const fn new(element: BigInteger) -> Self {
        Self(element, PhantomData)
    }
}

impl<P: Fp384Parameters> Fp384<P> {
    #[inline]
    pub(crate) fn is_valid(&self) -> bool {
        self.0 < P::MODULUS
    }

    #[inline]
    fn reduce(&mut self) {
        if !self.is_valid() {
            self.0.sub_noborrow(&P::MODULUS);
        }
    }
}

impl<P: Fp384Parameters> Zero for Fp384<P> {
    #[inline]
    fn zero() -> Self {
        Fp384::<P>(BigInteger::from(0), PhantomData)
    }

    #[inline]
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<P: Fp384Parameters> One for Fp384<P> {
    #[inline]
    fn one() -> Self {
        Fp384::<P>(P::R, PhantomData)
    }

    #[inline]
    fn is_one(&self) -> bool {
        self.0 == P::R
    }
}

impl<P: Fp384Parameters> Field for Fp384<P> {
    #[inline]
    fn double(&self) -> Self {
        let mut temp = *self;
        temp.double_in_place();
        temp
    }

    #[inline]
    fn double_in_place(&mut self) -> &mut Self {
        // This cannot exceed the backing capacity.
        self.0.mul2();
        // However, it may need to be reduced.
        self.reduce();
        self
    }

    #[inline]
    fn characteristic<'a>() -> &'a [u64] {
        P::MODULUS.as_ref()
    }

    #[inline]
    fn square(&self) -> Self {
        let mut temp = self.clone();
        temp.square_in_place();
        temp
    }

    impl_field_square_in_place!(6);

    #[inline]
    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // Guajardo Kumar Paar Pelzl
            // Efficient Software-Implementation of Finite Fields with Applications to
            // Cryptography
            // Algorithm 16 (BEA for Inversion in Fp)

            let one = BigInteger::from(1);

            let mut u = self.0;
            let mut v = P::MODULUS;
            let mut b = Fp384::<P>(P::R2, PhantomData); // Avoids unnecessary reduction step.
            let mut c = Self::zero();

            while u != one && v != one {
                while u.is_even() {
                    u.div2();

                    if b.0.is_even() {
                        b.0.div2();
                    } else {
                        b.0.add_nocarry(&P::MODULUS);
                        b.0.div2();
                    }
                }

                while v.is_even() {
                    v.div2();

                    if c.0.is_even() {
                        c.0.div2();
                    } else {
                        c.0.add_nocarry(&P::MODULUS);
                        c.0.div2();
                    }
                }

                if v < u {
                    u.sub_noborrow(&v);
                    b.sub_assign(&c);
                } else {
                    v.sub_noborrow(&u);
                    c.sub_assign(&b);
                }
            }

            if u == one {
                Some(b)
            } else {
                Some(c)
            }
        }
    }

    fn inverse_in_place(&mut self) -> Option<&mut Self> {
        if let Some(inverse) = self.inverse() {
            *self = inverse;
            Some(self)
        } else {
            None
        }
    }

    #[inline]
    fn frobenius_map(&mut self, _: usize) {
        // No-op: No effect in a prime field.
    }
}

impl<P: Fp384Parameters> PrimeField for Fp384<P> {
    type Params = P;
    type BigInt = BigInteger;

    #[inline]
    fn from_repr(r: BigInteger) -> Self {
        let mut r = Fp384(r, PhantomData);
        if r.is_valid() {
            r.mul_assign(&Fp384(P::R2, PhantomData));
            r
        } else {
            Self::zero()
        }
    }

    impl_field_into_repr!(6);

    #[inline]
    fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
        let mut result_bytes = vec![0u8; (Self::zero().0).0.len() * 8];
        for (result_byte, in_byte) in result_bytes.iter_mut().zip(bytes.iter()) {
            *result_byte = *in_byte;
        }
        BigInteger::read(result_bytes.as_slice())
            .ok()
            .and_then(|mut res| {
                res.as_mut()[5] &= 0xffffffffffffffff >> P::REPR_SHAVE_BITS;
                let result = Self::new(res);
                if result.is_valid() {
                    Some(result)
                } else {
                    None
                }
            })
    }

    #[inline]
    fn multiplicative_generator() -> Self {
        Fp384::<P>(P::GENERATOR, PhantomData)
    }

    #[inline]
    fn root_of_unity() -> Self {
        Fp384::<P>(P::ROOT_OF_UNITY, PhantomData)
    }
}

impl<P: Fp384Parameters> SquareRootField for Fp384<P> {
    #[inline]
    fn legendre(&self) -> LegendreSymbol {
        use crate::fields::LegendreSymbol::*;

        // s = self^((MODULUS - 1) // 2)
        let s = self.pow(P::MODULUS_MINUS_ONE_DIV_TWO);
        if s.is_zero() {
            Zero
        } else if s.is_one() {
            QuadraticResidue
        } else {
            QuadraticNonResidue
        }
    }

    #[inline]
    fn sqrt(&self) -> Option<Self> {
        sqrt_impl!(Self, P, self)
    }

    fn sqrt_in_place(&mut self) -> Option<&mut Self> {
        (*self).sqrt().map(|sqrt| {
            *self = sqrt;
            self
        })
    }
}

impl<P: Fp384Parameters> Ord for Fp384<P> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        self.into_repr().cmp(&other.into_repr())
    }
}

impl<P: Fp384Parameters> PartialOrd for Fp384<P> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl_prime_field_from_int!(Fp384, u128, Fp384Parameters);
impl_prime_field_from_int!(Fp384, u64, Fp384Parameters);
impl_prime_field_from_int!(Fp384, u32, Fp384Parameters);
impl_prime_field_from_int!(Fp384, u16, Fp384Parameters);
impl_prime_field_from_int!(Fp384, u8, Fp384Parameters);

impl_prime_field_standard_sample!(Fp384, Fp384Parameters);

impl<P: Fp384Parameters> ToBytes for Fp384<P> {
    #[inline]
    fn write<W: Write>(&self, writer: W) -> IoResult<()> {
        self.into_repr().write(writer)
    }
}

impl<P: Fp384Parameters> FromBytes for Fp384<P> {
    #[inline]
    fn read<R: Read>(reader: R) -> IoResult<Self> {
        BigInteger::read(reader).map(Fp384::from_repr)
    }
}

impl<P: Fp384Parameters> FromStr for Fp384<P> {
    type Err = ();

    /// Interpret a string of numbers as a (congruent) prime field element.
    /// Does not accept unnecessary leading zeroes or a blank string.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(());
        }

        if s == "0" {
            return Ok(Self::zero());
        }

        let mut res = Self::zero();

        let ten = Self::from_repr(<Self as PrimeField>::BigInt::from(10));

        let mut first_digit = true;

        for c in s.chars() {
            match c.to_digit(10) {
                Some(c) => {
                    if first_digit {
                        if c == 0 {
                            return Err(());
                        }

                        first_digit = false;
                    }

                    res.mul_assign(&ten);
                    res.add_assign(&Self::from_repr(<Self as PrimeField>::BigInt::from(
                        u64::from(c),
                    )));
                },
                None => {
                    return Err(());
                },
            }
        }
        if !res.is_valid() {
            Err(())
        } else {
            Ok(res)
        }
    }
}

impl<P: Fp384Parameters> Display for Fp384<P> {
    #[inline]
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "Fp384({})", self.into_repr())
    }
}

impl<P: Fp384Parameters> Neg for Fp384<P> {
    type Output = Self;
    #[inline]
    #[must_use]
    fn neg(self) -> Self {
        if !self.is_zero() {
            let mut tmp = P::MODULUS.clone();
            tmp.sub_noborrow(&self.0);
            Fp384::<P>(tmp, PhantomData)
        } else {
            self
        }
    }
}

impl<'a, P: Fp384Parameters> Add<&'a Fp384<P>> for Fp384<P> {
    type Output = Self;

    #[inline]
    fn add(self, other: &Self) -> Self {
        let mut result = self.clone();
        result.add_assign(other);
        result
    }
}

impl<'a, P: Fp384Parameters> Sub<&'a Fp384<P>> for Fp384<P> {
    type Output = Self;

    #[inline]
    fn sub(self, other: &Self) -> Self {
        let mut result = self.clone();
        result.sub_assign(other);
        result
    }
}

impl<'a, P: Fp384Parameters> Mul<&'a Fp384<P>> for Fp384<P> {
    type Output = Self;

    #[inline]
    fn mul(self, other: &Self) -> Self {
        let mut result = self.clone();
        result.mul_assign(other);
        result
    }
}

impl<'a, P: Fp384Parameters> Div<&'a Fp384<P>> for Fp384<P> {
    type Output = Self;

    #[inline]
    fn div(self, other: &Self) -> Self {
        let mut result = self.clone();
        result.mul_assign(&other.inverse().unwrap());
        result
    }
}

impl_additive_ops_from_ref!(Fp384, Fp384Parameters);
impl_multiplicative_ops_from_ref!(Fp384, Fp384Parameters);

impl<'a, P: Fp384Parameters> AddAssign<&'a Self> for Fp384<P> {
    #[inline]
    fn add_assign(&mut self, other: &Self) {
        // This cannot exceed the backing capacity.
        self.0.add_nocarry(&other.0);
        // However, it may need to be reduced
        self.reduce();
    }
}

impl<'a, P: Fp384Parameters> SubAssign<&'a Self> for Fp384<P> {
    #[inline]
    fn sub_assign(&mut self, other: &Self) {
        // If `other` is larger than `self`, add the modulus to self first.
        if other.0 > self.0 {
            self.0.add_nocarry(&P::MODULUS);
        }

        self.0.sub_noborrow(&other.0);
    }
}

impl<'a, P: Fp384Parameters> MulAssign<&'a Self> for Fp384<P> {
    impl_field_mul_assign!(6);
}

impl<'a, P: Fp384Parameters> DivAssign<&'a Self> for Fp384<P> {
    #[inline]
    fn div_assign(&mut self, other: &Self) {
        self.mul_assign(&other.inverse().unwrap());
    }
}
