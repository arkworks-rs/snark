use crate::UniformRand;
use std::{
    cmp::{Ord, Ordering, PartialOrd},
    io::{Read, Result as IoResult, Write},
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use rand::{Rng, distributions::{Standard, Distribution}};

use crate::{
    bytes::{FromBytes, ToBytes},
    fields::{Field, LegendreSymbol, PrimeField, SquareRootField},
};

pub trait Fp2Parameters: 'static + Send + Sync {
    type Fp: PrimeField;

    const NONRESIDUE: Self::Fp;

    const QUADRATIC_NONRESIDUE: (Self::Fp, Self::Fp);

    /// Coefficients for the Frobenius automorphism.
    const FROBENIUS_COEFF_FP2_C1: [Self::Fp; 2];

    #[inline(always)]
    fn mul_fp_by_nonresidue(fe: &Self::Fp) -> Self::Fp {
        Self::NONRESIDUE * fe
    }
}

#[derive(Derivative)]
#[derivative(
    Default(bound = "P: Fp2Parameters"),
    Hash(bound = "P: Fp2Parameters"),
    Clone(bound = "P: Fp2Parameters"),
    Copy(bound = "P: Fp2Parameters"),
    Debug(bound = "P: Fp2Parameters"),
    PartialEq(bound = "P: Fp2Parameters"),
    Eq(bound = "P: Fp2Parameters")
)]
pub struct Fp2<P: Fp2Parameters> {
    pub c0: P::Fp,
    pub c1: P::Fp,
    #[derivative(Debug = "ignore")]
    #[doc(hidden)]
    pub _parameters: PhantomData<P>,
}

impl<P: Fp2Parameters> Fp2<P> {
    pub fn new(c0: P::Fp, c1: P::Fp) -> Self {
        Fp2 {
            c0,
            c1,
            _parameters: PhantomData,
        }
    }

    /// Norm of Fp2 over Fp: Norm(a) = a.x^2 - beta * a.y^2
    pub fn norm(&self) -> P::Fp {
        let t0 = self.c0.square();
        let mut t1 = self.c1.square();
        t1 = -P::mul_fp_by_nonresidue(&t1);
        t1.add_assign(&t0);
        t1
    }

    pub fn mul_by_fp(&mut self, element: &P::Fp) {
        self.c0.mul_assign(&element);
        self.c1.mul_assign(&element);
    }
}

impl<P: Fp2Parameters> Field for Fp2<P> {
    fn zero() -> Self {
        Fp2::new(P::Fp::zero(), P::Fp::zero())
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn one() -> Self {
        Fp2::new(P::Fp::one(), P::Fp::zero())
    }

    fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }

    #[inline]
    fn characteristic<'a>() -> &'a [u64] {
        P::Fp::characteristic()
    }

    fn double(&self) -> Self {
        let mut result = self.clone();
        result.double_in_place();
        result
    }

    fn double_in_place(&mut self) -> &mut Self {
        self.c0.double_in_place();
        self.c1.double_in_place();
        self
    }

    fn square(&self) -> Self {
        let mut result = *self;
        result.square_in_place();
        result
    }

    fn square_in_place(&mut self) -> &mut Self {
        // v0 = c0 - c1
        let mut v0 = self.c0 - &self.c1;
        // v3 = c0 - beta * c1
        let v3 = self.c0 - &P::mul_fp_by_nonresidue(&self.c1);
        // v2 = c0 * c1
        let v2 = self.c0 * &self.c1;

        // v0 = (v0 * v3) + v2
        v0 *= &v3;
        v0 += &v2;

        self.c1 = v2.double();
        self.c0 = v0 + &P::mul_fp_by_nonresidue(&v2);

        self
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // Guide to Pairing-based Cryptography, Algorithm 5.19.
            // v0 = c0.square()
            let mut v0 = self.c0.square();
            // v1 = c1.square()
            let v1 = self.c1.square();
            // v0 = v0 - beta * v1
            v0 -= &P::mul_fp_by_nonresidue(&v1);
            v0.inverse().map(|v1| {
                let c0 = self.c0 * &v1;
                let c1 = -(self.c1 * &v1);
                Self::new(c0, c1)
            })
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

    fn frobenius_map(&mut self, power: usize) {
        self.c1.mul_assign(&P::FROBENIUS_COEFF_FP2_C1[power % 2]);
    }
}

impl<'a, P: Fp2Parameters> SquareRootField for Fp2<P> 
where P::Fp: SquareRootField
{
    fn legendre(&self) -> LegendreSymbol {
        self.norm().legendre()
    }

    fn sqrt(&self) -> Option<Self> {
        use crate::LegendreSymbol::*;
        if self.c1.is_zero() {
            return self.c0.sqrt().map(|c0| Self::new(c0, P::Fp::zero()));
        }
        match self.legendre() {
            // Square root based on the complex method. See
            // https://eprint.iacr.org/2012/685.pdf (page 15, algorithm 8)
            Zero => Some(*self),
            QuadraticNonResidue => None,
            QuadraticResidue => {
                let two_inv = P::Fp::one()
                    .double()
                    .inverse()
                    .expect("Two should always have an inverse");
                let alpha = self
                    .norm()
                    .sqrt()
                    .expect("We are in the QR case, the norm should have a square root");
                let mut delta = (alpha + &self.c0) * &two_inv;
                if delta.legendre().is_qnr() {
                    delta -= &alpha;
                }
                let c0 = delta.sqrt().expect("Delta must have a square root");
                let c0_inv = c0.inverse().expect("c0 must have an inverse");
                Some(Self::new(c0, self.c1 * &two_inv * &c0_inv))
            },
        }
    }

    fn sqrt_in_place(&mut self) -> Option<&mut Self> {
        (*self).sqrt().map(|sqrt| {
            *self = sqrt;
            self
        })
    }
}

/// `Fp2` elements are ordered lexicographically.
impl<P: Fp2Parameters> Ord for Fp2<P> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
        }
    }
}

impl<P: Fp2Parameters> PartialOrd for Fp2<P> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<P: Fp2Parameters> From<u128> for Fp2<P> {
    fn from(other: u128) -> Self {
        Self::new(other.into(), P::Fp::zero())
    }
}

impl<P: Fp2Parameters> From<u64> for Fp2<P> {
    fn from(other: u64) -> Self {
        Self::new(other.into(), P::Fp::zero())
    }
}

impl<P: Fp2Parameters> From<u32> for Fp2<P> {
    fn from(other: u32) -> Self {
        Self::new(other.into(), P::Fp::zero())
    }
}

impl<P: Fp2Parameters> From<u16> for Fp2<P> {
    fn from(other: u16) -> Self {
        Self::new(other.into(), P::Fp::zero())
    }
}

impl<P: Fp2Parameters> From<u8> for Fp2<P> {
    fn from(other: u8) -> Self {
        Self::new(other.into(), P::Fp::zero())
    }
}

impl<P: Fp2Parameters> ToBytes for Fp2<P> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.c0.write(&mut writer)?;
        self.c1.write(writer)
    }
}

impl<P: Fp2Parameters> FromBytes for Fp2<P> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let c0 = P::Fp::read(&mut reader)?;
        let c1 = P::Fp::read(reader)?;
        Ok(Fp2::new(c0, c1))
    }
}

impl<P: Fp2Parameters> Neg for Fp2<P> {
    type Output = Self;
    #[inline]
    #[must_use]
    fn neg(self) -> Self {
        let mut res = self.clone();
        res.c0 = res.c0.neg();
        res.c1 = res.c1.neg();
        res
    }
}

impl<P: Fp2Parameters> Distribution<Fp2<P>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Fp2<P> {
        Fp2::new(UniformRand::rand(rng), UniformRand::rand(rng))
    }
}

impl<'a, P: Fp2Parameters> Add<&'a Fp2<P>> for Fp2<P> {
    type Output = Self;

    #[inline]
    fn add(self, other: &Self) -> Self {
        let mut result = self;
        result.add_assign(&other);
        result
    }
}

impl<'a, P: Fp2Parameters> Sub<&'a Fp2<P>> for Fp2<P> {
    type Output = Self;

    #[inline]
    fn sub(self, other: &Self) -> Self {
        let mut result = self;
        result.sub_assign(&other);
        result
    }
}

impl<'a, P: Fp2Parameters> Mul<&'a Fp2<P>> for Fp2<P> {
    type Output = Self;

    #[inline]
    fn mul(self, other: &Self) -> Self {
        let mut result = self;
        result.mul_assign(&other);
        result
    }
}

impl<'a, P: Fp2Parameters> Div<&'a Fp2<P>> for Fp2<P> {
    type Output = Self;

    #[inline]
    fn div(self, other: &Self) -> Self {
        let mut result = self;
        result.mul_assign(&other.inverse().unwrap());
        result
    }
}

impl<'a, P: Fp2Parameters> AddAssign<&'a Self> for Fp2<P> {
    #[inline]
    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }
}

impl<'a, P: Fp2Parameters> SubAssign<&'a Self> for Fp2<P> {
    #[inline]
    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }
}

impl<'a, P: Fp2Parameters> MulAssign<&'a Self> for Fp2<P> {
    #[inline]
    fn mul_assign(&mut self, other: &Self) {
        // Karatsuba multiplication;
        // Guide to Pairing-based cryprography, Algorithm 5.16.
        let v0 = self.c0 * &other.c0;
        let v1 = self.c1 * &other.c1;

        self.c1 += &self.c0;
        self.c1 *= &(other.c0 + &other.c1);
        self.c1 -= &v0;
        self.c1 -= &v1;
        self.c0 = v0 + &P::mul_fp_by_nonresidue(&v1);
    }
}

impl<'a, P: Fp2Parameters> DivAssign<&'a Self> for Fp2<P> {
    #[inline]
    fn div_assign(&mut self, other: &Self) {
        self.mul_assign(&other.inverse().unwrap());
    }
}

impl<P: Fp2Parameters> std::fmt::Display for Fp2<P> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Fp2({} + {} * u)", self.c0, self.c1)
    }
}
