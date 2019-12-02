use rand::{Rng, distributions::{Standard, Distribution}};
use crate::{UniformRand, ToCompressed, FromCompressed};
use std::{
    cmp::Ordering,
    io::{Read, Result as IoResult, Write},
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{
    biginteger::BigInteger,
    bytes::{FromBytes, ToBytes},
    fields::{Field, SquareRootField, Fp3, Fp3Parameters},
    to_bytes,
};

pub trait Fp6Parameters: 'static + Send + Sync {
    type Fp3Params: Fp3Parameters;

    const NONRESIDUE: Fp3<Self::Fp3Params>;

    /// Coefficients for the Frobenius automorphism.
    const FROBENIUS_COEFF_FP6_C1: [<Self::Fp3Params as Fp3Parameters>::Fp; 6];

    #[inline(always)]
    fn mul_fp3_by_nonresidue(fe: &Fp3<Self::Fp3Params>) -> Fp3<Self::Fp3Params> {
        Self::NONRESIDUE * fe
    }
}

#[derive(Derivative)]
#[derivative(
    Default(bound = "P: Fp6Parameters"),
    Hash(bound = "P: Fp6Parameters"),
    Clone(bound = "P: Fp6Parameters"),
    Copy(bound = "P: Fp6Parameters"),
    Debug(bound = "P: Fp6Parameters"),
    PartialEq(bound = "P: Fp6Parameters"),
    Eq(bound = "P: Fp6Parameters")
)]
pub struct Fp6<P: Fp6Parameters> {
    pub c0: Fp3<P::Fp3Params>,
    pub c1: Fp3<P::Fp3Params>,
    #[derivative(Debug = "ignore")]
    #[doc(hidden)]
    pub _parameters: PhantomData<P>,
}

impl<P: Fp6Parameters> Fp6<P> {
    pub fn new(c0: Fp3<P::Fp3Params>, c1: Fp3<P::Fp3Params>) -> Self {
        Fp6 {
            c0,
            c1,
            _parameters: PhantomData,
        }
    }

    /// Multiply by quadratic nonresidue v.
    pub fn mul_by_nonresidue(value: &Fp3<P::Fp3Params>) -> Fp3<P::Fp3Params> {
        let mut res = *value;
        res.c0 = value.c2;
        res.c1 = value.c0;
        res.c2 = value.c1;
        res.c0
            .mul_assign(&<P::Fp3Params as Fp3Parameters>::NONRESIDUE);
        res
    }

    pub fn unitary_inverse(&self) -> Self {
        Self::new(self.c0, -self.c1)
    }

    pub fn cyclotomic_exp<B: BigInteger>(&self, exponent: &B) -> Self {
        let mut res = Self::one();
        let self_inverse = self.unitary_inverse();

        let mut found_nonzero = false;
        let naf = exponent.find_wnaf();

        for &value in naf.iter().rev() {
            if found_nonzero {
                res = res.square();
            }

            if value != 0 {
                found_nonzero = true;

                if value > 0 {
                    res = res * self;
                } else {
                    res = res * &self_inverse;
                }
            }
        }

        res
    }

    //Mul by an element of the form [c0: (0, 0, a), c1: (b, c, d)]
    pub fn mul_by_2345(self, other: &Self) -> Self
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Karatsuba) */
    {
        let v0 = {
            let t = other.c0.c2 * &<P::Fp3Params as Fp3Parameters>::NONRESIDUE;
            Fp3::<P::Fp3Params>::new(self.c0.c1 * &t, self.c0.c2 * &t, self.c0.c0 * &other.c0.c2)
        };
        let v1 = self.c1 * &other.c1;
        let beta_v1 = Self::mul_by_nonresidue(&v1);
        let c0 = v0 + &beta_v1;
        let c1 = (self.c0 + &self.c1) * &(other.c0 + &other.c1) -&v0 -&v1;
        Self::new(c0, c1)
    }
}

impl<P: Fp6Parameters> Field for Fp6<P> {
    fn zero() -> Self {
        Fp6 {
            c0:          Fp3::zero(),
            c1:          Fp3::zero(),
            _parameters: PhantomData,
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn one() -> Self {
        Fp6 {
            c0:          Fp3::one(),
            c1:          Fp3::zero(),
            _parameters: PhantomData,
        }
    }

    fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }

    #[inline]
    fn is_odd(&self) -> bool {
        self.c1.is_odd() || ( self.c1.is_zero() && self.c0.is_odd())
    }

    #[inline]
    fn characteristic<'a>() -> &'a [u64] {
        Fp3::<P::Fp3Params>::characteristic()
    }

    fn double(&self) -> Self {
        let mut result = *self;
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
        // Devegili OhEig Scott Dahab --- Multiplication and Squaring on
        // Pairing-Friendly
        // Fields.pdf; Section 3 (Complex)
        let a = self.c0;
        let mut b = self.c1;
        let ab_add = a + &b;
        let mut ab_mul = a * &b;

        let c0 = ab_add * &(a + &Self::mul_by_nonresidue(&mut b))
            - &ab_mul
            - &Self::mul_by_nonresidue(&mut ab_mul);
        let c1 = ab_mul.double();

        self.c0 = c0;
        self.c1 = c1;
        self
    }

    fn inverse(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // From "High-Speed Software Implementation of the Optimal Ate Pairing over
            // Barreto-Naehrig
            // Curves"; Algorithm 8
            let a = self.c0;
            let b = self.c1;

            let mut t1 = b.square();
            let t0 = a.square() - &Self::mul_by_nonresidue(&mut t1);
            let t2 = t0.inverse().unwrap();

            let c0 = a * &t2;
            let c1 = (b * &t2).neg();

            Some(Self::new(c0, c1))
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
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);
        self.c1
            .mul_assign_by_fp(&P::FROBENIUS_COEFF_FP6_C1[power % 6]);
    }
}

/// `Fp6` elements are ordered lexicographically.
impl<P: Fp6Parameters> Ord for Fp6<P> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        let c1_cmp = self.c1.cmp(&other.c1);
        if c1_cmp == Ordering::Equal {
            self.c0.cmp(&other.c0)
        } else {
            c1_cmp
        }
    }
}

impl<P: Fp6Parameters> PartialOrd for Fp6<P> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<P: Fp6Parameters> From<u128> for Fp6<P> {
    fn from(other: u128) -> Self {
        Self::new(other.into(), Fp3::zero())
    }
}

impl<P: Fp6Parameters> From<u64> for Fp6<P> {
    fn from(other: u64) -> Self {
        Self::new(other.into(), Fp3::zero())
    }
}

impl<P: Fp6Parameters> From<u32> for Fp6<P> {
    fn from(other: u32) -> Self {
        Self::new(other.into(), Fp3::zero())
    }
}

impl<P: Fp6Parameters> From<u16> for Fp6<P> {
    fn from(other: u16) -> Self {
        Self::new(other.into(), Fp3::zero())
    }
}

impl<P: Fp6Parameters> From<u8> for Fp6<P> {
    fn from(other: u8) -> Self {
        Self::new(other.into(), Fp3::zero())
    }
}

impl<P: Fp6Parameters> ToBytes for Fp6<P> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.c0.write(&mut writer)?;
        self.c1.write(&mut writer)
    }
}

impl<P: Fp6Parameters> FromBytes for Fp6<P> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let c0 = Fp3::read(&mut reader)?;
        let c1 = Fp3::read(&mut reader)?;
        Ok(Fp6::new(c0, c1))
    }
}

impl<P: Fp6Parameters> Neg for Fp6<P> {
    type Output = Self;
    #[inline]
    fn neg(mut self) -> Self {
        self.c0 = self.c0.neg();
        self.c1 = self.c1.neg();
        self
    }
}

impl<P: Fp6Parameters> Distribution<Fp6<P>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Fp6<P> {
        Fp6::new(UniformRand::rand(rng), UniformRand::rand(rng))
    }
}


impl<'a, P: Fp6Parameters> Add<&'a Fp6<P>> for Fp6<P> {
    type Output = Self;

    #[inline]
    fn add(self, other: &Self) -> Self {
        let mut result = self;
        result.add_assign(&other);
        result
    }
}

impl<'a, P: Fp6Parameters> Sub<&'a Fp6<P>> for Fp6<P> {
    type Output = Self;

    #[inline]
    fn sub(self, other: &Self) -> Self {
        let mut result = self;
        result.sub_assign(&other);
        result
    }
}

impl<'a, P: Fp6Parameters> Mul<&'a Fp6<P>> for Fp6<P> {
    type Output = Self;

    #[inline]
    fn mul(self, other: &Self) -> Self {
        let mut result = self;
        result.mul_assign(&other);
        result
    }
}

impl<'a, P: Fp6Parameters> Div<&'a Fp6<P>> for Fp6<P> {
    type Output = Self;

    #[inline]
    fn div(self, other: &Self) -> Self {
        let mut result = self;
        result.mul_assign(&other.inverse().unwrap());
        result
    }
}

impl<'a, P: Fp6Parameters> AddAssign<&'a Self> for Fp6<P> {
    #[inline]
    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }
}

impl<'a, P: Fp6Parameters> SubAssign<&'a Self> for Fp6<P> {
    #[inline]
    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }
}

impl<'a, P: Fp6Parameters> MulAssign<&'a Self> for Fp6<P> {
    #[inline]
    fn mul_assign(&mut self, other: &Self) {
        // Devegili OhEig Scott Dahab --- Multiplication and Squaring on
        // Pairing-Friendly
        // Fields.pdf; Section 3 (Karatsuba)
        let a0 = self.c0;
        let b0 = self.c1;
        let a1 = other.c0;
        let b1 = other.c1;

        let a0a1 = a0 * &a1;
        let mut b0b1 = b0 * &b1;
        let beta_b0b1 = Self::mul_by_nonresidue(&mut b0b1);

        let c0 = a0a1 + &beta_b0b1;
        let c1 = (a0 + &b0) * &(a1 + &b1) - &a0a1 - &b0b1;

        self.c0 = c0;
        self.c1 = c1;
    }
}

impl<'a, P: Fp6Parameters> DivAssign<&'a Self> for Fp6<P> {
    #[inline]
    fn div_assign(&mut self, other: &Self) {
        self.mul_assign(&other.inverse().unwrap());
    }
}

impl<'a, P: Fp6Parameters> From<&'a [bool]> for Fp6<P> {
    fn from(_bits: &[bool]) -> Self {
        unimplemented!()
    }
}

impl<P: Fp6Parameters> ::std::fmt::Display for Fp6<P> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
        write!(f, "Fp6_2over3({}, {})", self.c0, self.c1)
    }
}

/*  Note: compression and decompression of a Fqk element is possible thanks to a property of Ate pairing.
    if c0 + i*c1 is the output of an Ate pairing, then holds that c0^2 - nr * c1^2 = 1.
    Therefore, we can save c1 and compute c0 as sqrt(1 + nr*c1^2), dedicating a bit also for the sign
    of the result.
*/

//Fq6_2_over_3 is the output type of a pairing using MNT6_753 as pairing curve. Therefore we
//can compress/decompress it.

impl<P: Fp6Parameters> ToCompressed for Fp6<P> {

    #[inline]
    fn compress(&self) -> Vec<u8> {

        //Serialize c1
        let mut res = to_bytes!(self.c1).unwrap();
        let len = res.len() - 1;

        //Set the MSB to indicate the parity of c0
        let parity = if self.c0.is_odd() {1u8 << 7} else {0u8};
        res[len] |= parity;
        res
    }
}

impl<P: Fp6Parameters> FromCompressed for Fp6<P> {

    #[inline]
    fn decompress(compressed: Vec<u8>) -> Option<Self> {
        let len = compressed.len() - 1;
        let parity_flag_set = bool::read([(compressed[len] >> 7) & 1].as_ref()).unwrap();

        //Mask away the flag bits and try to get the c1 component
        let val = {
            let mut tmp = compressed;
            tmp[len] &= 0b0111_1111;
            Fp3::read(tmp.as_slice())
        };

        match val {
            Ok(c1) => {

                //Compute c0
                let c0 = {
                    let t = Fp3::one() + &Self::mul_by_nonresidue(&(c1.square()));
                    t.sqrt()
                };

                match c0 {

                    //Estabilish c0 sign
                    Some(c0_u) => {
                        let neg_c0u = c0_u.neg();
                        let c0_s = if c0_u.is_odd() ^ parity_flag_set {neg_c0u} else {c0_u};
                        Some(Self::new(c0_s, c1))
                    },

                    //sqrt(1 + nr*c1^2) doesn't exists in the field
                    _ => None,
                }
            },

            //Unable to deserialize c1
            _ => None
        }
    }
}