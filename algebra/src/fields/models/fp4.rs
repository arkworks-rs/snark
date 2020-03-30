use rand::{Rng, distributions::{Standard, Distribution}};
use crate::{UniformRand, ToBits, FromBits, PrimeField, Error, BitSerializationError};

use std::{
    cmp::Ordering,
    io::{Read, Result as IoResult, Write},
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{bytes::{FromBytes, ToBytes}, fields::{Field, Fp2, Fp2Parameters, FpParameters},
            biginteger::BigInteger, ToCompressedBits, FromCompressedBits};
use crate::fields::SquareRootField;

/// Model for quadratic extension field F4 as towered extension
///
//     F4 = F2[Y]/(Y^2-X),
//     F2 = Fp[X]/(X^2-alpha),
///
/// using a "non-residue" alpha mod p such that (X^4-alpha) is irreducible over Fp.
/// Its arithmetics includes pairing-relevant operations such as exponentiation and
/// squaring on the r-th unit roots of F4 (cyclotomic exp. and squ.).

pub trait Fp4Parameters: 'static + Send + Sync {
    type Fp2Params: Fp2Parameters;

    const NONRESIDUE: Fp2<Self::Fp2Params>;

    /// Coefficients for the Frobenius map.
    const FROBENIUS_COEFF_FP4_C1: [<Self::Fp2Params as Fp2Parameters>::Fp; 4];

    #[inline(always)]
    fn mul_fp2_by_nonresidue(fe: &Fp2<Self::Fp2Params>) -> Fp2<Self::Fp2Params> {
        Self::NONRESIDUE * fe
    }
}

#[derive(Derivative)]
#[derivative(
Default(bound = "P: Fp4Parameters"),
Hash(bound = "P: Fp4Parameters"),
Clone(bound = "P: Fp4Parameters"),
Copy(bound = "P: Fp4Parameters"),
Debug(bound = "P: Fp4Parameters"),
PartialEq(bound = "P: Fp4Parameters"),
Eq(bound = "P: Fp4Parameters")
)]
pub struct Fp4<P: Fp4Parameters> {
    pub c0: Fp2<P::Fp2Params>,
    pub c1: Fp2<P::Fp2Params>,
    #[derivative(Debug = "ignore")]
    _parameters: PhantomData<P>,
}

impl<P: Fp4Parameters> Fp4<P> {
    pub fn new(c0: Fp2<P::Fp2Params>, c1: Fp2<P::Fp2Params>) -> Self {
        Fp4 {
            c0,
            c1,
            _parameters: PhantomData,
        }
    }

    pub fn mul_by_nonresidue(value: &Fp2<P::Fp2Params>) -> Fp2<P::Fp2Params> {
        let new_c0 = P::Fp2Params::mul_fp_by_nonresidue(&value.c1);
        let new_c1 = value.c0;
        Fp2::new(new_c0, new_c1)
    }

    pub fn unitary_inverse(&self) -> Self {
        Self::new(self.c0, self.c1.neg())
    }

    pub fn cyclotomic_square(&self) -> Self{
        let a = self.c1.square();
        let b = self.c1 + &self.c0;
        let c = b.square() - &a;
        let d = Self::mul_by_nonresidue(&a);
        let e = c - &d;
        Fp4::new(d.double() + &Fp2::<P::Fp2Params>::one(), e - &Fp2::<P::Fp2Params>::one())
    }

    // (signed) binary square and multiply for r-th roots of unity
    // used for the final exponentiation in the Ate pairing
    pub fn cyclotomic_exp<B: BigInteger>(&self, exponent: &B) -> Self {
        let mut res = Self::one();
        let self_inverse = self.unitary_inverse();

        let mut found_nonzero = false;
        let naf = exponent.find_wnaf();

        for &value in naf.iter().rev() {
            if found_nonzero {
                res = res.cyclotomic_square();
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

    //Mul by an element of the form (c0: [c0, 0] c1: [c2, c3])
    pub fn mul_by_023(self, other: &Self) -> Self
    {
        let v0 =
            {
                let v0_c0 = self.c0.c0 * &other.c0.c0;
                let v0_c1 = self.c0.c1 * &other.c0.c0;
                Fp2::new(v0_c0, v0_c1)
            };
        let v1 = self.c1 * &other.c1;

        let c0 = v0 + &Self::mul_by_nonresidue(&v1);
        let c1 = (self.c0 + &self.c1) * &(other.c0 + &other.c1) - &v0 - &v1;

        Self::new(c0, c1)
    }

}

impl<P: Fp4Parameters> Field for Fp4<P> {
    fn zero() -> Self {
        Fp4 {
            c0:          Fp2::zero(),
            c1:          Fp2::zero(),
            _parameters: PhantomData,
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn one() -> Self {
        Fp4 {
            c0:          Fp2::one(),
            c1:          Fp2::zero(),
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
        Fp2::<P::Fp2Params>::characteristic()
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
        let mut ab = self.c0;
        ab.mul_assign(&self.c1);
        let mut c0c1 = self.c0;
        c0c1.add_assign(&self.c1);
        let mut c0 = self.c1;
        c0 = Self::mul_by_nonresidue(&c0);
        c0.add_assign(&self.c0);
        c0.mul_assign(&c0c1);
        c0.sub_assign(&ab);
        self.c1 = ab;
        self.c1.add_assign(&ab);
        ab = Self::mul_by_nonresidue(&ab);
        c0.sub_assign(&ab);
        self.c0 = c0;
        self
    }

    fn inverse(&self) -> Option<Self> {
        let mut c0s = self.c0;
        c0s.square_in_place();
        let mut c1s = self.c1;
        c1s.square_in_place();
        c1s = Self::mul_by_nonresidue(&c1s);
        c0s.sub_assign(&c1s);

        c0s.inverse().map(|t| {
            let mut tmp = Fp4::new(t, t);
            tmp.c0.mul_assign(&self.c0);
            tmp.c1.mul_assign(&self.c1);
            tmp.c1 = tmp.c1.neg();

            tmp
        })
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
        self.c1.c0.mul_assign(&P::FROBENIUS_COEFF_FP4_C1[power % 4]);
        self.c1.c1.mul_assign(&P::FROBENIUS_COEFF_FP4_C1[power % 4]);
    }
}

/// `Fp4` elements are ordered lexicographically.
impl<P: Fp4Parameters> Ord for Fp4<P> {
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

impl<P: Fp4Parameters> PartialOrd for Fp4<P> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<P: Fp4Parameters> From<u128> for Fp4<P> {
    fn from(other: u128) -> Self {
        Self::new(other.into(), Fp2::zero())
    }
}

impl<P: Fp4Parameters> From<u64> for Fp4<P> {
    fn from(other: u64) -> Self {
        Self::new(other.into(), Fp2::zero())
    }
}

impl<P: Fp4Parameters> From<u32> for Fp4<P> {
    fn from(other: u32) -> Self {
        Self::new(other.into(), Fp2::zero())
    }
}

impl<P: Fp4Parameters> From<u16> for Fp4<P> {
    fn from(other: u16) -> Self {
        Self::new(other.into(), Fp2::zero())
    }
}

impl<P: Fp4Parameters> From<u8> for Fp4<P> {
    fn from(other: u8) -> Self {
        Self::new(other.into(), Fp2::zero())
    }
}

impl<P: Fp4Parameters> ToBytes for Fp4<P> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.c0.write(&mut writer)?;
        self.c1.write(&mut writer)
    }
}

impl<P: Fp4Parameters> FromBytes for Fp4<P> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let c0 = Fp2::read(&mut reader)?;
        let c1 = Fp2::read(&mut reader)?;
        Ok(Fp4::new(c0, c1))
    }
}

impl<P: Fp4Parameters> ToBits for Fp4<P> {
    fn write_bits(&self) -> Vec<bool> {
        let mut bits = self.c0.write_bits();
        bits.extend_from_slice(self.c1.write_bits().as_slice());
        bits

    }
}

impl<P: Fp4Parameters> FromBits for Fp4<P> {
    fn read_bits(bits: Vec<bool>) -> Result<Self, Error> {
        let size = 2 * <<P::Fp2Params as Fp2Parameters>::Fp as PrimeField>::Params::MODULUS_BITS as usize;
        let c0 = Fp2::read_bits(bits[..size].to_vec())?;
        let c1 = Fp2::read_bits(bits[size..].to_vec())?;
        Ok(Fp4::new(c0, c1))
    }
}


/*  Note: compression and decompression of a Fqk element is possible thanks to a property of Ate pairing.
    if c0 + i*c1 is the output of an Ate pairing, then holds that c0^2 - nr * c1^2 = 1.
    Therefore, we can save c1 and compute c0 as sqrt(1 + nr*c1^2), dedicating a bit also for the sign
    of the result.
*/

//Fq4 is the output type of a pairing using MNT4_753 as pairing curve. Therefore we
//can compress/decompress it.

impl<P: Fp4Parameters> ToCompressedBits for Fp4<P> {

    #[inline]
    fn compress(&self) -> Vec<bool> {

        //Serialize c1
        let mut res = self.c1.write_bits();

        //Set the MSB to indicate the parity of c0
        let parity = self.c0.is_odd();
        res.push(parity);

        res
    }
}

impl<P: Fp4Parameters> FromCompressedBits for Fp4<P> {

    #[inline]
    fn decompress(compressed: Vec<bool>) -> Result<Self, Error> {
        let len = compressed.len() - 1;
        let parity_flag_set = compressed[len];

        //Mask away the flag bits and try to get the c1 component
        let c1 = Fp2::read_bits(compressed[..len].to_vec())?;

        //Compute c0
        let c0 = {
            let t = Fp2::one() + &Self::mul_by_nonresidue(&(c1.square()));
            t.sqrt()
        };

        match c0 {

            //Estabilish c0 parity
            Some(c0_u) => {
                let neg_c0u = c0_u.neg();
                let c0_s = if c0_u.is_odd() ^ parity_flag_set {neg_c0u} else {c0_u};
                Ok(Self::new(c0_s, c1))
            },

            //sqrt(1 + nr*c1^2) doesn't exists in the field
            _ => Err(Box::new(BitSerializationError::UndefinedSqrt)),
        }
    }
}

impl<P: Fp4Parameters> Neg for Fp4<P> {
    type Output = Self;
    #[inline]
    fn neg(mut self) -> Self {
        self.c0 = self.c0.neg();
        self.c1 = self.c1.neg();
        self
    }
}

impl<P: Fp4Parameters> Distribution<Fp4<P>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Fp4<P> {
        Fp4::new(UniformRand::rand(rng), UniformRand::rand(rng))
    }
}

impl<'a, P: Fp4Parameters> Add<&'a Fp4<P>> for Fp4<P> {
    type Output = Self;

    #[inline]
    fn add(self, other: &Self) -> Self {
        let mut result = self;
        result.add_assign(&other);
        result
    }
}

impl<'a, P: Fp4Parameters> Sub<&'a Fp4<P>> for Fp4<P> {
    type Output = Self;

    #[inline]
    fn sub(self, other: &Self) -> Self {
        let mut result = self;
        result.sub_assign(&other);
        result
    }
}

impl<'a, P: Fp4Parameters> Mul<&'a Fp4<P>> for Fp4<P> {
    type Output = Self;

    #[inline]
    fn mul(self, other: &Self) -> Self {
        let mut result = self;
        result.mul_assign(&other);
        result
    }
}

impl<'a, P: Fp4Parameters> Div<&'a Fp4<P>> for Fp4<P> {
    type Output = Self;

    #[inline]
    fn div(self, other: &Self) -> Self {
        let mut result = self;
        result.mul_assign(&other.inverse().unwrap());
        result
    }
}

impl<'a, P: Fp4Parameters> AddAssign<&'a Self> for Fp4<P> {
    #[inline]
    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }
}

impl<'a, P: Fp4Parameters> SubAssign<&'a Self> for Fp4<P> {
    #[inline]
    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }
}

impl<'a, P: Fp4Parameters> MulAssign<&'a Self> for Fp4<P> {
    #[inline]
    fn mul_assign(&mut self, other: &Self) {
        let mut aa = self.c0;
        aa.mul_assign(&other.c0);
        let mut bb = self.c1;
        bb.mul_assign(&other.c1);
        let mut o = other.c0;
        o.add_assign(&other.c1);
        self.c1.add_assign(&self.c0);
        self.c1.mul_assign(&o);
        self.c1.sub_assign(&aa);
        self.c1.sub_assign(&bb);
        self.c0 = bb;
        self.c0 = Self::mul_by_nonresidue(&self.c0);
        self.c0.add_assign(&aa);
    }
}

impl<'a, P: Fp4Parameters> DivAssign<&'a Self> for Fp4<P> {
    #[inline]
    fn div_assign(&mut self, other: &Self) {
        self.mul_assign(&other.inverse().unwrap());
    }
}

impl<P: Fp4Parameters> ::std::fmt::Display for Fp4<P> {
    fn fmt(&self, f: &mut ::std::fmt::Formatter<'_>) -> ::std::fmt::Result {
        write!(f, "Fp4({}, {})", self.c0, self.c1)
    }
}
