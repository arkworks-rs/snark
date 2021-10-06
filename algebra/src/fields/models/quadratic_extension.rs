use std::{
    cmp::{Ord, Ordering, PartialOrd},
    fmt,
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    io::{Read, Write, Result as IoResult},
};
use rand::{
    distributions::{Distribution, Standard},
    Rng,
};

use crate::{
    bytes::{FromBytes, ToBytes},
    bits::{FromBits, ToBits},
    fields::{Field, FpParameters, LegendreSymbol, PrimeField, SquareRootField},
    UniformRand, Error, SemanticallyValid,
    CanonicalSerialize, Flags,
    SerializationError, CanonicalSerializeWithFlags, CanonicalDeserialize,
    CanonicalDeserializeWithFlags, EmptyFlags
};
use crate::biginteger::arithmetic::find_wnaf;
use serde::{Serialize, Deserialize};

/// Model for quadratic extension field of prime field F=Fp
///     F2 = F[X]/(X^2-alpha),
/// with alpha being a (quadratic) "non-residue".
/// We implement the inversion and Karatsuba multiplication according to
/// Mrabet, Joye, Guide to Pairing-based Cryptography
/// https://dl.acm.org/doi/book/10.5555/3092800
/// and the square root algorithm from
/// Adj, et al., Square root computation over even extension fields,
/// https://eprint.iacr.org/2012/685.pdf
pub trait QuadExtParameters: 'static + Send + Sync + Sized {
    /// The prime field that this quadratic extension is eventually an extension of.
    type BasePrimeField: PrimeField;
    /// The base field that this field is a quadratic extension of.
    type BaseField: Field;
    /// The type of the coefficients for an efficient implemntation of the
    /// Frobenius endomorphism.
    type FrobCoeff: Field;

    /// The degree of the extension over the base prime field.
    const DEGREE_OVER_BASE_PRIME_FIELD: usize;

    /// The quadratic non-residue used to construct the extension.
    const NONRESIDUE: Self::BaseField;

    /// Coefficients for the Frobenius automorphism.
    const FROBENIUS_COEFF_C1: &'static [Self::FrobCoeff];

    /// A specializable method for multiplying an element of the base field by
    /// the quadratic non-residue. This is used in Karatsuba multiplication
    /// and in complex squaring.
    #[inline(always)]
    fn mul_base_field_by_nonresidue(fe: &Self::BaseField) -> Self::BaseField {
        Self::NONRESIDUE * fe
    }

    /// A specializable method for multiplying an element of the base field by
    /// the appropriate Frobenius coefficient.
    fn mul_base_field_by_frob_coeff(fe: &mut Self::BaseField, power: usize);

    fn cyclotomic_square(fe: &QuadExtField<Self>) -> QuadExtField<Self> {
        fe.square()
    }
}

#[derive(Derivative)]
#[derivative(
Default(bound = "P: QuadExtParameters"),
Hash(bound = "P: QuadExtParameters"),
Clone(bound = "P: QuadExtParameters"),
Copy(bound = "P: QuadExtParameters"),
Debug(bound = "P: QuadExtParameters"),
PartialEq(bound = "P: QuadExtParameters"),
Eq(bound = "P: QuadExtParameters"),
)]
#[derive(Serialize, Deserialize)]
pub struct QuadExtField<P: QuadExtParameters> {
    pub c0: P::BaseField,
    pub c1: P::BaseField,
    #[derivative(Debug = "ignore")]
    #[serde(skip)]
    #[doc(hidden)]
    pub _parameters: PhantomData<P>,
}

impl<P: QuadExtParameters> QuadExtField<P> {
    pub fn new(c0: P::BaseField, c1: P::BaseField) -> Self {
        QuadExtField {
            c0,
            c1,
            _parameters: PhantomData,
        }
    }

    /// This is only to be used when the element is *known* to be in the cyclotomic subgroup.
    pub fn conjugate(&mut self) {
        self.c1 = -self.c1;
    }

    /// This is only to be used when the element is *known* to be in the cyclotomic subgroup.
    pub fn unitary_inverse(&self) -> Self {
        Self::new(self.c0, -self.c1)
    }


    // (signed) binary square and multiply for r-th roots of unity
    // used for the final exponentiation in the Ate pairing
    pub fn cyclotomic_exp<S: AsRef<[u64]>>(&self, exponent: S) -> Self {
        let mut res = Self::one();
        let self_inverse = self.unitary_inverse();

        let mut found_nonzero = false;
        let naf = find_wnaf(exponent.as_ref());

        for &value in naf.iter().rev() {
            if found_nonzero {
                res = P::cyclotomic_square(&res);
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

    /// Norm of QuadExtField over P::BaseField: Norm(a) = a.x^2 - P::NON_RESIDUE * a.y^2
    pub fn norm(&self) -> P::BaseField {
        let t0 = self.c0.square();
        let mut t1 = self.c1.square();
        t1 = -P::mul_base_field_by_nonresidue(&t1);
        t1.add_assign(&t0);
        t1
    }

    pub fn mul_assign_by_basefield(&mut self, element: &P::BaseField) {
        self.c0.mul_assign(element);
        self.c1.mul_assign(element);
    }
}

impl<P: QuadExtParameters> Field for QuadExtField<P> {
    type BasePrimeField = P::BasePrimeField;

    fn zero() -> Self {
        QuadExtField::new(P::BaseField::zero(), P::BaseField::zero())
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn one() -> Self {
        QuadExtField::new(P::BaseField::one(), P::BaseField::zero())
    }

    fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }

    fn is_odd(&self) -> bool {
        self.c1.is_odd() || ( self.c1.is_zero() && self.c0.is_odd())
    }

    #[inline]
    fn characteristic<'a>() -> &'a [u64] {
        P::BaseField::characteristic()
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
        let v3 = self.c0 - &P::mul_base_field_by_nonresidue(&self.c1);
        // v2 = c0 * c1
        let v2 = self.c0 * &self.c1;

        // v0 = (v0 * v3) + v2
        v0 *= &v3;
        v0 += &v2;

        self.c1 = v2.double();
        self.c0 = v0 + &P::mul_base_field_by_nonresidue(&v2);

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
            v0 -= &P::mul_base_field_by_nonresidue(&v1);
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
        self.c0.frobenius_map(power);
        self.c1.frobenius_map(power);
        P::mul_base_field_by_frob_coeff(&mut self.c1, power);
    }

    #[inline]
    fn from_random_bytes_with_flags<F: Flags>(bytes: &[u8]) -> Option<(Self, F)> {
        let split_at = bytes.len() / 2;
        if let Some(c0) = P::BaseField::from_random_bytes(&bytes[..split_at]) {
            if let Some((c1, flags)) =
            P::BaseField::from_random_bytes_with_flags(&bytes[split_at..])
            {
                return Some((QuadExtField::new(c0, c1), flags));
            }
        }
        None
    }

    #[inline]
    fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
        Self::from_random_bytes_with_flags::<EmptyFlags>(bytes).map(|f| f.0)
    }
}

impl<'a, P: QuadExtParameters> SquareRootField for QuadExtField<P>
    where
        P::BaseField: SquareRootField,
{
    fn legendre(&self) -> LegendreSymbol {
        self.norm().legendre()
    }

    fn sqrt(&self) -> Option<Self> {
        use crate::LegendreSymbol::*;
        if self.c1.is_zero() {
            return self.c0.sqrt().map(|c0| Self::new(c0, P::BaseField::zero()));
        }
        match self.legendre() {
            // Square root based on the complex method. See
            // https://eprint.iacr.org/2012/685.pdf (page 15, algorithm 8)
            Zero => Some(*self),
            QuadraticNonResidue => None,
            QuadraticResidue => {
                let two_inv = P::BaseField::one()
                    .double()
                    .inverse();
                let alpha = self
                    .norm()
                    .sqrt();
                if two_inv.is_none() || alpha.is_none() {
                    return None;
                }
                let mut delta = (alpha.unwrap() + &self.c0) * &two_inv.unwrap();
                if delta.legendre().is_qnr() {
                    delta -= &alpha.unwrap();
                }
                let c0 = delta.sqrt();
                if c0.is_none() {
                    return None;
                }
                let c0_inv = c0.unwrap().inverse();
                if c0_inv.is_none() {
                    return None;
                }
                Some(Self::new(c0.unwrap(), self.c1 * &two_inv.unwrap() * &c0_inv.unwrap()))
            }
        }
    }

    fn sqrt_in_place(&mut self) -> Option<&mut Self> {
        (*self).sqrt().map(|sqrt| {
            *self = sqrt;
            self
        })
    }
}

/// `QuadExtField` elements are ordered lexicographically.
impl<P: QuadExtParameters> Ord for QuadExtField<P> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
        }
    }
}

impl<P: QuadExtParameters> PartialOrd for QuadExtField<P> {
    #[inline(always)]
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<P: QuadExtParameters> From<u128> for QuadExtField<P>
    where
        P::BaseField: From<u128>,
{
    fn from(other: u128) -> Self {
        Self::new(other.into(), P::BaseField::zero())
    }
}

impl<P: QuadExtParameters> From<u64> for QuadExtField<P>
    where
        P::BaseField: From<u64>,
{
    fn from(other: u64) -> Self {
        Self::new(other.into(), P::BaseField::zero())
    }
}

impl<P: QuadExtParameters> From<u32> for QuadExtField<P>
    where
        P::BaseField: From<u32>,
{
    fn from(other: u32) -> Self {
        Self::new(other.into(), P::BaseField::zero())
    }
}

impl<P: QuadExtParameters> From<u16> for QuadExtField<P>
    where
        P::BaseField: From<u16>,
{
    fn from(other: u16) -> Self {
        Self::new(other.into(), P::BaseField::zero())
    }
}

impl<P: QuadExtParameters> From<u8> for QuadExtField<P>
    where
        P::BaseField: From<u8>,
{
    fn from(other: u8) -> Self {
        Self::new(other.into(), P::BaseField::zero())
    }
}

impl<P: QuadExtParameters> ToBytes for QuadExtField<P> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.c0.write(&mut writer)?;
        self.c1.write(writer)
    }
}

impl<P: QuadExtParameters> FromBytes for QuadExtField<P> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let c0 = P::BaseField::read(&mut reader)?;
        let c1 = P::BaseField::read(reader)?;
        Ok(QuadExtField::new(c0, c1))
    }
}

impl<P: QuadExtParameters> ToBits for QuadExtField<P> {
    fn write_bits(&self) -> Vec<bool> {
        let mut bits = self.c0.write_bits();
        bits.extend_from_slice(self.c1.write_bits().as_slice());
        bits

    }
}

impl<P: QuadExtParameters> FromBits for QuadExtField<P> {
    fn read_bits(bits: Vec<bool>) -> Result<Self, Error> {
        let size = (P::DEGREE_OVER_BASE_PRIME_FIELD/2) * <P::BasePrimeField as PrimeField>::Params::MODULUS_BITS as usize;
        let c0 = P::BaseField::read_bits(bits[..size].to_vec())?;
        let c1 = P::BaseField::read_bits(bits[size..].to_vec())?;
        Ok(QuadExtField::new(c0, c1))
    }
}

impl<P: QuadExtParameters> CanonicalSerializeWithFlags for QuadExtField<P> {
    #[inline]
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        mut writer: W,
        flags: F,
    ) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(&self.c0, &mut writer)?;
        self.c1.serialize_with_flags(&mut writer, flags)?;
        Ok(())
    }

    #[inline]
    fn serialized_size_with_flags<F: Flags>(&self) -> usize {
        self.c0.serialized_size() + self.c1.serialized_size_with_flags::<F>()
    }
}

impl<P: QuadExtParameters> CanonicalSerialize for QuadExtField<P> {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.serialize_with_flags(writer, EmptyFlags)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        self.serialized_size_with_flags::<EmptyFlags>()
    }
}

impl<P: QuadExtParameters> CanonicalDeserializeWithFlags for QuadExtField<P> {
    #[inline]
    fn deserialize_with_flags<R: Read, F: Flags>(
        mut reader: R,
    ) -> Result<(Self, F), SerializationError> {
        let c0: P::BaseField = CanonicalDeserialize::deserialize(&mut reader)?;
        let (c1, flags): (P::BaseField, _) =
            CanonicalDeserializeWithFlags::deserialize_with_flags(&mut reader)?;
        Ok((QuadExtField::new(c0, c1), flags))
    }
}

impl<P: QuadExtParameters> CanonicalDeserialize for QuadExtField<P> {
    #[inline]
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let c0: P::BaseField = CanonicalDeserialize::deserialize(&mut reader)?;
        let c1: P::BaseField = CanonicalDeserialize::deserialize(&mut reader)?;
        Ok(QuadExtField::new(c0, c1))
    }
}

impl<P: QuadExtParameters> SemanticallyValid for QuadExtField<P> {
    #[inline]
    fn is_valid(&self) -> bool {
        self.c0.is_valid() && self.c1.is_valid()
    }
}

impl<P: QuadExtParameters> Neg for QuadExtField<P> {
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

impl<P: QuadExtParameters> Distribution<QuadExtField<P>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> QuadExtField<P> {
        QuadExtField::new(UniformRand::rand(rng), UniformRand::rand(rng))
    }
}

impl<'a, P: QuadExtParameters> Add<&'a QuadExtField<P>> for QuadExtField<P> {
    type Output = Self;

    #[inline]
    fn add(self, other: &Self) -> Self {
        let mut result = self;
        result.add_assign(other);
        result
    }
}

impl<'a, P: QuadExtParameters> Sub<&'a QuadExtField<P>> for QuadExtField<P> {
    type Output = Self;

    #[inline]
    fn sub(self, other: &Self) -> Self {
        let mut result = self;
        result.sub_assign(other);
        result
    }
}

impl<'a, P: QuadExtParameters> Mul<&'a QuadExtField<P>> for QuadExtField<P> {
    type Output = Self;

    #[inline]
    fn mul(self, other: &Self) -> Self {
        let mut result = self;
        result.mul_assign(other);
        result
    }
}

impl<'a, P: QuadExtParameters> Div<&'a QuadExtField<P>> for QuadExtField<P> {
    type Output = Self;

    #[inline]
    fn div(self, other: &Self) -> Self {
        let mut result = self;
        result.mul_assign(&other.inverse().unwrap());
        result
    }
}

impl<'a, P: QuadExtParameters> AddAssign<&'a Self> for QuadExtField<P> {
    #[inline]
    fn add_assign(&mut self, other: &Self) {
        self.c0.add_assign(&other.c0);
        self.c1.add_assign(&other.c1);
    }
}

impl<'a, P: QuadExtParameters> SubAssign<&'a Self> for QuadExtField<P> {
    #[inline]
    fn sub_assign(&mut self, other: &Self) {
        self.c0.sub_assign(&other.c0);
        self.c1.sub_assign(&other.c1);
    }
}

impl_additive_ops_from_ref!(QuadExtField, QuadExtParameters);
impl_multiplicative_ops_from_ref!(QuadExtField, QuadExtParameters);

impl<'a, P: QuadExtParameters> MulAssign<&'a Self> for QuadExtField<P> {
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
        self.c0 = v0 + &P::mul_base_field_by_nonresidue(&v1);
    }
}

impl<'a, P: QuadExtParameters> DivAssign<&'a Self> for QuadExtField<P> {
    #[inline]
    fn div_assign(&mut self, other: &Self) {
        self.mul_assign(&other.inverse().unwrap());
    }
}

impl<P: QuadExtParameters> fmt::Display for QuadExtField<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "QuadExtField({} + {} * u)", self.c0, self.c1)
    }
}