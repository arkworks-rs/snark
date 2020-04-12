use crate::{
    io::{Read, Result as IoResult, Write},
    BigInteger, CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags, ConstantSerializedSize, EmptyFlags, Flags, SerializationError,
    UniformRand,
};
use core::{
    cmp::{Ord, Ordering, PartialOrd},
    fmt,
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};
use num_traits::{One, Zero};
use rand::{
    distributions::{Distribution, Standard},
    Rng,
};

use crate::{
    bytes::{FromBytes, ToBytes},
    fields::{Field, Fp2, Fp2Parameters},
};

pub trait Fp4Parameters: 'static + Send + Sync {
    type Fp2Params: Fp2Parameters;

    const NONRESIDUE: <Self::Fp2Params as Fp2Parameters>::Fp;

    /// Coefficients for the Frobenius automorphism.
    /// non_residue^((modulus^i-1)/4) for i=0,1,2,3
    const FROBENIUS_COEFF_FP4_C1: [<Self::Fp2Params as Fp2Parameters>::Fp; 4];

    #[inline(always)]
    fn mul_fp2_by_nonresidue(fe: &Fp2<Self::Fp2Params>) -> Fp2<Self::Fp2Params> {
        Fp2::new(Self::NONRESIDUE * &fe.c1, fe.c0)
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
    #[doc(hidden)]
    pub _parameters: PhantomData<P>,
}

impl<P: Fp4Parameters> Fp4<P> {
    pub fn new(c0: Fp2<P::Fp2Params>, c1: Fp2<P::Fp2Params>) -> Self {
        Fp4 {
            c0,
            c1,
            _parameters: PhantomData,
        }
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
                    res *= self;
                } else {
                    res *= &self_inverse;
                }
            }
        }

        res
    }

    pub fn mul_by_fp(&mut self, element: &<P::Fp2Params as Fp2Parameters>::Fp) {
        self.c0.mul_assign_by_fp(element);
        self.c1.mul_assign_by_fp(element);
    }

    pub fn mul_by_fp2(&mut self, element: &Fp2<P::Fp2Params>) {
        self.c0 *= element;
        self.c1 *= element;
    }
}

impl<P: Fp4Parameters> Zero for Fp4<P> {
    fn zero() -> Self {
        Fp4::new(Fp2::zero(), Fp2::zero())
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }
}

impl<P: Fp4Parameters> One for Fp4<P> {
    fn one() -> Self {
        Fp4::new(Fp2::one(), Fp2::zero())
    }

    fn is_one(&self) -> bool {
        self.c0.is_one() && self.c1.is_zero()
    }
}

impl<P: Fp4Parameters> Field for Fp4<P> {
    #[inline]
    fn characteristic<'a>() -> &'a [u64] {
        <P::Fp2Params as Fp2Parameters>::Fp::characteristic()
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

    #[inline]
    fn from_random_bytes_with_flags(bytes: &[u8]) -> Option<(Self, u8)> {
        let split_at = bytes.len() / 2;
        if let Some(c0) = Fp2::<P::Fp2Params>::from_random_bytes(&bytes[..split_at]) {
            if let Some((c1, flags)) =
                Fp2::<P::Fp2Params>::from_random_bytes_with_flags(&bytes[split_at..])
            {
                return Some((Fp4::new(c0, c1), flags));
            }
        }
        None
    }

    #[inline]
    fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
        Self::from_random_bytes_with_flags(bytes).map(|f| f.0)
    }

    fn square_in_place(&mut self) -> &mut Self {
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        // Devegili, OhEigeartaigh, Scott, Dahab
        // Section 5 (Complex squaring)

        // v1 = c0 + c1
        let mut v1 = self.c0 + &self.c1;
        // v2 = c0 + beta * c1
        let v2 = self.c0 + &P::mul_fp2_by_nonresidue(&self.c1);
        // v0 = c0 * c1
        let v0 = self.c0 * &self.c1;

        // v1 = (v1 * v2) - v0
        v1 *= &v2;
        v1 -= &v0;

        self.c1 = v0.double();
        self.c0 = v1 - &P::mul_fp2_by_nonresidue(&v0);

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
            v0 -= &P::mul_fp2_by_nonresidue(&v1);
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
        self.c1
            .mul_assign_by_fp(&P::FROBENIUS_COEFF_FP4_C1[power % 4]);
    }
}

/// `Fp4` elements are ordered lexicographically.
impl<P: Fp4Parameters> Ord for Fp4<P> {
    #[inline(always)]
    fn cmp(&self, other: &Self) -> Ordering {
        match self.c1.cmp(&other.c1) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self.c0.cmp(&other.c0),
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
        Self::new(other.into(), <Fp2<P::Fp2Params>>::zero())
    }
}

impl<P: Fp4Parameters> From<u64> for Fp4<P> {
    fn from(other: u64) -> Self {
        Self::new(other.into(), <Fp2<P::Fp2Params>>::zero())
    }
}

impl<P: Fp4Parameters> From<u32> for Fp4<P> {
    fn from(other: u32) -> Self {
        Self::new(other.into(), <Fp2<P::Fp2Params>>::zero())
    }
}

impl<P: Fp4Parameters> From<u16> for Fp4<P> {
    fn from(other: u16) -> Self {
        Self::new(other.into(), <Fp2<P::Fp2Params>>::zero())
    }
}

impl<P: Fp4Parameters> From<u8> for Fp4<P> {
    fn from(other: u8) -> Self {
        Self::new(other.into(), <Fp2<P::Fp2Params>>::zero())
    }
}

impl<P: Fp4Parameters> ToBytes for Fp4<P> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.c0.write(&mut writer)?;
        self.c1.write(writer)
    }
}

impl<P: Fp4Parameters> FromBytes for Fp4<P> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let c0 = <Fp2<P::Fp2Params>>::read(&mut reader)?;
        let c1 = <Fp2<P::Fp2Params>>::read(reader)?;
        Ok(Fp4::new(c0, c1))
    }
}

impl<P: Fp4Parameters> Neg for Fp4<P> {
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
        result.add_assign(other);
        result
    }
}

impl<'a, P: Fp4Parameters> Sub<&'a Fp4<P>> for Fp4<P> {
    type Output = Self;

    #[inline]
    fn sub(self, other: &Self) -> Self {
        let mut result = self;
        result.sub_assign(other);
        result
    }
}

impl<'a, P: Fp4Parameters> Mul<&'a Fp4<P>> for Fp4<P> {
    type Output = Self;

    #[inline]
    fn mul(self, other: &Self) -> Self {
        let mut result = self;
        result.mul_assign(other);
        result
    }
}

impl<'a, P: Fp4Parameters> Div<&'a Fp4<P>> for Fp4<P> {
    type Output = Self;

    #[inline]
    fn div(self, other: &Self) -> Self {
        let mut result = self;
        result.div_assign(other);
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

impl_additive_ops_from_ref!(Fp4, Fp4Parameters);
impl_multiplicative_ops_from_ref!(Fp4, Fp4Parameters);

impl<'a, P: Fp4Parameters> MulAssign<&'a Self> for Fp4<P> {
    #[inline]
    fn mul_assign(&mut self, other: &Self) {
        // Karatsuba multiplication;
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        // Devegili, OhEigeartaigh, Scott, Dahab
        let v0 = self.c0 * &other.c0;
        let v1 = self.c1 * &other.c1;

        self.c1 += &self.c0;
        self.c1 *= &(other.c0 + &other.c1);
        self.c1 -= &v0;
        self.c1 -= &v1;
        self.c0 = v0 + &P::mul_fp2_by_nonresidue(&v1);
    }
}

impl<'a, P: Fp4Parameters> DivAssign<&'a Self> for Fp4<P> {
    #[inline]
    fn div_assign(&mut self, other: &Self) {
        self.mul_assign(&other.inverse().unwrap());
    }
}

impl<P: Fp4Parameters> fmt::Display for Fp4<P> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Fp4({} + {} * u)", self.c0, self.c1)
    }
}

impl<P: Fp4Parameters> CanonicalSerializeWithFlags for Fp4<P> {
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        writer: &mut W,
        flags: F,
    ) -> Result<(), SerializationError> {
        self.c0.serialize(writer)?;
        self.c1.serialize_with_flags(writer, flags)?;
        Ok(())
    }
}

impl<P: Fp4Parameters> CanonicalSerialize for Fp4<P> {
    fn serialize<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        self.serialize_with_flags(writer, EmptyFlags)
    }

    fn serialized_size(&self) -> usize {
        self.c0.serialized_size() + self.c1.serialized_size()
    }
}

impl<P: Fp4Parameters> ConstantSerializedSize for Fp4<P> {
    const SERIALIZED_SIZE: usize =
        2 * <Fp2<P::Fp2Params> as ConstantSerializedSize>::SERIALIZED_SIZE;
    const UNCOMPRESSED_SIZE: usize = Self::SERIALIZED_SIZE;
}

impl<P: Fp4Parameters> CanonicalDeserializeWithFlags for Fp4<P> {
    fn deserialize_with_flags<R: Read, F: Flags>(
        reader: &mut R,
    ) -> Result<(Self, F), SerializationError> {
        let c0: Fp2<P::Fp2Params> = CanonicalDeserialize::deserialize(reader)?;
        let (c1, flags): (Fp2<P::Fp2Params>, _) =
            CanonicalDeserializeWithFlags::deserialize_with_flags(reader)?;
        Ok((Fp4::new(c0, c1), flags))
    }
}

impl<P: Fp4Parameters> CanonicalDeserialize for Fp4<P> {
    fn deserialize<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        let c0: Fp2<P::Fp2Params> = CanonicalDeserialize::deserialize(reader)?;
        let c1: Fp2<P::Fp2Params> = CanonicalDeserialize::deserialize(reader)?;
        Ok(Fp4::new(c0, c1))
    }
}
