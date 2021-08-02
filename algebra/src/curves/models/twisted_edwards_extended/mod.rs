use rand::{Rng, distributions::{Standard, Distribution}};
use std::{
    fmt::{Display, Formatter, Result as FmtResult},
    io::{Read, Result as IoResult, Error as IoError, Write, ErrorKind},
    marker::PhantomData,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{
    bytes::{FromBytes, ToBytes},
    curves::{
        models::TEModelParameters as Parameters, models::MontgomeryModelParameters as MontgomeryParameters,
        AffineCurve, ProjectiveCurve
    }, fields::{BitIterator, Field, PrimeField, SquareRootField}, UniformRand, SemanticallyValid,
    Error, FromBytesChecked, ToCompressedBits, FromCompressedBits, BitSerializationError,
    CanonicalSerialize, SerializationError, CanonicalSerializeWithFlags, CanonicalDeserialize,
    CanonicalDeserializeWithFlags, EdwardsFlags
};
use serde::{Serialize, Deserialize};

#[cfg(test)]
pub mod tests;

#[derive(Derivative)]
#[derivative(
    Copy(bound = "P: Parameters"),
    Clone(bound = "P: Parameters"),
    PartialEq(bound = "P: Parameters"),
    Eq(bound = "P: Parameters"),
    Debug(bound = "P: Parameters"),
    Hash(bound = "P: Parameters")
)]
#[derive(Serialize, Deserialize)]
pub struct GroupAffine<P: Parameters> {
    pub x: P::BaseField,
    pub y: P::BaseField,
    #[derivative(Debug = "ignore")]
    #[serde(skip)]
    _params: PhantomData<P>,
}

impl<P: Parameters> PartialEq<GroupProjective<P>> for GroupAffine<P> {
    fn eq(&self, other: &GroupProjective<P>) -> bool {
        self.into_projective() == *other
    }
}

impl<P: Parameters> PartialEq<GroupAffine<P>> for GroupProjective<P> {
    fn eq(&self, other: &GroupAffine<P>) -> bool {
        *self == other.into_projective()
    }
}

impl<P: Parameters> Display for GroupAffine<P> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "GroupAffine(x={}, y={})", self.x, self.y)
    }
}

impl<P: Parameters> GroupAffine<P> {
    pub fn new(x: P::BaseField, y: P::BaseField) -> Self {
        Self {
            x,
            y,
            _params: PhantomData,
        }
    }

    #[must_use]
    pub fn scale_by_cofactor(&self) -> <Self as AffineCurve>::Projective {
        self.mul_bits(BitIterator::new(P::COFACTOR))
    }

    /// WARNING: This implementation doesn't take costant time with respect
    /// to the exponent, and therefore is susceptible to side-channel attacks.
    /// Be sure to use it in applications where timing (or similar) attacks
    /// are not possible.
    /// TODO: Add a side-channel secure variant.
    #[must_use]
    pub(crate) fn mul_bits<S: AsRef<[u64]>>(
        &self,
        bits: BitIterator<S>,
    ) -> <Self as AffineCurve>::Projective {
        let mut res = GroupProjective::zero();
        for i in bits {
            res.double_in_place();
            if i {
                res.add_assign_mixed(self)
            }
        }
        res
    }

    /// Attempts to construct an affine point given an x-coordinate. The
    /// point is not guaranteed to be in the prime order subgroup.
    ///
    /// If and only if `greatest` is set will the lexicographically
    /// largest y-coordinate be selected.
    #[allow(dead_code)]
    pub(crate) fn get_point_from_x(x: P::BaseField, greatest: bool) -> Option<Self> {
        let x2 = x.square();
        let one = P::BaseField::one();
        let numerator = P::mul_by_a(&x2) - &one;
        let denominator = P::COEFF_D * &x2 - &one;
        let y2 = denominator.inverse().map(|denom| denom * &numerator);
        y2.and_then(|y2| y2.sqrt()).map(|y| {
            let negy = -y;
            let y = if (y < negy) ^ greatest { y } else { negy };
            Self::new(x, y)
        })
    }

    /// Attempts to construct an affine point given an x-coordinate. The
    /// point is not guaranteed to be in the prime order subgroup.
    ///
    /// If and only if `parity` is set will the odd y-coordinate be selected.
    #[allow(dead_code)]
    pub(crate) fn get_point_from_x_and_parity(x: P::BaseField, parity: bool) -> Option<Self> {
        let x2 = x.square();
        let one = P::BaseField::one();
        let numerator = P::mul_by_a(&x2) - &one;
        let denominator = P::COEFF_D * &x2 - &one;
        let y2 = denominator.inverse().map(|denom| denom * &numerator);
        y2.and_then(|y2| y2.sqrt()).map(|y| {
            let negy = -y;
            let y = if y.is_odd() ^ parity { negy } else { y };
            Self::new(x, y)
        })
    }

    /// Checks that the current point is on the elliptic curve.
    pub fn is_on_curve(&self) -> bool {
        let x2 = self.x.square();
        let y2 = self.y.square();

        let lhs = y2 + &P::mul_by_a(&x2);
        let rhs = P::BaseField::one() + &(P::COEFF_D * &(x2 * &y2));

        lhs == rhs
    }

    #[inline]
    pub fn is_in_correct_subgroup_assuming_on_curve(&self) -> bool {
        self.mul_bits(BitIterator::new(P::ScalarField::characteristic()))
            .is_zero()
    }
}

impl<P: Parameters> AffineCurve for GroupAffine<P> {
    type ScalarField = P::ScalarField;
    type BaseField = P::BaseField;
    type Projective = GroupProjective<P>;

    fn zero() -> Self {
        Self::new(Self::BaseField::zero(), Self::BaseField::one())
    }

    fn prime_subgroup_generator() -> Self {
        Self::new(P::AFFINE_GENERATOR_COEFFS.0, P::AFFINE_GENERATOR_COEFFS.1)
    }

    fn from_random_bytes(bytes: &[u8]) -> Option<Self> {
        P::BaseField::from_random_bytes_with_flags::<EdwardsFlags>(bytes).and_then(|(x, flags)| {
            // if x is valid and is zero, then parse this
            // point as infinity.
            if x.is_zero() {
                Some(Self::zero())
            } else {
                Self::get_point_from_x_and_parity(x, flags.is_odd())
            }
        })
    }

    fn is_zero(&self) -> bool {
        self.x.is_zero() & self.y.is_one()
    }

    fn group_membership_test(&self) -> bool {
        self.is_on_curve() && if !self.is_zero() {
            self.is_in_correct_subgroup_assuming_on_curve()
        } else {
            true
        }
    }

    fn add_points(_: &mut [Vec<Self>]) {
        unimplemented!()
    }

    fn mul<S: Into<<Self::ScalarField as PrimeField>::BigInt>>(&self, by: S) -> GroupProjective<P> {
        self.mul_bits(BitIterator::new(by.into()))
    }

    fn into_projective(&self) -> GroupProjective<P> {
        (*self).into()
    }

    fn mul_by_cofactor(&self) -> Self {
        self.scale_by_cofactor().into()
    }

    fn mul_by_cofactor_inv(&self) -> Self {
        self.mul(P::COFACTOR_INV).into()
    }
}

impl<P: Parameters> SemanticallyValid for GroupAffine<P>
{
    fn is_valid(&self) -> bool {

        self.x.is_valid() &&
        self.y.is_valid() &&
        self.group_membership_test()
    }
}

impl<P: Parameters> Neg for GroupAffine<P> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.x, self.y)
    }
}

impl<'a, P: Parameters> Add<&'a Self> for GroupAffine<P> {
    type Output = Self;
    fn add(self, other: &'a Self) -> Self {
        let mut copy = self;
        copy += other;
        copy
    }
}

impl<'a, P: Parameters> AddAssign<&'a Self> for GroupAffine<P> {
    fn add_assign(&mut self, other: &'a Self) {
        let y1y2 = self.y * &other.y;
        let x1x2 = self.x * &other.x;
        let dx1x2y1y2 = P::COEFF_D * &y1y2 * &x1x2;

        let d1 = P::BaseField::one() + &dx1x2y1y2;
        let d2 = P::BaseField::one() - &dx1x2y1y2;

        let x1y2 = self.x * &other.y;
        let y1x2 = self.y * &other.x;

        self.x = (x1y2 + &y1x2) / &d1;
        self.y = (y1y2 - &P::mul_by_a(&x1x2)) / &d2;
    }
}

impl<'a, P: Parameters> Sub<&'a Self> for GroupAffine<P> {
    type Output = Self;
    fn sub(self, other: &'a Self) -> Self {
        let mut copy = self;
        copy -= other;
        copy
    }
}

impl<'a, P: Parameters> SubAssign<&'a Self> for GroupAffine<P> {
    fn sub_assign(&mut self, other: &'a Self) {
        *self += &(-(*other));
    }
}

impl<'a, P: Parameters> Mul<&'a P::ScalarField> for GroupAffine<P> {
    type Output = Self;
    fn mul(self, other: &'a P::ScalarField) -> Self {
        let mut copy = self;
        copy *= other;
        copy
    }
}

impl<'a, P: Parameters> MulAssign<&'a P::ScalarField> for GroupAffine<P> {
    fn mul_assign(&mut self, other: &'a P::ScalarField) {
        *self = <Self as AffineCurve>::mul(self, other.into_repr()).into_affine();
    }
}

impl<P: Parameters> ToBytes for GroupAffine<P> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.x.write(&mut writer)?;
        self.y.write(&mut writer)
    }
}

impl<P: Parameters> FromBytes for GroupAffine<P> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let x = P::BaseField::read(&mut reader)?;
        let y = P::BaseField::read(&mut reader)?;
        Ok(Self::new(x, y))
    }
}

impl<P: Parameters> FromBytesChecked for GroupAffine<P> {
    #[inline]
    fn read_checked<R: Read>(mut reader: R) -> IoResult<Self> {
        let x = P::BaseField::read_checked(&mut reader)?;
        let y = P::BaseField::read_checked(reader)?;
        let p = Self::new(x, y);
        if !p.group_membership_test() {
            return Err(IoError::new(ErrorKind::InvalidData, "invalid point: group membership test failed"));
        }
        Ok(p)
    }
}

use crate::{ToBits, FromBits};
impl<P: Parameters> ToCompressedBits for GroupAffine<P>
{
    #[inline]
    fn compress(&self) -> Vec<bool> {

        let mut res = self.x.write_bits();

        // Is the y-coordinate the odd one of the two associated with the
        // x-coordinate?
        res.push(self.y.is_odd());

        res
    }
}

impl<P: Parameters> FromCompressedBits for GroupAffine<P>
{
    #[inline]
    fn decompress(compressed: Vec<bool>) -> Result<Self, Error> {
        let len = compressed.len() - 1;
        let parity_flag_set = compressed[len];

        //Mask away the flag bits and try to get the x coordinate
        let x = P::BaseField::read_bits(compressed[0..(len - 1)].to_vec())?;

        //Attempt to get the y coordinate from its parity and x
        match Self::get_point_from_x_and_parity(x, parity_flag_set) {

            //Check p belongs to the subgroup we expect
            Some(p) => {
                if p.is_zero() || p.is_in_correct_subgroup_assuming_on_curve() {
                    Ok(p)
                }
                else {
                    let e = BitSerializationError::NotPrimeOrder;
                    Err(Box::new(e))
                }
            }
            _ => Err(Box::new(BitSerializationError::NotOnCurve)),
        }
    }
}

impl<P: Parameters> Default for GroupAffine<P> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl<P: Parameters> Distribution<GroupAffine<P>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> GroupAffine<P> {
        loop {
            let x = P::BaseField::rand(rng);
            let greatest = rng.gen();

            if let Some(p) = GroupAffine::get_point_from_x(x, greatest) {
                return p.scale_by_cofactor().into();
            }
        }
    }
}

mod group_impl {
    use super::*;
    use crate::groups::Group;

    impl<P: Parameters> Group for GroupAffine<P> {
        type ScalarField = P::ScalarField;
        fn zero() -> Self {
            <Self as AffineCurve>::zero()
        }

        fn is_zero(&self) -> bool {
            <Self as AffineCurve>::is_zero(&self)
        }

        #[inline]
        #[must_use]
        fn double(&self) -> Self {
            let mut tmp = *self;
            tmp += self;
            tmp
        }

        #[inline]
        fn double_in_place(&mut self) -> &mut Self {
            let mut tmp = *self;
            tmp += self;
            *self = tmp;
            self
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

#[derive(Derivative)]
#[derivative(
    Copy(bound = "P: Parameters"),
    Clone(bound = "P: Parameters"),
    Eq(bound = "P: Parameters"),
    Debug(bound = "P: Parameters"),
    Hash(bound = "P: Parameters")
)]
#[derive(Serialize, Deserialize)]
pub struct GroupProjective<P: Parameters> {
    pub x: P::BaseField,
    pub y: P::BaseField,
    pub t: P::BaseField,
    pub z: P::BaseField,
    #[derivative(Debug = "ignore")]
    #[serde(skip)]
    _params: PhantomData<P>,
}

impl<P: Parameters> Display for GroupProjective<P> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{}", self.into_affine())
    }
}

impl<P: Parameters> PartialEq for GroupProjective<P> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() {
            return other.is_zero();
        }

        if other.is_zero() {
            return false;
        }

        // x1/z1 == x2/z2  <==> x1 * z2 == x2 * z1
        (self.x * &other.z) == (other.x * &self.z) && (self.y * &other.z) == (other.y * &self.z)
    }
}

impl<P: Parameters> Distribution<GroupProjective<P>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> GroupProjective<P> {
        loop {
            let x = P::BaseField::rand(rng);
            let greatest = rng.gen();

            if let Some(p) = GroupAffine::get_point_from_x(x, greatest) {
                return p.scale_by_cofactor();
            }
        }
    }
}

impl<P: Parameters> ToBytes for GroupProjective<P> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.x.write(&mut writer)?;
        self.y.write(&mut writer)?;
        self.t.write(&mut writer)?;
        self.z.write(writer)
    }
}

impl<P: Parameters> FromBytes for GroupProjective<P> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let x = P::BaseField::read(&mut reader)?;
        let y = P::BaseField::read(&mut reader)?;
        let t = P::BaseField::read(&mut reader)?;
        let z = P::BaseField::read(reader)?;
        Ok(Self::new(x, y, t, z))
    }
}

impl<P: Parameters> FromBytesChecked for GroupProjective<P> {
    #[inline]
    fn read_checked<R: Read>(mut reader: R) -> IoResult<Self> {
        let x = P::BaseField::read_checked(&mut reader)?;
        let y = P::BaseField::read_checked(&mut reader)?;
        let t = P::BaseField::read_checked(&mut reader)?;
        let z = P::BaseField::read_checked(reader)?;
        let p = Self::new(x, y, t, z);
        if !p.group_membership_test() {
            return Err(IoError::new(ErrorKind::InvalidData, "invalid point: group membership test failed"));
        }
        Ok(p)
    }
}

impl<P: Parameters> Default for GroupProjective<P> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl<P: Parameters> GroupProjective<P> {
    pub fn new(x: P::BaseField, y: P::BaseField, t: P::BaseField, z: P::BaseField) -> Self {
        Self {
            x,
            y,
            t,
            z,
            _params: PhantomData,
        }
    }
}

impl<P: Parameters> ProjectiveCurve for GroupProjective<P> {
    type BaseField = P::BaseField;
    type ScalarField = P::ScalarField;
    type Affine = GroupAffine<P>;

    fn zero() -> Self {
        Self::new(
            P::BaseField::zero(),
            P::BaseField::one(),
            P::BaseField::zero(),
            P::BaseField::one(),
        )
    }

    fn prime_subgroup_generator() -> Self {
        GroupAffine::prime_subgroup_generator().into()
    }

    fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y == self.z && !self.y.is_zero() && self.t.is_zero()
    }

    #[inline]
    fn group_membership_test(&self) -> bool {
        self.into_affine().group_membership_test()
    }

    fn is_normalized(&self) -> bool {
        self.z.is_one()
    }

    fn batch_normalization(v: &mut [Self]) {
        // Montgomery’s Trick and Fast Implementation of Masked AES
        // Genelle, Prouff and Quisquater
        // Section 3.2

        // First pass: compute [a, ab, abc, ...]
        let mut prod = Vec::with_capacity(v.len());
        let mut tmp = P::BaseField::one();
        for g in v.iter_mut()
            // Ignore normalized elements
            .filter(|g| !g.is_normalized())
        {
            tmp.mul_assign(&g.z);
            prod.push(tmp);
        }

        // Invert `tmp`.
        tmp = tmp.inverse().unwrap(); // Guaranteed to be nonzero.

        // Second pass: iterate backwards to compute inverses
        for (g, s) in v.iter_mut()
            // Backwards
            .rev()
                // Ignore normalized elements
                .filter(|g| !g.is_normalized())
                // Backwards, skip last element, fill in one for last term.
                .zip(prod.into_iter().rev().skip(1).chain(Some(P::BaseField::one())))
        {
            // tmp := tmp * g.z; g.z := tmp * s = 1/z
            let newtmp = tmp * &g.z;
            g.z = tmp * &s;
            tmp = newtmp;
        }

        // Perform affine transformations
        for g in v.iter_mut().filter(|g| !g.is_normalized()) {
            g.x *= &g.z; // x/z
            g.y *= &g.z;
            g.t *= &g.z;
            g.z = P::BaseField::one(); // z = 1
        }
    }

    fn double_in_place(&mut self) -> &mut Self {
        let tmp = *self;
        *self += &tmp;
        self
    }

    fn add_assign_mixed(&mut self, other: &Self::Affine) {
        // A = X1*X2
        let a = self.x * &other.x;
        // B = Y1*Y2
        let b = self.y * &other.y;
        // C = T1*d*T2
        let c = P::COEFF_D * &self.t * &other.x * &other.y;
        // D = Z1
        let d = self.z;
        // E = (X1+Y1)*(X2+Y2)-A-B
        let e = (self.x + &self.y) * &(other.x + &other.y) - &a - &b;
        // F = D-C
        let f = d - &c;
        // G = D+C
        let g = d + &c;
        // H = B-a*A
        let h = b - &P::mul_by_a(&a);
        // X3 = E*F
        self.x = e * &f;
        // Y3 = G*H
        self.y = g * &h;
        // T3 = E*H
        self.t = e * &h;
        // Z3 = F*G
        self.z = f * &g;
    }

    /// WARNING: This implementation doesn't take costant time with respect
    /// to the exponent, and therefore is susceptible to side-channel attacks.
    /// Be sure to use it in applications where timing (or similar) attacks
    /// are not possible.
    /// TODO: Add a side-channel secure variant.
    fn mul_assign<S: Into<<Self::ScalarField as PrimeField>::BigInt>>(&mut self, other: S) {
        let mut res = Self::zero();

        let mut found_one = false;

        for i in BitIterator::new(other.into()) {
            if found_one {
                res.double_in_place();
            } else {
                found_one = i;
            }

            if i {
                res.add_assign(self);
            }
        }

        *self = res;
    }

    fn into_affine(&self) -> GroupAffine<P> {
        (*self).into()
    }

    fn recommended_wnaf_for_scalar(scalar: <Self::ScalarField as PrimeField>::BigInt) -> usize {
        P::empirical_recommended_wnaf_for_scalar(scalar)
    }

    fn recommended_wnaf_for_num_scalars(num_scalars: usize) -> usize {
        P::empirical_recommended_wnaf_for_num_scalars(num_scalars)
    }
}

impl<P: Parameters> SemanticallyValid for GroupProjective<P>
{
    fn is_valid(&self) -> bool {

        self.x.is_valid() &&
        self.y.is_valid() &&
        self.z.is_valid() &&
        self.t.is_valid() &&
        self.group_membership_test()
    }
}

impl<P: Parameters> Neg for GroupProjective<P> {
    type Output = Self;
    fn neg(mut self) -> Self {
        self.x = -self.x;
        self.t = -self.t;
        self
    }
}

impl<'a, P: Parameters> Add<&'a Self> for GroupProjective<P> {
    type Output = Self;
    fn add(self, other: &'a Self) -> Self {
        let mut copy = self;
        copy += other;
        copy
    }
}

impl<'a, P: Parameters> AddAssign<&'a Self> for GroupProjective<P> {
    fn add_assign(&mut self, other: &'a Self) {
        // See "Twisted Edwards Curves Revisited"
        // Huseyin Hisil, Kenneth Koon-Ho Wong, Gary Carter, and Ed Dawson
        // 3.1 Unified Addition in E^e

        // A = x1 * x2
        let a = self.x * &other.x;

        // B = y1 * y2
        let b = self.y * &other.y;

        // C = d * t1 * t2
        let c = P::COEFF_D * &self.t * &other.t;

        // D = z1 * z2
        let d = self.z * &other.z;

        // H = B - aA
        let h = b - &P::mul_by_a(&a);

        // E = (x1 + y1) * (x2 + y2) - A - B
        let e = (self.x + &self.y) * &(other.x + &other.y) - &a - &b;

        // F = D - C
        let f = d - &c;

        // G = D + C
        let g = d + &c;

        // x3 = E * F
        self.x = e * &f;

        // y3 = G * H
        self.y = g * &h;

        // t3 = E * H
        self.t = e * &h;

        // z3 = F * G
        self.z = f * &g;
    }
}

impl<'a, P: Parameters> Sub<&'a Self> for GroupProjective<P> {
    type Output = Self;
    fn sub(self, other: &'a Self) -> Self {
        let mut copy = self;
        copy -= other;
        copy
    }
}

impl<'a, P: Parameters> SubAssign<&'a Self> for GroupProjective<P> {
    fn sub_assign(&mut self, other: &'a Self) {
        *self += &(-(*other));
    }
}

impl<'a, P: Parameters> Mul<&'a P::ScalarField> for GroupProjective<P> {
    type Output = Self;
    fn mul(self, other: &'a P::ScalarField) -> Self {
        let mut copy = self;
        copy *= other;
        copy
    }
}

impl<'a, P: Parameters> MulAssign<&'a P::ScalarField> for GroupProjective<P> {
    fn mul_assign(&mut self, other: &'a P::ScalarField) {
        <Self as ProjectiveCurve>::mul_assign(self, other.into_repr());
    }
}

// The affine point (X, Y) is represented in the Extended Projective coordinates
// with Z = 1.
impl<P: Parameters> From<GroupAffine<P>> for GroupProjective<P> {
    fn from(p: GroupAffine<P>) -> GroupProjective<P> {
        Self::new(p.x, p.y, p.x * &p.y, P::BaseField::one())
    }
}

// The projective point X, Y, T, Z is represented in the affine
// coordinates as X/Z, Y/Z.
impl<P: Parameters> From<GroupProjective<P>> for GroupAffine<P> {
    fn from(p: GroupProjective<P>) -> GroupAffine<P> {
        if p.is_zero() {
            GroupAffine::zero()
        } else if p.z.is_one() {
            // If Z is one, the point is already normalized.
            GroupAffine::new(p.x, p.y)
        } else {
            // Z is nonzero, so it must have an inverse in a field.
            let z_inv = p.z.inverse().unwrap();
            let x = p.x * &z_inv;
            let y = p.y * &z_inv;
            GroupAffine::new(x, y)
        }
    }
}

impl<P: Parameters> CanonicalSerialize for GroupAffine<P> {
    #[allow(unused_qualifications)]
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        if self.is_zero() {
            let flags = EdwardsFlags::default();
            // Serialize 0.
            P::BaseField::zero().serialize_with_flags(writer, flags)
        } else {
            let flags = EdwardsFlags::from_y_parity(self.y.is_odd());
            self.x.serialize_with_flags(writer, flags)
        }
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        P::BaseField::zero().serialized_size_with_flags::<EdwardsFlags>()
    }

    #[allow(unused_qualifications)]
    #[inline]
    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.x.serialize_uncompressed(&mut writer)?;
        self.y.serialize_uncompressed(&mut writer)?;
        Ok(())
    }

    #[inline]
    fn uncompressed_size(&self) -> usize {
        // x  + y
        self.x.serialized_size() + self.y.serialized_size()
    }
}

impl<P: Parameters> CanonicalSerialize for GroupProjective<P> {
    #[allow(unused_qualifications)]
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        let aff = GroupAffine::<P>::from(self.clone());
        CanonicalSerialize::serialize(&aff, writer)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        let aff = GroupAffine::<P>::from(self.clone());
        aff.serialized_size()
    }

    #[allow(unused_qualifications)]
    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        let aff = GroupAffine::<P>::from(self.clone());
        aff.serialize_uncompressed(writer)
    }

    #[inline]
    fn uncompressed_size(&self) -> usize {
        let aff = GroupAffine::<P>::from(self.clone());
        aff.uncompressed_size()
    }
}

impl<P: Parameters> CanonicalDeserialize for GroupAffine<P> {
    #[allow(unused_qualifications)]
    fn deserialize<R: Read>(reader: R) -> Result<Self, SerializationError> {
        let p = Self::deserialize_unchecked(reader)?;
        if !p.is_zero() && !p.is_in_correct_subgroup_assuming_on_curve() {
            return Err(SerializationError::InvalidData);
        }
        Ok(p)
    }

    #[allow(unused_qualifications)]
    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let (x, flags): (P::BaseField, EdwardsFlags) =
            CanonicalDeserializeWithFlags::deserialize_with_flags(&mut reader)?;
        if x == P::BaseField::zero() {
            Ok(Self::zero())
        } else {
            let p = GroupAffine::<P>::get_point_from_x_and_parity(x, flags.is_odd())
                .ok_or(SerializationError::InvalidData)?;
            Ok(p)
        }
    }

    #[allow(unused_qualifications)]
    fn deserialize_uncompressed<R: Read>(reader: R) -> Result<Self, SerializationError> {
        let p = Self::deserialize_uncompressed_unchecked(reader)?;

        if !p.group_membership_test() {
            return Err(SerializationError::InvalidData);
        }
        Ok(p)
    }

    #[allow(unused_qualifications)]
    fn deserialize_uncompressed_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let x: P::BaseField = CanonicalDeserialize::deserialize(&mut reader)?;
        let y: P::BaseField = CanonicalDeserialize::deserialize(&mut reader)?;

        let p = GroupAffine::<P>::new(x, y);
        Ok(p)
    }
}

impl<P: Parameters> CanonicalDeserialize for GroupProjective<P> {
    #[allow(unused_qualifications)]
    fn deserialize<R: Read>(reader: R) -> Result<Self, SerializationError> {
        let aff = <GroupAffine<P> as CanonicalDeserialize>::deserialize(reader)?;
        Ok(aff.into())
    }

    #[allow(unused_qualifications)]
    fn deserialize_unchecked<R: Read>(reader: R) -> Result<Self, SerializationError> {
        let aff = <GroupAffine<P> as CanonicalDeserialize>::deserialize_unchecked(reader)?;
        Ok(aff.into())
    }

    #[allow(unused_qualifications)]
    fn deserialize_uncompressed<R: Read>(reader: R) -> Result<Self, SerializationError> {
        let aff = <GroupAffine<P> as CanonicalDeserialize>::deserialize_uncompressed(reader)?;
        Ok(aff.into())
    }

    #[allow(unused_qualifications)]
    fn deserialize_uncompressed_unchecked<R: Read>(reader: R) -> Result<Self, SerializationError> {
        let aff = <GroupAffine<P> as CanonicalDeserialize>::deserialize_uncompressed_unchecked(reader)?;
        Ok(aff.into())
    }
}

#[derive(Derivative)]
#[derivative(
Copy(bound = "P: MontgomeryParameters"),
Clone(bound = "P: MontgomeryParameters"),
PartialEq(bound = "P: MontgomeryParameters"),
Eq(bound = "P: MontgomeryParameters"),
Debug(bound = "P: MontgomeryParameters"),
Hash(bound = "P: MontgomeryParameters")
)]
pub struct MontgomeryGroupAffine<P: MontgomeryParameters> {
    pub x: P::BaseField,
    pub y: P::BaseField,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

impl<P: MontgomeryParameters> Display for MontgomeryGroupAffine<P> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "MontgomeryGroupAffine(x={}, y={})", self.x, self.y)
    }
}

impl<P: MontgomeryParameters> MontgomeryGroupAffine<P> {
    pub fn new(x: P::BaseField, y: P::BaseField) -> Self {
        Self {
            x,
            y,
            _params: PhantomData,
        }
    }
}

