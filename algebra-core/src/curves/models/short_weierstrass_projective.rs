use crate::{
    curves::models::SWModelParameters as Parameters,
    io::{Read, Result as IoResult, Write},
    serialize::{Flags, SWFlags},
    UniformRand, Vec,
};
use core::{
    fmt::{Display, Formatter, Result as FmtResult},
    marker::PhantomData,
    ops::{Add, AddAssign, MulAssign, Neg, Sub, SubAssign},
};
use num_traits::{One, Zero};
use rand::{
    distributions::{Distribution, Standard},
    Rng,
};

use crate::{
    bytes::{FromBytes, ToBytes},
    curves::{AffineCurve, BatchGroupArithmetic, ProjectiveCurve},
    fields::{BitIteratorBE, Field, PrimeField, SquareRootField},
};

use crate::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags, ConstantSerializedSize,
};

#[derive(Derivative)]
#[derivative(
    Copy(bound = "P: Parameters"),
    Clone(bound = "P: Parameters"),
    Eq(bound = "P: Parameters"),
    Debug(bound = "P: Parameters"),
    Hash(bound = "P: Parameters")
)]
#[must_use]
pub struct GroupProjective<P: Parameters> {
    pub x: P::BaseField,
    pub y: P::BaseField,
    pub z: P::BaseField,
    _params: PhantomData<P>,
}

specialise_affine_to_proj!(GroupProjective);

impl<P: Parameters> Display for GroupProjective<P> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{}", GroupAffine::from(*self))
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
        if (self.x * &other.z) != (other.x * &self.z) {
            false
        } else {
            (self.y * &other.z) == (other.y * &self.z)
        }
    }
}

impl<P: Parameters> Distribution<GroupProjective<P>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> GroupProjective<P> {
        let mut res = GroupProjective::prime_subgroup_generator();
        res.mul_assign(P::ScalarField::rand(rng));
        debug_assert!(GroupAffine::from(res).is_in_correct_subgroup_assuming_on_curve());
        res
    }
}

impl<P: Parameters> ToBytes for GroupProjective<P> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.x.write(&mut writer)?;
        self.y.write(&mut writer)?;
        self.z.write(writer)
    }
}

impl<P: Parameters> FromBytes for GroupProjective<P> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let x = P::BaseField::read(&mut reader)?;
        let y = P::BaseField::read(&mut reader)?;
        let z = P::BaseField::read(reader)?;
        Ok(Self::new(x, y, z))
    }
}

impl<P: Parameters> Default for GroupProjective<P> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl<P: Parameters> GroupProjective<P> {
    pub fn new(x: P::BaseField, y: P::BaseField, z: P::BaseField) -> Self {
        Self {
            x,
            y,
            z,
            _params: PhantomData,
        }
    }
}

impl<P: Parameters> Zero for GroupProjective<P> {
    // The point at infinity is always represented by Z = 0.
    #[inline]
    fn zero() -> Self {
        Self::new(
            P::BaseField::zero(),
            P::BaseField::one(),
            P::BaseField::zero(),
        )
    }

    // The point at infinity is always represented by
    // Z = 0.
    #[inline]
    fn is_zero(&self) -> bool {
        self.z.is_zero()
    }
}

impl<P: Parameters> ProjectiveCurve for GroupProjective<P> {
    const COFACTOR: &'static [u64] = P::COFACTOR;
    type BaseField = P::BaseField;
    type ScalarField = P::ScalarField;
    type Affine = GroupAffine<P>;

    fn get_x(&mut self) -> &mut Self::BaseField {
        &mut self.x
    }

    #[inline]
    fn prime_subgroup_generator() -> Self {
        GroupAffine::prime_subgroup_generator().into()
    }

    #[inline]
    fn is_normalized(&self) -> bool {
        self.is_zero() || self.z.is_one()
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
            tmp *= &g.z;
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
            g.x *= &g.z; // x/z^2
            g.y *= &g.z;
            g.z = P::BaseField::one(); // z = 1
        }
    }

    fn double_in_place(&mut self) -> &mut Self {
        if self.is_zero() {
            self
        } else {
            // https://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2007-bl

            // XX = X1^2
            let xx = self.x.square();
            // ZZ = Z1^2
            let zz = self.z.square();
            // w = a*ZZ + 3*XX
            let w = P::mul_by_a(&zz) + &(xx + &xx.double());
            // s = 2*Y1*Z1
            let mut s = self.y * &(self.z);
            s.double_in_place();
            // sss = s^3
            let mut sss = s.square();
            sss *= &s;
            // R = Y1*s
            let r = self.y * &s;
            // RR = R2
            let rr = r.square();
            // B = (X1+R)^2-XX-RR
            let b = (self.x + &r).square() - &xx - &rr;
            // h = w2-2*B
            let h = w.square() - &(b + &b);
            // X3 = h*s
            self.x = h * &s;
            // Y3 = w*(B-h)-2*RR
            self.y = w * &(b - &h) - &(rr + &rr);
            // Z3 = sss
            self.z = sss;

            self
        }
    }

    fn add_assign_mixed(&mut self, other: &GroupAffine<P>) {
        if other.is_zero() {
            return;
        } else if self.is_zero() {
            self.x = other.x;
            self.y = other.y;
            self.z = P::BaseField::one();
            return;
        }
        let mut v = other.x * &self.z;
        let mut u = other.y * &self.z;
        if u == self.y && v == self.x {
            // x1 / z1 == x2 / z2 <==> x1 * z2 == x2 * z1;
            // Here, z2 = 1, so we have x1 == x2 * z1;
            self.double_in_place();
        } else {
            // https://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-madd-1998-cmo
            // u = Y2*Z1-Y1
            u -= &self.y;
            // uu = u^2
            let uu = u.square();
            // v = X2*Z1-X1
            v -= &self.x;
            // vv = v2
            let vv = v.square();
            // vvv = v*vv
            let vvv = v * &vv;
            // r = vv*X1
            let r = vv * &self.x;
            // a = uu*Z1-vvv-2*r
            let a = uu * &self.z - &vvv - &r.double();
            // X3 = v*a
            self.x = v * &a;
            // Y3 = u*(R-A)-vvv*Y1
            self.y = u * &(r - &a) - &(vvv * &self.y);
            // Z3 = vvv*Z1
            self.z = vvv * &self.z;
        }
    }

    fn mul<S: Into<<Self::ScalarField as PrimeField>::BigInt>>(mut self, other: S) -> Self {
        if P::has_glv() {
            let w = P::glv_window_size();
            let mut res = Self::zero();
            impl_glv_mul!(Self, P, w, self, res, other);
            res
        } else {
            let mut res = Self::zero();
            for b in BitIteratorBE::without_leading_zeros(other.into()) {
                res.double_in_place();
                if b {
                    res += self;
                }
            }

            self = res;
            self
        }
    }
}

impl<P: Parameters> Neg for GroupProjective<P> {
    type Output = Self;
    fn neg(self) -> Self {
        if !self.is_zero() {
            Self::new(self.x, -self.y, self.z)
        } else {
            self
        }
    }
}

crate::impl_additive_ops_from_ref!(GroupProjective, Parameters);

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
        if self.is_zero() {
            *self = *other;
            return;
        }

        if other.is_zero() {
            return;
        }
        // https://www.hyperelliptic.org/EFD/g1p/data/shortw/projective/addition/add-1998-cmo-2

        if self == other {
            self.double_in_place();
        } else {
            // Y1Z2 = Y1*Z2
            let y1z2 = self.y * &other.z;
            // X1Z2 = X1*Z2
            let x1z2 = self.x * &other.z;
            // Z1Z2 = Z1*Z2
            let z1z2 = self.z * &other.z;
            // u = Y2*Z1-Y1Z2
            let u = (self.z * &other.y) - &y1z2;
            // uu = u^2
            let uu = u.square();
            // v = X2*Z1-X1Z2
            let v = (self.z * &other.x) - &x1z2;
            // vv = v^2
            let vv = v.square();
            // vvv = v*vv
            let vvv = v * &vv;
            // R = vv*X1Z2
            let r = vv * &x1z2;
            // A = uu*Z1Z2-vvv-2*R
            let a = (uu * &z1z2) - &(vvv + &r + &r);
            // X3 = v*A
            self.x = v * &a;
            // Y3 = u*(R-A)-vvv*Y1Z2
            self.y = ((r - &a) * &u) - &(vvv * &y1z2);
            // Z3 = vvv*Z1Z2
            self.z = vvv * &z1z2;
        }
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

impl<P: Parameters> MulAssign<P::ScalarField> for GroupProjective<P> {
    fn mul_assign(&mut self, other: P::ScalarField) {
        *self = self.mul(other.into_repr())
    }
}

// The affine point X, Y is represented in the jacobian
// coordinates with Z = 1.
impl<P: Parameters> From<GroupAffine<P>> for GroupProjective<P> {
    fn from(p: GroupAffine<P>) -> GroupProjective<P> {
        if p.is_zero() {
            Self::zero()
        } else {
            Self::new(p.x, p.y, P::BaseField::one())
        }
    }
}

// The projective point X, Y, Z is represented in the affine
// coordinates as X/Z, Y/Z.
impl<P: Parameters> From<GroupProjective<P>> for GroupAffine<P> {
    fn from(p: GroupProjective<P>) -> GroupAffine<P> {
        if p.is_zero() {
            GroupAffine::zero()
        } else if p.z.is_one() {
            // If Z is one, the point is already normalized.
            GroupAffine::new(p.x, p.y, false)
        } else {
            // Z is nonzero, so it must have an inverse in a field.
            let z_inv = p.z.inverse().unwrap();
            let x = p.x * &z_inv;
            let y = p.y * &z_inv;
            GroupAffine::new(x, y, false)
        }
    }
}
