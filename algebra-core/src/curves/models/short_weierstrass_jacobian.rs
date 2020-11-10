use crate::{
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

#[cfg(not(feature = "cuda"))]
use crate::accel_dummy::*;
#[cfg(feature = "cuda")]
use accel::*;

#[cfg(feature = "cuda")]
use {
    crate::curves::BatchGroupArithmeticSlice, closure::closure, log::debug, peekmore::PeekMore,
    std::sync::Mutex,
};

use crate::{
    bytes::{FromBytes, ToBytes},
    cfg_chunks_mut, cfg_iter,
    curves::{
        cuda::scalar_mul::{internal::GPUScalarMulInternal, ScalarMulProfiler},
        AffineCurve, BatchGroupArithmetic, ModelParameters, ProjectiveCurve,
    },
    fields::{BitIteratorBE, Field, FpParameters, PrimeField, SquareRootField},
    impl_gpu_cpu_run_kernel, impl_gpu_sw_projective, impl_run_kernel,
};

use crate::{
    CanonicalDeserialize, CanonicalDeserializeWithFlags, CanonicalSerialize,
    CanonicalSerializeWithFlags, ConstantSerializedSize,
};

specialise_affine_to_proj!(GroupProjective);

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub trait SWModelParameters: ModelParameters + Sized {
    const COEFF_A: Self::BaseField;
    const COEFF_B: Self::BaseField;
    const COFACTOR: &'static [u64];
    const COFACTOR_INV: Self::ScalarField;
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField);

    #[inline(always)]
    fn mul_by_a(elem: &Self::BaseField) -> Self::BaseField {
        let mut copy = *elem;
        copy *= &Self::COEFF_A;
        copy
    }

    #[inline(always)]
    fn glv_window_size() -> usize {
        4
    }

    #[inline(always)]
    fn add_b(elem: &Self::BaseField) -> Self::BaseField {
        let mut copy = *elem;
        copy += &Self::COEFF_B;
        copy
    }

    #[inline(always)]
    fn has_glv() -> bool {
        false
    }

    #[inline(always)]
    fn glv_endomorphism_in_place(_elem: &mut Self::BaseField) {
        unimplemented!()
    }

    #[inline(always)]
    fn glv_scalar_decomposition(
        _k: <Self::ScalarField as PrimeField>::BigInt,
    ) -> (
        (bool, <Self::ScalarField as PrimeField>::BigInt),
        (bool, <Self::ScalarField as PrimeField>::BigInt),
    ) {
        unimplemented!()
    }

    fn scalar_mul_kernel(
        ctx: &Context,
        grid: usize,
        block: usize,
        table: *const GroupProjective<Self>,
        exps: *const u8,
        out: *mut GroupProjective<Self>,
        n: isize,
    ) -> error::Result<()>;

    fn scalar_mul_static_profiler() -> ScalarMulProfiler;

    fn namespace() -> &'static str;
}

impl_gpu_sw_projective!(SWModelParameters);

#[derive(Derivative)]
#[derivative(
    Copy(bound = "P: SWModelParameters"),
    Clone(bound = "P: SWModelParameters"),
    Eq(bound = "P: SWModelParameters"),
    Debug(bound = "P: SWModelParameters"),
    Hash(bound = "P: SWModelParameters")
)]
pub struct GroupProjective<P: SWModelParameters> {
    pub x: P::BaseField,
    pub y: P::BaseField,
    pub z: P::BaseField,
    _params: PhantomData<P>,
}

impl<P: SWModelParameters> GroupProjective<P> {
    #[inline(always)]
    pub fn has_glv() -> bool {
        P::has_glv()
    }

    #[inline(always)]
    pub fn glv_endomorphism_in_place(elem: &mut <Self as ProjectiveCurve>::BaseField) {
        P::glv_endomorphism_in_place(elem);
    }

    #[inline]
    pub fn glv_scalar_decomposition(
        k: <<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt,
    ) -> (
        (
            bool,
            <<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt,
        ),
        (
            bool,
            <<Self as ProjectiveCurve>::ScalarField as PrimeField>::BigInt,
        ),
    ) {
        P::glv_scalar_decomposition(k)
    }
}

impl<P: SWModelParameters> Display for GroupProjective<P> {
    fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
        write!(f, "{}", GroupAffine::from(*self))
    }
}

impl<P: SWModelParameters> PartialEq for GroupProjective<P> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() {
            return other.is_zero();
        }

        if other.is_zero() {
            return false;
        }

        // The points (X, Y, Z) and (X', Y', Z')
        // are equal when (X * Z^2) = (X' * Z'^2)
        // and (Y * Z^3) = (Y' * Z'^3).
        let z1 = self.z.square();
        let z2 = other.z.square();

        if self.x * &z2 != other.x * &z1 {
            false
        } else {
            self.y * &(z2 * &other.z) == other.y * &(z1 * &self.z)
        }
    }
}

impl<P: SWModelParameters> Distribution<GroupProjective<P>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> GroupProjective<P> {
        let mut res = GroupProjective::prime_subgroup_generator();
        res.mul_assign(P::ScalarField::rand(rng));
        debug_assert!(GroupAffine::from(res).is_in_correct_subgroup_assuming_on_curve());
        res
    }
}

impl<P: SWModelParameters> ToBytes for GroupProjective<P> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.x.write(&mut writer)?;
        self.y.write(&mut writer)?;
        self.z.write(writer)
    }
}

impl<P: SWModelParameters> FromBytes for GroupProjective<P> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let x = P::BaseField::read(&mut reader)?;
        let y = P::BaseField::read(&mut reader)?;
        let z = P::BaseField::read(reader)?;
        Ok(Self::new(x, y, z))
    }
}

impl<P: SWModelParameters> Default for GroupProjective<P> {
    #[inline]
    fn default() -> Self {
        Self::zero()
    }
}

impl<P: SWModelParameters> GroupProjective<P> {
    pub fn new(x: P::BaseField, y: P::BaseField, z: P::BaseField) -> Self {
        Self {
            x,
            y,
            z,
            _params: PhantomData,
        }
    }
}

impl<P: SWModelParameters> Zero for GroupProjective<P> {
    // The point at infinity is always represented by
    // Z = 0.
    #[inline]
    fn zero() -> Self {
        Self::new(
            P::BaseField::one(),
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

impl<P: SWModelParameters> ProjectiveCurve for GroupProjective<P> {
    const COFACTOR: &'static [u64] = P::COFACTOR;
    type BaseField = P::BaseField;
    type ScalarField = P::ScalarField;
    type Affine = GroupAffine<P>;

    #[inline]
    fn prime_subgroup_generator() -> Self {
        GroupAffine::prime_subgroup_generator().into()
    }

    #[inline]
    fn is_normalized(&self) -> bool {
        self.is_zero() || self.z.is_one()
    }

    #[inline]
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

        #[cfg(not(feature = "parallel"))]
        let v_iter = v.iter_mut();
        #[cfg(feature = "parallel")]
        let v_iter = v.par_iter_mut();

        // Perform affine transformations
        v_iter.filter(|g| !g.is_normalized()).for_each(|g| {
            let z2 = g.z.square(); // 1/z
            g.x *= &z2; // x/z^2
            g.y *= &(z2 * &g.z); // y/z^3
            g.z = P::BaseField::one(); // z = 1
        });
    }

    fn double_in_place(&mut self) -> &mut Self {
        if self.is_zero() {
            return self;
        }

        if P::COEFF_A.is_zero() {
            // A = X1^2
            let mut a = self.x.square();

            // B = Y1^2
            let b = self.y.square();

            // C = B^2
            let mut c = b.square();

            // D = 2*((X1+B)2-A-C)
            let d = ((self.x + &b).square() - &a - &c).double();

            // E = 3*A
            let e = a + &*a.double_in_place();

            // F = E^2
            let f = e.square();

            // Z3 = 2*Y1*Z1
            self.z *= &self.y;
            self.z.double_in_place();

            // X3 = F-2*D
            self.x = f - &d - &d;

            // Y3 = E*(D-X3)-8*C
            self.y = (d - &self.x) * &e - &*c.double_in_place().double_in_place().double_in_place();
            self
        } else {
            // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
            // XX = X1^2
            let xx = self.x.square();

            // YY = Y1^2
            let yy = self.y.square();

            // YYYY = YY^2
            let mut yyyy = yy.square();

            // ZZ = Z1^2
            let zz = self.z.square();

            // S = 2*((X1+YY)^2-XX-YYYY)
            let s = ((self.x + &yy).square() - &xx - &yyyy).double();

            // M = 3*XX+a*ZZ^2
            let m = xx + &xx + &xx + &P::mul_by_a(&zz.square());

            // T = M^2-2*S
            let t = m.square() - &s.double();

            // X3 = T
            self.x = t;
            // Y3 = M*(S-T)-8*YYYY
            let old_y = self.y;
            self.y = m * &(s - &t) - &*yyyy.double_in_place().double_in_place().double_in_place();
            // Z3 = (Y1+Z1)^2-YY-ZZ
            self.z = (old_y + &self.z).square() - &yy - &zz;
            self
        }
    }

    fn add_assign_mixed(&mut self, other: &GroupAffine<P>) {
        if other.is_zero() {
            return;
        }

        if self.is_zero() {
            self.x = other.x;
            self.y = other.y;
            self.z = P::BaseField::one();
            return;
        }

        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-madd-2007-bl
        // Works for all curves.

        // Z1Z1 = Z1^2
        let z1z1 = self.z.square();

        // U2 = X2*Z1Z1
        let u2 = other.x * &z1z1;

        // S2 = Y2*Z1*Z1Z1
        let s2 = (other.y * &self.z) * &z1z1;

        if self.x == u2 && self.y == s2 {
            // The two points are equal, so we double.
            self.double_in_place();
        } else {
            // If we're adding -a and a together, self.z becomes zero as H becomes zero.

            // H = U2-X1
            let h = u2 - &self.x;

            // HH = H^2
            let hh = h.square();

            // I = 4*HH
            let mut i = hh;
            i.double_in_place().double_in_place();

            // J = H*I
            let mut j = h * &i;

            // r = 2*(S2-Y1)
            let r = (s2 - &self.y).double();

            // V = X1*I
            let v = self.x * &i;

            // X3 = r^2 - J - 2*V
            self.x = r.square();
            self.x -= &j;
            self.x -= &v;
            self.x -= &v;

            // Y3 = r*(V-X3)-2*Y1*J
            j *= &self.y; // J = 2*Y1*J
            j.double_in_place();
            self.y = v - &self.x;
            self.y *= &r;
            self.y -= &j;

            // Z3 = (Z1+H)^2-Z1Z1-HH
            self.z += &h;
            self.z.square_in_place();
            self.z -= &z1z1;
            self.z -= &hh;
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

    fn get_x(&mut self) -> &mut Self::BaseField {
        &mut self.x
    }
}

impl<P: SWModelParameters> Neg for GroupProjective<P> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self {
        if !self.is_zero() {
            Self::new(self.x, -self.y, self.z)
        } else {
            self
        }
    }
}

crate::impl_additive_ops_from_ref!(GroupProjective, SWModelParameters);

impl<'a, P: SWModelParameters> Add<&'a Self> for GroupProjective<P> {
    type Output = Self;

    #[inline]
    fn add(self, other: &'a Self) -> Self {
        let mut copy = self;
        copy += other;
        copy
    }
}

impl<'a, P: SWModelParameters> AddAssign<&'a Self> for GroupProjective<P> {
    fn add_assign(&mut self, other: &'a Self) {
        if self.is_zero() {
            *self = *other;
            return;
        }

        if other.is_zero() {
            return;
        }

        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
        // Works for all curves.

        // Z1Z1 = Z1^2
        let z1z1 = self.z.square();

        // Z2Z2 = Z2^2
        let z2z2 = other.z.square();

        // U1 = X1*Z2Z2
        let u1 = self.x * &z2z2;

        // U2 = X2*Z1Z1
        let u2 = other.x * &z1z1;

        // S1 = Y1*Z2*Z2Z2
        let s1 = self.y * &other.z * &z2z2;

        // S2 = Y2*Z1*Z1Z1
        let s2 = other.y * &self.z * &z1z1;

        if u1 == u2 && s1 == s2 {
            // The two points are equal, so we double.
            self.double_in_place();
        } else {
            // If we're adding -a and a together, self.z becomes zero as H becomes zero.

            // H = U2-U1
            let h = u2 - &u1;

            // I = (2*H)^2
            let i = (h.double()).square();

            // J = H*I
            let j = h * &i;

            // r = 2*(S2-S1)
            let r = (s2 - &s1).double();

            // V = U1*I
            let v = u1 * &i;

            // X3 = r^2 - J - 2*V
            self.x = r.square() - &j - &(v.double());

            // Y3 = r*(V - X3) - 2*S1*J
            self.y = r * &(v - &self.x) - &(s1 * &j).double();

            // Z3 = ((Z1+Z2)^2 - Z1Z1 - Z2Z2)*H
            self.z = ((self.z + &other.z).square() - &z1z1 - &z2z2) * &h;
        }
    }
}

impl<'a, P: SWModelParameters> Sub<&'a Self> for GroupProjective<P> {
    type Output = Self;

    #[inline]
    fn sub(self, other: &'a Self) -> Self {
        let mut copy = self;
        copy -= other;
        copy
    }
}

impl<'a, P: SWModelParameters> SubAssign<&'a Self> for GroupProjective<P> {
    fn sub_assign(&mut self, other: &'a Self) {
        *self += &(-(*other));
    }
}

impl<P: SWModelParameters> MulAssign<P::ScalarField> for GroupProjective<P> {
    fn mul_assign(&mut self, other: P::ScalarField) {
        *self = self.mul(other.into_repr())
    }
}

// The affine point X, Y is represented in the Jacobian
// coordinates with Z = 1.
impl<P: SWModelParameters> From<GroupAffine<P>> for GroupProjective<P> {
    #[inline]
    fn from(p: GroupAffine<P>) -> GroupProjective<P> {
        if p.is_zero() {
            Self::zero()
        } else {
            Self::new(p.x, p.y, P::BaseField::one())
        }
    }
}

// The projective point X, Y, Z is represented in the affine
// coordinates as X/Z^2, Y/Z^3.
impl<P: SWModelParameters> From<GroupProjective<P>> for GroupAffine<P> {
    #[inline]
    fn from(p: GroupProjective<P>) -> GroupAffine<P> {
        if p.is_zero() {
            GroupAffine::zero()
        } else if p.z.is_one() {
            // If Z is one, the point is already normalized.
            GroupAffine::new(p.x, p.y, false)
        } else {
            // Z is nonzero, so it must have an inverse in a field.
            let zinv = p.z.inverse().unwrap();
            let zinv_squared = zinv.square();

            // X/Z^2
            let x = p.x * &zinv_squared;

            // Y/Z^3
            let y = p.y * &(zinv_squared * &zinv);

            GroupAffine::new(x, y, false)
        }
    }
}
