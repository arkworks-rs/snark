use algebra::{
    curves::{
        short_weierstrass_jacobian::{GroupAffine as SWAffine, GroupProjective as SWProjective},
        SWModelParameters,
    },
    AffineCurve, BigInteger, BitIteratorBE, Field, One, PrimeField, ProjectiveCurve, Zero,
};
use core::{borrow::Borrow, marker::PhantomData};
use r1cs_core::{ConstraintSystemRef, Namespace, SynthesisError};

use crate::fields::fp::FpVar;
use crate::{prelude::*, ToConstraintFieldGadget, Vec};

/// This module provides a generic implementation of G1 and G2 for
/// the [[BLS12]](https://eprint.iacr.org/2002/088.pdf) family of bilinear groups.
pub mod bls12;

/// This module provides a generic implementation of G1 and G2 for
/// the [[MNT4]](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.20.8113&rep=rep1&type=pdf)
///  family of bilinear groups.
pub mod mnt4;
/// This module provides a generic implementation of G1 and G2 for
/// the [[MNT6]](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.20.8113&rep=rep1&type=pdf)
///  family of bilinear groups.
pub mod mnt6;

/// An implementation of arithmetic for Short Weierstrass curves that relies on
/// the complete formulae derived in the paper of
/// [[Renes, Costello, Batina 2015]](https://eprint.iacr.org/2015/1060).
#[derive(Derivative)]
#[derivative(Debug, Clone)]
#[must_use]
pub struct ProjectiveVar<
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
> where
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    /// The x-coordinate.
    pub x: F,
    /// The y-coordinate.
    pub y: F,
    /// The z-coordinate.
    pub z: F,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

/// An affine representation of a curve point.
#[derive(Derivative)]
#[derivative(Debug, Clone)]
#[must_use]
pub struct AffineVar<
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
> where
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    /// The x-coordinate.
    pub x: F,
    /// The y-coordinate.
    pub y: F,
    /// Is `self` the point at infinity.
    pub infinity: Boolean<<P::BaseField as Field>::BasePrimeField>,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

impl<P, F> AffineVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    fn new(x: F, y: F, infinity: Boolean<<P::BaseField as Field>::BasePrimeField>) -> Self {
        Self {
            x,
            y,
            infinity,
            _params: PhantomData,
        }
    }

    /// Returns the value assigned to `self` in the underlying
    /// constraint system.
    pub fn value(&self) -> Result<SWAffine<P>, SynthesisError> {
        Ok(SWAffine::new(
            self.x.value()?,
            self.y.value()?,
            self.infinity.value()?,
        ))
    }
}

impl<P, F> ToConstraintFieldGadget<<P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
    F: ToConstraintFieldGadget<<P::BaseField as Field>::BasePrimeField>,
{
    fn to_constraint_field(
        &self,
    ) -> Result<Vec<FpVar<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let mut res = Vec::<FpVar<<P::BaseField as Field>::BasePrimeField>>::new();

        res.extend_from_slice(&self.x.to_constraint_field()?);
        res.extend_from_slice(&self.y.to_constraint_field()?);

        Ok(res)
    }
}

impl<P, F> R1CSVar<<P::BaseField as Field>::BasePrimeField> for ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    type Value = SWProjective<P>;

    fn cs(&self) -> ConstraintSystemRef<<P::BaseField as Field>::BasePrimeField> {
        self.x.cs().or(self.y.cs()).or(self.z.cs())
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let (x, y, z) = (self.x.value()?, self.y.value()?, self.z.value()?);
        let result = if let Some(z_inv) = z.inverse() {
            SWAffine::new(x * &z_inv, y * &z_inv, false)
        } else {
            SWAffine::zero()
        };
        Ok(result.into())
    }
}

impl<P: SWModelParameters, F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>>
    ProjectiveVar<P, F>
where
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    /// Constructs `Self` from an `(x, y, z)` coordinate triple.
    pub fn new(x: F, y: F, z: F) -> Self {
        Self {
            x,
            y,
            z,
            _params: PhantomData,
        }
    }

    /// Convert this point into affine form.
    #[tracing::instrument(target = "r1cs")]
    pub fn to_affine(&self) -> Result<AffineVar<P, F>, SynthesisError> {
        let cs = self.cs();
        let mode = if self.is_constant() {
            let point = self.value()?.into_affine();
            let x = F::new_constant(ConstraintSystemRef::None, point.x)?;
            let y = F::new_constant(ConstraintSystemRef::None, point.y)?;
            let infinity = Boolean::constant(point.infinity);
            return Ok(AffineVar::new(x, y, infinity));
        } else {
            AllocationMode::Witness
        };

        let infinity = self.is_zero()?;
        let zero_x = F::zero();
        let zero_y = F::one();

        let non_zero_x = F::new_variable(
            r1cs_core::ns!(cs, "non-zero x"),
            || {
                let z_inv = self.z.value()?.inverse().unwrap_or(P::BaseField::zero());
                Ok(self.x.value()? * &z_inv)
            },
            mode,
        )?;
        let non_zero_y = F::new_variable(
            r1cs_core::ns!(cs, "non-zero y"),
            || {
                let z_inv = self.z.value()?.inverse().unwrap_or(P::BaseField::zero());
                Ok(self.y.value()? * &z_inv)
            },
            mode,
        )?;
        let x = infinity.select(&zero_x, &non_zero_x)?;
        let y = infinity.select(&zero_y, &non_zero_y)?;
        Ok(AffineVar::new(x, y, infinity))
    }

    /// Allocates a new variable without performing an on-curve check, which is
    /// useful if the variable is known to be on the curve (eg., if the point
    /// is a constant or is a public input).
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    pub fn new_variable_omit_on_curve_check(
        cs: impl Into<Namespace<<P::BaseField as Field>::BasePrimeField>>,
        f: impl FnOnce() -> Result<SWProjective<P>, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let (x, y, z) = match f() {
            Ok(ge) => {
                let ge = ge.into_affine();
                if ge.is_zero() {
                    (
                        Ok(P::BaseField::zero()),
                        Ok(P::BaseField::one()),
                        Ok(P::BaseField::zero()),
                    )
                } else {
                    (Ok(ge.x), Ok(ge.y), Ok(P::BaseField::one()))
                }
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let x = F::new_variable(r1cs_core::ns!(cs, "x"), || x, mode)?;
        let y = F::new_variable(r1cs_core::ns!(cs, "y"), || y, mode)?;
        let z = F::new_variable(r1cs_core::ns!(cs, "z"), || z, mode)?;

        Ok(Self::new(x, y, z))
    }
}

impl<P, F> CurveVar<SWProjective<P>, <P::BaseField as Field>::BasePrimeField>
    for ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    fn constant(g: SWProjective<P>) -> Self {
        let cs = ConstraintSystemRef::None;
        Self::new_variable_omit_on_curve_check(cs, || Ok(g), AllocationMode::Constant).unwrap()
    }

    fn zero() -> Self {
        Self::new(F::zero(), F::one(), F::zero())
    }

    fn is_zero(&self) -> Result<Boolean<<P::BaseField as Field>::BasePrimeField>, SynthesisError> {
        self.z.is_zero()
    }

    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable_omit_prime_order_check(
        cs: impl Into<Namespace<<P::BaseField as Field>::BasePrimeField>>,
        f: impl FnOnce() -> Result<SWProjective<P>, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        // Curve equation in projective form:
        // E: Y² * Z = X³ + aX * Z² + bZ³
        //
        // This can be re-written as
        // E: Y² * Z - bZ³ = X³ + aX * Z²
        // E: Z * (Y² - bZ²) = X * (X² + aZ²)
        // so, compute X², Y², Z²,
        //     compute temp = X * (X² + aZ²)
        //     check Z.mul_equals((Y² - bZ²), temp)
        //
        //     A total of 5 multiplications

        let g = Self::new_variable_omit_on_curve_check(cs, f, mode)?;

        if mode != AllocationMode::Constant {
            // Perform on-curve check.
            let b = P::COEFF_B;
            let a = P::COEFF_A;

            let x2 = g.x.square()?;
            let y2 = g.y.square()?;
            let z2 = g.z.square()?;
            let t = &g.x * (x2 + &z2 * a);

            g.z.mul_equals(&(y2 - z2 * b), &t)?;
        }
        Ok(g)
    }

    /// Enforce that `self` is in the prime-order subgroup.
    ///
    /// Does so by multiplying by the prime order, and checking that the result
    /// is unchanged.
    // TODO: at the moment this doesn't work, because the addition and doubling
    // formulae are incomplete for even-order points.
    #[tracing::instrument(target = "r1cs")]
    fn enforce_prime_order(&self) -> Result<(), SynthesisError> {
        let r_minus_1 = (-P::ScalarField::one()).into_repr();

        let mut result = Self::zero();
        for b in BitIteratorBE::without_leading_zeros(r_minus_1) {
            result.double_in_place()?;

            if b {
                result += self;
            }
        }
        self.negate()?.enforce_equal(&result)?;
        Ok(())
    }

    #[inline]
    #[tracing::instrument(target = "r1cs")]
    fn double_in_place(&mut self) -> Result<(), SynthesisError> {
        // Complete doubling formula from Renes-Costello-Batina 2015
        // Algorithm 3
        // (https://eprint.iacr.org/2015/1060).
        //
        // Adapted from code in
        // https://github.com/RustCrypto/elliptic-curves/blob/master/p256/src/arithmetic.rs
        let three_b = P::COEFF_B.double() + &P::COEFF_B;

        let xx = self.x.square()?; // 1
        let yy = self.y.square()?; // 2
        let zz = self.z.square()?; // 3
        let xy2 = (&self.x * &self.y).double()?; // 4, 5
        let xz2 = (&self.x * &self.z).double()?; // 6, 7

        let axz2 = mul_by_coeff_a::<P, F>(&xz2); // 8

        let bzz3_part = &axz2 + &zz * three_b; // 9, 10
        let yy_m_bzz3 = &yy - &bzz3_part; // 11
        let yy_p_bzz3 = &yy + &bzz3_part; // 12
        let y_frag = yy_p_bzz3 * &yy_m_bzz3; // 13
        let x_frag = yy_m_bzz3 * &xy2; // 14

        let bxz3 = xz2 * three_b; // 15
        let azz = mul_by_coeff_a::<P, F>(&zz); // 16
        let b3_xz_pairs = mul_by_coeff_a::<P, F>(&(&xx - &azz)) + &bxz3; // 15, 16, 17, 18, 19
        let xx3_p_azz = (xx.double()? + &xx + &azz) * &b3_xz_pairs; // 23, 24, 25

        let y = y_frag + &xx3_p_azz; // 26, 27
        let yz2 = (&self.y * &self.z).double()?; // 28, 29
        let x = x_frag - &(b3_xz_pairs * &yz2); // 30, 31
        let z = (yz2 * &yy).double()?.double()?; // 32, 33, 34
        self.x = x;
        self.y = y;
        self.z = z;
        Ok(())
    }

    #[tracing::instrument(target = "r1cs")]
    fn negate(&self) -> Result<Self, SynthesisError> {
        Ok(Self::new(self.x.clone(), self.y.negate()?, self.z.clone()))
    }
}

fn mul_by_coeff_a<
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
>(
    f: &F,
) -> F
where
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    if !P::COEFF_A.is_zero() {
        f * P::COEFF_A
    } else {
        F::zero()
    }
}

impl_bounded_ops!(
    ProjectiveVar<P, F>,
    SWProjective<P>,
    Add,
    add,
    AddAssign,
    add_assign,
    |this: &'a ProjectiveVar<P, F>, other: &'a ProjectiveVar<P, F>| {
        // Complete addition formula from Renes-Costello-Batina 2015
        // Algorithm 1
        // (https://eprint.iacr.org/2015/1060).
        //
        // Adapted from code in
        // https://github.com/RustCrypto/elliptic-curves/blob/master/p256/src/arithmetic.rs

        let three_b = P::COEFF_B.double() + &P::COEFF_B;

        let xx = &this.x * &other.x; // 1
        let yy = &this.y * &other.y; // 2
        let zz = &this.z * &other.z; // 3
        let xy_pairs = ((&this.x + &this.y) * &(&other.x + &other.y)) - (&xx + &yy); // 4, 5, 6, 7, 8
        let xz_pairs = ((&this.x + &this.z) * &(&other.x + &other.z)) - (&xx + &zz); // 9, 10, 11, 12, 13
        let yz_pairs = ((&this.y + &this.z) * &(&other.y + &other.z)) - (&yy + &zz); // 14, 15, 16, 17, 18

        let axz = mul_by_coeff_a::<P, F>(&xz_pairs); // 19

        let bzz3_part = &axz + &zz * three_b; // 20, 21

        let yy_m_bzz3 = &yy - &bzz3_part; // 22
        let yy_p_bzz3 = &yy + &bzz3_part; // 23

        let azz = mul_by_coeff_a::<P, F>(&zz);
        let xx3_p_azz = xx.double().unwrap() + &xx + &azz; // 25, 26, 27, 29

        let bxz3 = &xz_pairs * three_b; // 28
        let b3_xz_pairs = mul_by_coeff_a::<P, F>(&(&xx - &azz)) + &bxz3; // 30, 31, 32

        let x = (&yy_m_bzz3 * &xy_pairs) - &yz_pairs * &b3_xz_pairs; // 35, 39, 40
        let y = (&yy_p_bzz3 * &yy_m_bzz3) + &xx3_p_azz * b3_xz_pairs; // 24, 36, 37, 38
        let z = (&yy_p_bzz3 * &yz_pairs) + xy_pairs * xx3_p_azz; // 41, 42, 43

        ProjectiveVar::new(x, y, z)
    },
    |this: &'a ProjectiveVar<P, F>, other: SWProjective<P>| {
        this + ProjectiveVar::constant(other)
    },
    (F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>, P: SWModelParameters),
    for <'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
);

impl_bounded_ops!(
    ProjectiveVar<P, F>,
    SWProjective<P>,
    Sub,
    sub,
    SubAssign,
    sub_assign,
    |this: &'a ProjectiveVar<P, F>, other: &'a ProjectiveVar<P, F>| this + other.negate().unwrap(),
    |this: &'a ProjectiveVar<P, F>, other: SWProjective<P>| this - ProjectiveVar::constant(other),
    (F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>, P: SWModelParameters),
    for <'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>
);

impl<'a, P, F> GroupOpsBounds<'a, SWProjective<P>, ProjectiveVar<P, F>> for ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
{
}

impl<'a, P, F> GroupOpsBounds<'a, SWProjective<P>, ProjectiveVar<P, F>> for &'a ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
{
}

impl<P, F> CondSelectGadget<<P::BaseField as Field>::BasePrimeField> for ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    #[inline]
    #[tracing::instrument(target = "r1cs")]
    fn conditionally_select(
        cond: &Boolean<<P::BaseField as Field>::BasePrimeField>,
        true_value: &Self,
        false_value: &Self,
    ) -> Result<Self, SynthesisError> {
        let x = cond.select(&true_value.x, &false_value.x)?;
        let y = cond.select(&true_value.y, &false_value.y)?;
        let z = cond.select(&true_value.z, &false_value.z)?;

        Ok(Self::new(x, y, z))
    }
}

impl<P, F> EqGadget<<P::BaseField as Field>::BasePrimeField> for ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    #[tracing::instrument(target = "r1cs")]
    fn is_eq(
        &self,
        other: &Self,
    ) -> Result<Boolean<<P::BaseField as Field>::BasePrimeField>, SynthesisError> {
        let x_equal = (&self.x * &other.z).is_eq(&(&other.x * &self.z))?;
        let y_equal = (&self.y * &other.z).is_eq(&(&other.y * &self.z))?;
        let coordinates_equal = x_equal.and(&y_equal)?;
        let both_are_zero = self.is_zero()?.and(&other.is_zero()?)?;
        both_are_zero.or(&coordinates_equal)
    }

    #[inline]
    #[tracing::instrument(target = "r1cs")]
    fn conditional_enforce_equal(
        &self,
        other: &Self,
        condition: &Boolean<<P::BaseField as Field>::BasePrimeField>,
    ) -> Result<(), SynthesisError> {
        let x_equal = (&self.x * &other.z).is_eq(&(&other.x * &self.z))?;
        let y_equal = (&self.y * &other.z).is_eq(&(&other.y * &self.z))?;
        let coordinates_equal = x_equal.and(&y_equal)?;
        let both_are_zero = self.is_zero()?.and(&other.is_zero()?)?;
        both_are_zero
            .or(&coordinates_equal)?
            .conditional_enforce_equal(&Boolean::Constant(true), condition)?;
        Ok(())
    }

    #[inline]
    #[tracing::instrument(target = "r1cs")]
    fn conditional_enforce_not_equal(
        &self,
        other: &Self,
        condition: &Boolean<<P::BaseField as Field>::BasePrimeField>,
    ) -> Result<(), SynthesisError> {
        let is_equal = self.is_eq(other)?;
        is_equal
            .and(condition)?
            .enforce_equal(&Boolean::Constant(false))
    }
}

impl<P, F> AllocVar<SWAffine<P>, <P::BaseField as Field>::BasePrimeField> for ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    fn new_variable<T: Borrow<SWAffine<P>>>(
        cs: impl Into<Namespace<<P::BaseField as Field>::BasePrimeField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        Self::new_variable(cs, || f().map(|b| b.borrow().into_projective()), mode)
    }
}

impl<P, F> AllocVar<SWProjective<P>, <P::BaseField as Field>::BasePrimeField>
    for ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    fn new_variable<T: Borrow<SWProjective<P>>>(
        cs: impl Into<Namespace<<P::BaseField as Field>::BasePrimeField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        let f = || Ok(*f()?.borrow());
        match mode {
            AllocationMode::Constant => Self::new_variable_omit_prime_order_check(cs, f, mode),
            AllocationMode::Input => Self::new_variable_omit_prime_order_check(cs, f, mode),
            AllocationMode::Witness => {
                // if cofactor.is_even():
                //   divide until you've removed all even factors
                // else:
                //   just directly use double and add.
                let mut power_of_2: u32 = 0;
                let mut cofactor = P::COFACTOR.to_vec();
                while cofactor[0] % 2 == 0 {
                    div2(&mut cofactor);
                    power_of_2 += 1;
                }

                let cofactor_weight = BitIteratorBE::new(cofactor.as_slice())
                    .filter(|b| *b)
                    .count();
                let modulus_minus_1 = (-P::ScalarField::one()).into_repr(); // r - 1
                let modulus_minus_1_weight =
                    BitIteratorBE::new(modulus_minus_1).filter(|b| *b).count();

                // We pick the most efficient method of performing the prime order check:
                // If the cofactor has lower hamming weight than the scalar field's modulus,
                // we first multiply by the inverse of the cofactor, and then, after allocating,
                // multiply by the cofactor. This ensures the resulting point has no cofactors
                //
                // Else, we multiply by the scalar field's modulus and ensure that the result
                // equals the identity.

                let (mut ge, iter) = if cofactor_weight < modulus_minus_1_weight {
                    let ge = Self::new_variable_omit_prime_order_check(
                        r1cs_core::ns!(cs, "Witness without subgroup check with cofactor mul"),
                        || f().map(|g| g.borrow().into_affine().mul_by_cofactor_inv().into()),
                        mode,
                    )?;
                    (
                        ge,
                        BitIteratorBE::without_leading_zeros(cofactor.as_slice()),
                    )
                } else {
                    let ge = Self::new_variable_omit_prime_order_check(
                        r1cs_core::ns!(cs, "Witness without subgroup check with `r` check"),
                        || {
                            f().map(|g| {
                                let g = g.into_affine();
                                let mut power_of_two = P::ScalarField::one().into_repr();
                                power_of_two.muln(power_of_2);
                                let power_of_two_inv = P::ScalarField::from_repr(power_of_two)
                                    .and_then(|n| n.inverse())
                                    .unwrap();
                                g.mul(power_of_two_inv)
                            })
                        },
                        mode,
                    )?;

                    (
                        ge,
                        BitIteratorBE::without_leading_zeros(modulus_minus_1.as_ref()),
                    )
                };
                // Remove the even part of the cofactor
                for _ in 0..power_of_2 {
                    ge.double_in_place()?;
                }

                let mut result = Self::zero();
                for b in iter {
                    result.double_in_place()?;

                    if b {
                        result += &ge
                    }
                }
                if cofactor_weight < modulus_minus_1_weight {
                    Ok(result)
                } else {
                    ge.enforce_equal(&ge)?;
                    Ok(ge)
                }
            }
        }
    }
}

#[inline]
fn div2(limbs: &mut [u64]) {
    let mut t = 0;
    for i in limbs.iter_mut().rev() {
        let t2 = *i << 63;
        *i >>= 1;
        *i |= t;
        t = t2;
    }
}

impl<P, F> ToBitsGadget<<P::BaseField as Field>::BasePrimeField> for ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    #[tracing::instrument(target = "r1cs")]
    fn to_bits_le(
        &self,
    ) -> Result<Vec<Boolean<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let g = self.to_affine()?;
        let mut bits = g.x.to_bits_le()?;
        let y_bits = g.y.to_bits_le()?;
        bits.extend_from_slice(&y_bits);
        bits.push(g.infinity);
        Ok(bits)
    }

    #[tracing::instrument(target = "r1cs")]
    fn to_non_unique_bits_le(
        &self,
    ) -> Result<Vec<Boolean<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let g = self.to_affine()?;
        let mut bits = g.x.to_non_unique_bits_le()?;
        let y_bits = g.y.to_non_unique_bits_le()?;
        bits.extend_from_slice(&y_bits);
        bits.push(g.infinity);
        Ok(bits)
    }
}

impl<P, F> ToBytesGadget<<P::BaseField as Field>::BasePrimeField> for ProjectiveVar<P, F>
where
    P: SWModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    #[tracing::instrument(target = "r1cs")]
    fn to_bytes(
        &self,
    ) -> Result<Vec<UInt8<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let g = self.to_affine()?;
        let mut bytes = g.x.to_bytes()?;
        let y_bytes = g.y.to_bytes()?;
        let inf_bytes = g.infinity.to_bytes()?;
        bytes.extend_from_slice(&y_bytes);
        bytes.extend_from_slice(&inf_bytes);
        Ok(bytes)
    }

    #[tracing::instrument(target = "r1cs")]
    fn to_non_unique_bytes(
        &self,
    ) -> Result<Vec<UInt8<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let g = self.to_affine()?;
        let mut bytes = g.x.to_non_unique_bytes()?;
        let y_bytes = g.y.to_non_unique_bytes()?;
        let inf_bytes = g.infinity.to_non_unique_bytes()?;
        bytes.extend_from_slice(&y_bytes);
        bytes.extend_from_slice(&inf_bytes);
        Ok(bytes)
    }
}

#[cfg(test)]
#[allow(dead_code)]
pub(crate) fn test<P, GG>() -> Result<(), SynthesisError>
where
    P: SWModelParameters,
    GG: CurveVar<SWProjective<P>, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a GG: GroupOpsBounds<'a, SWProjective<P>, GG>,
{
    use crate::prelude::*;
    use algebra::{test_rng, BitIteratorLE, Group, UniformRand};
    use r1cs_core::ConstraintSystem;

    crate::groups::test::group_test::<SWProjective<P>, _, GG>()?;

    let mut rng = test_rng();

    let cs = ConstraintSystem::<<P::BaseField as Field>::BasePrimeField>::new_ref();

    let a = SWProjective::<P>::rand(&mut rng);
    let b = SWProjective::<P>::rand(&mut rng);
    let a_affine = a.into_affine();
    let b_affine = b.into_affine();

    println!("Allocating things");
    let ns = r1cs_core::ns!(cs, "allocating variables");
    let mut gadget_a = GG::new_witness(cs.clone(), || Ok(a))?;
    let gadget_b = GG::new_witness(cs.clone(), || Ok(b))?;
    drop(ns);
    println!("Done Allocating things");
    assert_eq!(gadget_a.value()?.into_affine().x, a_affine.x);
    assert_eq!(gadget_a.value()?.into_affine().y, a_affine.y);
    assert_eq!(gadget_b.value()?.into_affine().x, b_affine.x);
    assert_eq!(gadget_b.value()?.into_affine().y, b_affine.y);
    assert_eq!(cs.which_is_unsatisfied().unwrap(), None);

    println!("Checking addition");
    // Check addition
    let ab = a + &b;
    let ab_affine = ab.into_affine();
    let gadget_ab = &gadget_a + &gadget_b;
    let gadget_ba = &gadget_b + &gadget_a;
    gadget_ba.enforce_equal(&gadget_ab)?;

    let ab_val = gadget_ab.value()?.into_affine();
    assert_eq!(ab_val, ab_affine, "Result of addition is unequal");
    assert!(cs.is_satisfied().unwrap());
    println!("Done checking addition");

    println!("Checking doubling");
    // Check doubling
    let aa = Group::double(&a);
    let aa_affine = aa.into_affine();
    gadget_a.double_in_place()?;
    let aa_val = gadget_a.value()?.into_affine();
    assert_eq!(
        aa_val, aa_affine,
        "Gadget and native values are unequal after double."
    );
    assert!(cs.is_satisfied().unwrap());
    println!("Done checking doubling");

    println!("Checking mul_bits");
    // Check mul_bits
    let scalar = P::ScalarField::rand(&mut rng);
    let native_result = aa.into_affine().mul(scalar);
    let native_result = native_result.into_affine();

    let scalar: Vec<bool> = BitIteratorLE::new(scalar.into_repr()).collect();
    let input: Vec<Boolean<_>> =
        Vec::new_witness(r1cs_core::ns!(cs, "bits"), || Ok(scalar)).unwrap();
    let result = gadget_a.scalar_mul_le(input.iter())?;
    let result_val = result.value()?.into_affine();
    assert_eq!(
        result_val, native_result,
        "gadget & native values are diff. after scalar mul"
    );
    assert!(cs.is_satisfied().unwrap());
    println!("Done checking mul_bits");

    if !cs.is_satisfied().unwrap() {
        println!("Not satisfied");
        println!("{:?}", cs.which_is_unsatisfied().unwrap());
    }

    assert!(cs.is_satisfied().unwrap());
    Ok(())
}
