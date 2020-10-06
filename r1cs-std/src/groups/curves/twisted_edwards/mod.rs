use algebra::{
    curves::{
        twisted_edwards_extended::{GroupAffine as TEAffine, GroupProjective as TEProjective},
        AffineCurve, MontgomeryModelParameters, ProjectiveCurve, TEModelParameters,
    },
    BigInteger, BitIteratorBE, Field, One, PrimeField, Zero,
};

use r1cs_core::{ConstraintSystemRef, Namespace, SynthesisError};

use crate::{prelude::*, ToConstraintFieldGadget, Vec};

use crate::fields::fp::FpVar;
use core::{borrow::Borrow, marker::PhantomData};

/// An implementation of arithmetic for Montgomery curves that relies on
/// incomplete addition formulae for the affine model, as outlined in the
/// [EFD](https://www.hyperelliptic.org/EFD/g1p/auto-montgom.html).
///
/// This is intended for use primarily for implementing efficient
/// multi-scalar-multiplication in the Bowe-Hopwood-Pedersen hash.
#[derive(Derivative)]
#[derivative(Debug, Clone)]
#[must_use]
pub struct MontgomeryAffineVar<
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
> where
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    /// The x-coordinate.
    pub x: F,
    /// The y-coordinate.
    pub y: F,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

mod montgomery_affine_impl {
    use super::*;
    use algebra::{twisted_edwards_extended::GroupAffine, Field};
    use core::ops::Add;

    impl<P, F> R1CSVar<<P::BaseField as Field>::BasePrimeField> for MontgomeryAffineVar<P, F>
    where
        P: TEModelParameters,
        F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
        for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
    {
        type Value = (P::BaseField, P::BaseField);

        fn cs(&self) -> ConstraintSystemRef<<P::BaseField as Field>::BasePrimeField> {
            self.x.cs().or(self.y.cs())
        }

        fn value(&self) -> Result<Self::Value, SynthesisError> {
            let x = self.x.value()?;
            let y = self.y.value()?;
            Ok((x, y))
        }
    }

    impl<
            P: TEModelParameters,
            F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
        > MontgomeryAffineVar<P, F>
    where
        for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
    {
        /// Constructs `Self` from an `(x, y)` coordinate pair.
        pub fn new(x: F, y: F) -> Self {
            Self {
                x,
                y,
                _params: PhantomData,
            }
        }

        /// Converts a Twisted Edwards curve point to coordinates for the corresponding affine
        /// Montgomery curve point.
        #[tracing::instrument(target = "r1cs")]
        pub fn from_edwards_to_coords(
            p: &TEAffine<P>,
        ) -> Result<(P::BaseField, P::BaseField), SynthesisError> {
            let montgomery_point: GroupAffine<P> = if p.y == P::BaseField::one() {
                GroupAffine::zero()
            } else if p.x == P::BaseField::zero() {
                GroupAffine::new(P::BaseField::zero(), P::BaseField::zero())
            } else {
                let u =
                    (P::BaseField::one() + &p.y) * &(P::BaseField::one() - &p.y).inverse().unwrap();
                let v = u * &p.x.inverse().unwrap();
                GroupAffine::new(u, v)
            };

            Ok((montgomery_point.x, montgomery_point.y))
        }

        /// Converts a Twisted Edwards curve point to coordinates for the corresponding affine
        /// Montgomery curve point.
        #[tracing::instrument(target = "r1cs")]
        pub fn new_witness_from_edwards(
            cs: ConstraintSystemRef<<P::BaseField as Field>::BasePrimeField>,
            p: &TEAffine<P>,
        ) -> Result<Self, SynthesisError> {
            let montgomery_coords = Self::from_edwards_to_coords(p)?;
            let u = F::new_witness(r1cs_core::ns!(cs, "u"), || Ok(montgomery_coords.0))?;
            let v = F::new_witness(r1cs_core::ns!(cs, "v"), || Ok(montgomery_coords.1))?;
            Ok(Self::new(u, v))
        }

        /// Converts `self` into a Twisted Edwards curve point variable.
        #[tracing::instrument(target = "r1cs")]
        pub fn into_edwards(&self) -> Result<AffineVar<P, F>, SynthesisError> {
            let cs = self.cs();
            // Compute u = x / y
            let u = F::new_witness(r1cs_core::ns!(cs, "u"), || {
                let y_inv = self
                    .y
                    .value()?
                    .inverse()
                    .ok_or(SynthesisError::DivisionByZero)?;
                Ok(self.x.value()? * &y_inv)
            })?;

            u.mul_equals(&self.y, &self.x)?;

            let v = F::new_witness(r1cs_core::ns!(cs, "v"), || {
                let mut t0 = self.x.value()?;
                let mut t1 = t0;
                t0 -= &P::BaseField::one();
                t1 += &P::BaseField::one();

                Ok(t0 * &t1.inverse().ok_or(SynthesisError::DivisionByZero)?)
            })?;

            let xplusone = &self.x + P::BaseField::one();
            let xminusone = &self.x - P::BaseField::one();
            v.mul_equals(&xplusone, &xminusone)?;

            Ok(AffineVar::new(u, v))
        }
    }

    impl<'a, P, F> Add<&'a MontgomeryAffineVar<P, F>> for MontgomeryAffineVar<P, F>
    where
        P: TEModelParameters,
        F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
        for<'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
    {
        type Output = MontgomeryAffineVar<P, F>;

        #[tracing::instrument(target = "r1cs")]
        fn add(self, other: &'a Self) -> Self::Output {
            let cs = [&self, other].cs();
            let mode = if cs.is_none() {
                AllocationMode::Constant
            } else {
                AllocationMode::Witness
            };

            let coeff_b = P::MontgomeryModelParameters::COEFF_B;
            let coeff_a = P::MontgomeryModelParameters::COEFF_A;

            let lambda = F::new_variable(
                r1cs_core::ns!(cs, "lambda"),
                || {
                    let n = other.y.value()? - &self.y.value()?;
                    let d = other.x.value()? - &self.x.value()?;
                    Ok(n * &d.inverse().ok_or(SynthesisError::DivisionByZero)?)
                },
                mode,
            )
            .unwrap();
            let lambda_n = &other.y - &self.y;
            let lambda_d = &other.x - &self.x;
            lambda_d.mul_equals(&lambda, &lambda_n).unwrap();

            // Compute x'' = B*lambda^2 - A - x - x'
            let xprime = F::new_variable(
                r1cs_core::ns!(cs, "xprime"),
                || {
                    Ok(lambda.value()?.square() * &coeff_b
                        - &coeff_a
                        - &self.x.value()?
                        - &other.x.value()?)
                },
                mode,
            )
            .unwrap();

            let xprime_lc = &self.x + &other.x + &xprime + coeff_a;
            // (lambda) * (lambda) = (A + x + x' + x'')
            let lambda_b = &lambda * coeff_b;
            lambda_b.mul_equals(&lambda, &xprime_lc).unwrap();

            let yprime = F::new_variable(
                r1cs_core::ns!(cs, "yprime"),
                || {
                    Ok(-(self.y.value()?
                        + &(lambda.value()? * &(xprime.value()? - &self.x.value()?))))
                },
                mode,
            )
            .unwrap();

            let xres = &self.x - &xprime;
            let yres = &self.y + &yprime;
            lambda.mul_equals(&xres, &yres).unwrap();
            MontgomeryAffineVar::new(xprime, yprime)
        }
    }
}

/// An implementation of arithmetic for Twisted Edwards curves that relies on
/// the complete formulae for the affine model, as outlined in the
/// [EFD](https://www.hyperelliptic.org/EFD/g1p/auto-twisted.html).
#[derive(Derivative)]
#[derivative(Debug, Clone)]
#[must_use]
pub struct AffineVar<
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
> where
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    /// The x-coordinate.
    pub x: F,
    /// The y-coordinate.
    pub y: F,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

impl<P: TEModelParameters, F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>>
    AffineVar<P, F>
where
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    /// Constructs `Self` from an `(x, y)` coordinate triple.
    pub fn new(x: F, y: F) -> Self {
        Self {
            x,
            y,
            _params: PhantomData,
        }
    }

    /// Allocates a new variable without performing an on-curve check, which is
    /// useful if the variable is known to be on the curve (eg., if the point
    /// is a constant or is a public input).
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    pub fn new_variable_omit_on_curve_check<T: Into<TEAffine<P>>>(
        cs: impl Into<Namespace<<P::BaseField as Field>::BasePrimeField>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let (x, y) = match f() {
            Ok(ge) => {
                let ge: TEAffine<P> = ge.into();
                (Ok(ge.x), Ok(ge.y))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let x = F::new_variable(r1cs_core::ns!(cs, "x"), || x, mode)?;
        let y = F::new_variable(r1cs_core::ns!(cs, "y"), || y, mode)?;

        Ok(Self::new(x, y))
    }
}

impl<P: TEModelParameters, F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>>
    AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>
        + TwoBitLookupGadget<<P::BaseField as Field>::BasePrimeField, TableConstant = P::BaseField>
        + ThreeBitCondNegLookupGadget<
            <P::BaseField as Field>::BasePrimeField,
            TableConstant = P::BaseField,
        >,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    /// Compute a scalar multiplication of `bases` with respect to `scalars`,
    /// where the elements of `scalars` are length-three slices of bits, and which
    /// such that the first two bits are use to select one of the bases,
    /// while the third bit is used to conditionally negate the selection.
    #[tracing::instrument(target = "r1cs", skip(bases, scalars))]
    pub fn precomputed_base_3_bit_signed_digit_scalar_mul<J>(
        bases: &[impl Borrow<[TEProjective<P>]>],
        scalars: &[impl Borrow<[J]>],
    ) -> Result<Self, SynthesisError>
    where
        J: Borrow<[Boolean<<P::BaseField as Field>::BasePrimeField>]>,
    {
        const CHUNK_SIZE: usize = 3;
        let mut ed_result: Option<AffineVar<P, F>> = None;
        let mut result: Option<MontgomeryAffineVar<P, F>> = None;

        let mut process_segment_result = |result: &MontgomeryAffineVar<P, F>| {
            let sgmt_result = result.into_edwards()?;
            ed_result = match ed_result.as_ref() {
                None => Some(sgmt_result),
                Some(r) => Some(sgmt_result + r),
            };
            Ok::<(), SynthesisError>(())
        };

        // Compute ∏(h_i^{m_i}) for all i.
        for (segment_bits_chunks, segment_powers) in scalars.iter().zip(bases) {
            for (bits, base_power) in segment_bits_chunks
                .borrow()
                .iter()
                .zip(segment_powers.borrow())
            {
                let base_power = base_power.borrow();
                let mut acc_power = *base_power;
                let mut coords = vec![];
                for _ in 0..4 {
                    coords.push(acc_power);
                    acc_power += base_power;
                }

                let bits = bits.borrow().to_bits_le()?;
                if bits.len() != CHUNK_SIZE {
                    return Err(SynthesisError::Unsatisfiable);
                }

                let coords = coords
                    .iter()
                    .map(|p| MontgomeryAffineVar::from_edwards_to_coords(&p.into_affine()))
                    .collect::<Result<Vec<_>, _>>()?;

                let x_coeffs = coords.iter().map(|p| p.0).collect::<Vec<_>>();
                let y_coeffs = coords.iter().map(|p| p.1).collect::<Vec<_>>();

                let precomp = bits[0].and(&bits[1])?;

                let x = F::zero()
                    + x_coeffs[0]
                    + F::from(bits[0].clone()) * (x_coeffs[1] - &x_coeffs[0])
                    + F::from(bits[1].clone()) * (x_coeffs[2] - &x_coeffs[0])
                    + F::from(precomp.clone())
                        * (x_coeffs[3] - &x_coeffs[2] - &x_coeffs[1] + &x_coeffs[0]);

                let y = F::three_bit_cond_neg_lookup(&bits, &precomp, &y_coeffs)?;

                let tmp = MontgomeryAffineVar::new(x, y);
                result = match result.as_ref() {
                    None => Some(tmp),
                    Some(r) => Some(tmp + r),
                };
            }

            process_segment_result(&result.unwrap())?;
            result = None;
        }
        if result.is_some() {
            process_segment_result(&result.unwrap())?;
        }
        Ok(ed_result.unwrap())
    }
}

impl<P, F> R1CSVar<<P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    type Value = TEProjective<P>;

    fn cs(&self) -> ConstraintSystemRef<<P::BaseField as Field>::BasePrimeField> {
        self.x.cs().or(self.y.cs())
    }

    #[inline]
    fn value(&self) -> Result<TEProjective<P>, SynthesisError> {
        let (x, y) = (self.x.value()?, self.y.value()?);
        let result = TEAffine::new(x, y);
        Ok(result.into())
    }
}

impl<P, F> CurveVar<TEProjective<P>, <P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>
        + TwoBitLookupGadget<<P::BaseField as Field>::BasePrimeField, TableConstant = P::BaseField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    fn constant(g: TEProjective<P>) -> Self {
        let cs = ConstraintSystemRef::None;
        Self::new_variable_omit_on_curve_check(cs, || Ok(g), AllocationMode::Constant).unwrap()
    }

    fn zero() -> Self {
        Self::new(F::zero(), F::one())
    }

    fn is_zero(&self) -> Result<Boolean<<P::BaseField as Field>::BasePrimeField>, SynthesisError> {
        self.x.is_zero()?.and(&self.x.is_one()?)
    }

    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable_omit_prime_order_check(
        cs: impl Into<Namespace<<P::BaseField as Field>::BasePrimeField>>,
        f: impl FnOnce() -> Result<TEProjective<P>, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let g = Self::new_variable_omit_on_curve_check(cs, f, mode)?;

        if mode != AllocationMode::Constant {
            let d = P::COEFF_D;
            let a = P::COEFF_A;
            // Check that ax^2 + y^2 = 1 + dx^2y^2
            // We do this by checking that ax^2 - 1 = y^2 * (dx^2 - 1)
            let x2 = g.x.square()?;
            let y2 = g.y.square()?;

            let one = P::BaseField::one();
            let d_x2_minus_one = &x2 * d - one;
            let a_x2_minus_one = &x2 * a - one;

            d_x2_minus_one.mul_equals(&y2, &a_x2_minus_one)?;
        }
        Ok(g)
    }

    /// Enforce that `self` is in the prime-order subgroup.
    ///
    /// Does so by multiplying by the prime order, and checking that the result
    /// is unchanged.
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
        if self.is_constant() {
            let value = self.value()?;
            *self = Self::constant(value.double());
        } else {
            let cs = self.cs();
            let a = P::COEFF_A;

            // xy
            let xy = &self.x * &self.y;
            let x2 = self.x.square()?;
            let y2 = self.y.square()?;

            let a_x2 = &x2 * a;

            // Compute x3 = (2xy) / (ax^2 + y^2)
            let x3 = F::new_witness(r1cs_core::ns!(cs, "x3"), || {
                let t0 = xy.value()?.double();
                let t1 = a * &x2.value()? + &y2.value()?;
                Ok(t0 * &t1.inverse().ok_or(SynthesisError::DivisionByZero)?)
            })?;

            let a_x2_plus_y2 = &a_x2 + &y2;
            let two_xy = xy.double()?;
            x3.mul_equals(&a_x2_plus_y2, &two_xy)?;

            // Compute y3 = (y^2 - ax^2) / (2 - ax^2 - y^2)
            let two = P::BaseField::one().double();
            let y3 = F::new_witness(r1cs_core::ns!(cs, "y3"), || {
                let a_x2 = a * &x2.value()?;
                let t0 = y2.value()? - &a_x2;
                let t1 = two - &a_x2 - &y2.value()?;
                Ok(t0 * &t1.inverse().ok_or(SynthesisError::DivisionByZero)?)
            })?;
            let y2_minus_a_x2 = &y2 - &a_x2;
            let two_minus_ax2_minus_y2 = (&a_x2 + &y2).negate()? + two;

            y3.mul_equals(&two_minus_ax2_minus_y2, &y2_minus_a_x2)?;
            self.x = x3;
            self.y = y3;
        }
        Ok(())
    }

    #[tracing::instrument(target = "r1cs")]
    fn negate(&self) -> Result<Self, SynthesisError> {
        Ok(Self::new(self.x.negate()?, self.y.clone()))
    }

    #[tracing::instrument(target = "r1cs", skip(scalar_bits_with_base_powers))]
    fn precomputed_base_scalar_mul_le<'a, I, B>(
        &mut self,
        scalar_bits_with_base_powers: I,
    ) -> Result<(), SynthesisError>
    where
        I: Iterator<Item = (B, &'a TEProjective<P>)>,
        B: Borrow<Boolean<<P::BaseField as Field>::BasePrimeField>>,
    {
        let scalar_bits_with_base_powers = scalar_bits_with_base_powers
            .map(|(bit, base)| (bit.borrow().clone(), (*base).into()))
            .collect::<Vec<(_, TEProjective<P>)>>();
        let zero = TEProjective::zero();
        for bits_base_powers in scalar_bits_with_base_powers.chunks(2) {
            if bits_base_powers.len() == 2 {
                let bits = [bits_base_powers[0].0.clone(), bits_base_powers[1].0.clone()];
                let base_powers = [&bits_base_powers[0].1, &bits_base_powers[1].1];

                let mut table = [
                    zero,
                    *base_powers[0],
                    *base_powers[1],
                    *base_powers[0] + base_powers[1],
                ];

                TEProjective::batch_normalization(&mut table);
                let x_s = [table[0].x, table[1].x, table[2].x, table[3].x];
                let y_s = [table[0].y, table[1].y, table[2].y, table[3].y];

                let x = F::two_bit_lookup(&bits, &x_s)?;
                let y = F::two_bit_lookup(&bits, &y_s)?;
                *self += Self::new(x, y);
            } else if bits_base_powers.len() == 1 {
                let bit = bits_base_powers[0].0.clone();
                let base_power = bits_base_powers[0].1;
                let new_encoded = &*self + base_power;
                *self = bit.select(&new_encoded, &self)?;
            }
        }

        Ok(())
    }
}

impl<P, F> AllocVar<TEProjective<P>, <P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>
        + TwoBitLookupGadget<<P::BaseField as Field>::BasePrimeField, TableConstant = P::BaseField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<Point: Borrow<TEProjective<P>>>(
        cs: impl Into<Namespace<<P::BaseField as Field>::BasePrimeField>>,
        f: impl FnOnce() -> Result<Point, SynthesisError>,
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
                        result += &ge;
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

impl<P, F> AllocVar<TEAffine<P>, <P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>
        + TwoBitLookupGadget<<P::BaseField as Field>::BasePrimeField, TableConstant = P::BaseField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<Point: Borrow<TEAffine<P>>>(
        cs: impl Into<Namespace<<P::BaseField as Field>::BasePrimeField>>,
        f: impl FnOnce() -> Result<Point, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        Self::new_variable(cs, || f().map(|b| b.borrow().into_projective()), mode)
    }
}

impl<P, F> ToConstraintFieldGadget<<P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a F: FieldOpsBounds<'a, P::BaseField, F>,
    F: ToConstraintFieldGadget<<P::BaseField as Field>::BasePrimeField>,
{
    fn to_constraint_field(
        &self,
    ) -> Result<Vec<FpVar<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let mut res = Vec::new();

        res.extend_from_slice(&self.x.to_constraint_field()?);
        res.extend_from_slice(&self.y.to_constraint_field()?);

        Ok(res)
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

impl_bounded_ops!(
    AffineVar<P, F>,
    TEProjective<P>,
    Add,
    add,
    AddAssign,
    add_assign,
    |this: &'a AffineVar<P, F>, other: &'a AffineVar<P, F>| {

        if [this, other].is_constant() {
            assert!(this.is_constant() && other.is_constant());
            AffineVar::constant(this.value().unwrap() + &other.value().unwrap())
        } else {
            let cs = [this, other].cs();
            let a = P::COEFF_A;
            let d = P::COEFF_D;

            // Compute U = (x1 + y1) * (x2 + y2)
            let u1 = (&this.x * -a) + &this.y;
            let u2 = &other.x + &other.y;

            let u = u1 *  &u2;

            // Compute v0 = x1 * y2
            let v0 = &other.y * &this.x;

            // Compute v1 = x2 * y1
            let v1 = &other.x * &this.y;

            // Compute C = d*v0*v1
            let v2 = &v0 * &v1 * d;

            // Compute x3 = (v0 + v1) / (1 + v2)
            let x3 = F::new_witness(r1cs_core::ns!(cs, "x3"), || {
                let t0 = v0.value()? + &v1.value()?;
                let t1 = P::BaseField::one() + &v2.value()?;
                Ok(t0 * &t1.inverse().ok_or(SynthesisError::DivisionByZero)?)
            }).unwrap();

            let v2_plus_one = &v2 + P::BaseField::one();
            let v0_plus_v1 = &v0 + &v1;
            x3.mul_equals(&v2_plus_one, &v0_plus_v1).unwrap();

            // Compute y3 = (U + a * v0 - v1) / (1 - v2)
            let y3 = F::new_witness(r1cs_core::ns!(cs, "y3"), || {
                let t0 = u.value()? + &(a * &v0.value()?) - &v1.value()?;
                let t1 = P::BaseField::one() - &v2.value()?;
                Ok(t0 * &t1.inverse().ok_or(SynthesisError::DivisionByZero)?)
            }).unwrap();

            let one_minus_v2 = (&v2 - P::BaseField::one()).negate().unwrap();
            let a_v0 = &v0 * a;
            let u_plus_a_v0_minus_v1 = &u + &a_v0 - &v1;

            y3.mul_equals(&one_minus_v2, &u_plus_a_v0_minus_v1).unwrap();

            AffineVar::new(x3, y3)
        }
    },
    |this: &'a AffineVar<P, F>, other: TEProjective<P>| this + AffineVar::constant(other),
    (
        F :FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>
            + TwoBitLookupGadget<<P::BaseField as Field>::BasePrimeField, TableConstant = P::BaseField>,
        P: TEModelParameters,
    ),
    for <'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
);

impl_bounded_ops!(
    AffineVar<P, F>,
    TEProjective<P>,
    Sub,
    sub,
    SubAssign,
    sub_assign,
    |this: &'a AffineVar<P, F>, other: &'a AffineVar<P, F>| this + other.negate().unwrap(),
    |this: &'a AffineVar<P, F>, other: TEProjective<P>| this - AffineVar::constant(other),
    (
        F :FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>
            + TwoBitLookupGadget<<P::BaseField as Field>::BasePrimeField, TableConstant = P::BaseField>,
        P: TEModelParameters,
    ),
    for <'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>
);

impl<'a, P, F> GroupOpsBounds<'a, TEProjective<P>, AffineVar<P, F>> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>
        + TwoBitLookupGadget<<P::BaseField as Field>::BasePrimeField, TableConstant = P::BaseField>,
    for<'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
{
}

impl<'a, P, F> GroupOpsBounds<'a, TEProjective<P>, AffineVar<P, F>> for &'a AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>
        + TwoBitLookupGadget<<P::BaseField as Field>::BasePrimeField, TableConstant = P::BaseField>,
    for<'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
{
}

impl<P, F> CondSelectGadget<<P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
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

        Ok(Self::new(x, y))
    }
}

impl<P, F> EqGadget<<P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
{
    #[tracing::instrument(target = "r1cs")]
    fn is_eq(
        &self,
        other: &Self,
    ) -> Result<Boolean<<P::BaseField as Field>::BasePrimeField>, SynthesisError> {
        let x_equal = self.x.is_eq(&other.x)?;
        let y_equal = self.y.is_eq(&other.y)?;
        x_equal.and(&y_equal)
    }

    #[inline]
    #[tracing::instrument(target = "r1cs")]
    fn conditional_enforce_equal(
        &self,
        other: &Self,
        condition: &Boolean<<P::BaseField as Field>::BasePrimeField>,
    ) -> Result<(), SynthesisError> {
        self.x.conditional_enforce_equal(&other.x, condition)?;
        self.y.conditional_enforce_equal(&other.y, condition)?;
        Ok(())
    }

    #[inline]
    #[tracing::instrument(target = "r1cs")]
    fn conditional_enforce_not_equal(
        &self,
        other: &Self,
        condition: &Boolean<<P::BaseField as Field>::BasePrimeField>,
    ) -> Result<(), SynthesisError> {
        self.is_eq(other)?
            .and(condition)?
            .enforce_equal(&Boolean::Constant(false))
    }
}

impl<P, F> ToBitsGadget<<P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
{
    #[tracing::instrument(target = "r1cs")]
    fn to_bits_le(
        &self,
    ) -> Result<Vec<Boolean<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let mut x_bits = self.x.to_bits_le()?;
        let y_bits = self.y.to_bits_le()?;
        x_bits.extend_from_slice(&y_bits);
        Ok(x_bits)
    }

    #[tracing::instrument(target = "r1cs")]
    fn to_non_unique_bits_le(
        &self,
    ) -> Result<Vec<Boolean<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let mut x_bits = self.x.to_non_unique_bits_le()?;
        let y_bits = self.y.to_non_unique_bits_le()?;
        x_bits.extend_from_slice(&y_bits);

        Ok(x_bits)
    }
}

impl<P, F> ToBytesGadget<<P::BaseField as Field>::BasePrimeField> for AffineVar<P, F>
where
    P: TEModelParameters,
    F: FieldVar<P::BaseField, <P::BaseField as Field>::BasePrimeField>,
    for<'b> &'b F: FieldOpsBounds<'b, P::BaseField, F>,
{
    #[tracing::instrument(target = "r1cs")]
    fn to_bytes(
        &self,
    ) -> Result<Vec<UInt8<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let mut x_bytes = self.x.to_bytes()?;
        let y_bytes = self.y.to_bytes()?;
        x_bytes.extend_from_slice(&y_bytes);
        Ok(x_bytes)
    }

    #[tracing::instrument(target = "r1cs")]
    fn to_non_unique_bytes(
        &self,
    ) -> Result<Vec<UInt8<<P::BaseField as Field>::BasePrimeField>>, SynthesisError> {
        let mut x_bytes = self.x.to_non_unique_bytes()?;
        let y_bytes = self.y.to_non_unique_bytes()?;
        x_bytes.extend_from_slice(&y_bytes);

        Ok(x_bytes)
    }
}

#[cfg(test)]
#[allow(dead_code)]
pub(crate) fn test<P, GG>() -> Result<(), SynthesisError>
where
    P: TEModelParameters,
    GG: CurveVar<TEProjective<P>, <P::BaseField as Field>::BasePrimeField>,
    for<'a> &'a GG: GroupOpsBounds<'a, TEProjective<P>, GG>,
{
    use crate::prelude::*;
    use algebra::{test_rng, BitIteratorLE, Group, UniformRand};
    use r1cs_core::ConstraintSystem;

    crate::groups::test::group_test::<TEProjective<P>, _, GG>()?;

    let mut rng = test_rng();

    let cs = ConstraintSystem::<<P::BaseField as Field>::BasePrimeField>::new_ref();

    let a = TEProjective::<P>::rand(&mut rng);
    let b = TEProjective::<P>::rand(&mut rng);
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
    assert_eq!(cs.which_is_unsatisfied()?, None);

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
    let native_result = AffineCurve::mul(&aa.into_affine(), scalar);
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
