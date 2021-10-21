use algebra::{
    curves::{
        twisted_edwards_extended::GroupAffine as TEAffine, MontgomeryModelParameters,
        TEModelParameters,
    },
    BitIterator, Field,
};

use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::prelude::*;

use std::{borrow::Borrow, marker::PhantomData};

#[derive(Derivative)]
#[derivative(Debug, Clone)]
#[derivative(Debug(bound = "P: TEModelParameters, ConstraintF: Field"))]
#[must_use]
pub struct MontgomeryAffineGadget<
    P: TEModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
> {
    pub x: F,
    pub y: F,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
    #[derivative(Debug = "ignore")]
    _engine: PhantomData<ConstraintF>,
}

mod montgomery_affine_impl {
    use super::*;
    use crate::Assignment;
    use algebra::{twisted_edwards_extended::GroupAffine, AffineCurve, Field};
    use std::ops::{AddAssign, MulAssign, SubAssign};

    impl<P: TEModelParameters, ConstraintF: Field, F: FieldGadget<P::BaseField, ConstraintF>>
        MontgomeryAffineGadget<P, ConstraintF, F>
    {
        pub fn new(x: F, y: F) -> Self {
            Self {
                x,
                y,
                _params: PhantomData,
                _engine: PhantomData,
            }
        }

        pub fn from_edwards_to_coords(
            p: &TEAffine<P>,
        ) -> Result<(P::BaseField, P::BaseField), SynthesisError> {
            let montgomery_point: GroupAffine<P> = if p.y == P::BaseField::one() {
                GroupAffine::zero()
            } else {
                if p.x == P::BaseField::zero() {
                    GroupAffine::new(P::BaseField::zero(), P::BaseField::zero())
                } else {
                    let u = (P::BaseField::one() + &p.y)
                        * &(P::BaseField::one() - &p.y).inverse().unwrap();
                    let v = u * &p.x.inverse().unwrap();
                    GroupAffine::new(u, v)
                }
            };

            Ok((montgomery_point.x, montgomery_point.y))
        }

        pub fn from_edwards<CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            p: &TEAffine<P>,
        ) -> Result<Self, SynthesisError> {
            let montgomery_coords = Self::from_edwards_to_coords(p)?;

            let u = F::alloc(cs.ns(|| "u"), || Ok(montgomery_coords.0))?;

            let v = F::alloc(cs.ns(|| "v"), || Ok(montgomery_coords.1))?;

            Ok(Self::new(u, v))
        }

        pub fn into_edwards<CS: ConstraintSystem<ConstraintF>>(
            &self,
            mut cs: CS,
        ) -> Result<AffineGadget<P, ConstraintF, F>, SynthesisError> {
            // Compute u = x / y
            let u = F::alloc(cs.ns(|| "u"), || {
                let mut t0 = self.x.get_value().get()?;

                match self.y.get_value().get()?.inverse() {
                    Some(invy) => {
                        t0.mul_assign(&invy);

                        Ok(t0)
                    }
                    None => Err(SynthesisError::DivisionByZero),
                }
            })?;

            u.mul_equals(cs.ns(|| "u equals"), &self.y, &self.x)?;

            let v = F::alloc(cs.ns(|| "v"), || {
                let mut t0 = self.x.get_value().get()?;
                let mut t1 = t0.clone();
                t0.sub_assign(&P::BaseField::one());
                t1.add_assign(&P::BaseField::one());

                match t1.inverse() {
                    Some(t1) => {
                        t0.mul_assign(&t1);

                        Ok(t0)
                    }
                    None => Err(SynthesisError::DivisionByZero),
                }
            })?;

            let xplusone = self
                .x
                .add_constant(cs.ns(|| "x plus one"), &P::BaseField::one())?;
            let xminusone = self
                .x
                .sub_constant(cs.ns(|| "x minus one"), &P::BaseField::one())?;
            v.mul_equals(cs.ns(|| "v equals"), &xplusone, &xminusone)?;

            Ok(AffineGadget::new(u, v))
        }

        pub fn add<CS: ConstraintSystem<ConstraintF>>(
            &self,
            mut cs: CS,
            other: &Self,
        ) -> Result<Self, SynthesisError> {
            let lambda = F::alloc(cs.ns(|| "lambda"), || {
                let mut n = other.y.get_value().get()?;
                n.sub_assign(&self.y.get_value().get()?);

                let mut d = other.x.get_value().get()?;
                d.sub_assign(&self.x.get_value().get()?);

                match d.inverse() {
                    Some(d) => {
                        n.mul_assign(&d);
                        Ok(n)
                    }
                    None => Err(SynthesisError::DivisionByZero),
                }
            })?;
            let lambda_n = other.y.sub(cs.ns(|| "other.y - self.y"), &self.y)?;
            let lambda_d = other.x.sub(cs.ns(|| "other.x - self.x"), &self.x)?;
            lambda_d.mul_equals(cs.ns(|| "lambda equals"), &lambda, &lambda_n)?;

            // Compute x'' = B*lambda^2 - A - x - x'
            let xprime = F::alloc(cs.ns(|| "xprime"), || {
                Ok(
                    lambda.get_value().get()?.square() * &P::MontgomeryModelParameters::COEFF_B
                        - &P::MontgomeryModelParameters::COEFF_A
                        - &self.x.get_value().get()?
                        - &other.x.get_value().get()?,
                )
            })?;

            let xprime_lc = self
                .x
                .add(cs.ns(|| "self.x + other.x"), &other.x)?
                .add(cs.ns(|| "+ xprime"), &xprime)?
                .add_constant(cs.ns(|| "+ A"), &P::MontgomeryModelParameters::COEFF_A)?;
            // (lambda) * (lambda) = (A + x + x' + x'')
            let lambda_b = lambda.mul_by_constant(
                cs.ns(|| "lambda * b"),
                &P::MontgomeryModelParameters::COEFF_B,
            )?;
            lambda_b.mul_equals(cs.ns(|| "xprime equals"), &lambda, &xprime_lc)?;

            let yprime = F::alloc(cs.ns(|| "yprime"), || {
                Ok(-(self.y.get_value().get()?
                    + &(lambda.get_value().get()?
                        * &(xprime.get_value().get()? - &self.x.get_value().get()?))))
            })?;

            let xres = self.x.sub(cs.ns(|| "xres"), &xprime)?;
            let yres = self.y.add(cs.ns(|| "yres"), &yprime)?;
            lambda.mul_equals(cs.ns(|| "yprime equals"), &xres, &yres)?;
            Ok(MontgomeryAffineGadget::new(xprime, yprime))
        }
    }
}

#[derive(Derivative)]
#[derivative(Debug, Clone)]
#[derivative(Debug(bound = "P: TEModelParameters, ConstraintF: Field"))]
#[must_use]
pub struct AffineGadget<
    P: TEModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
> {
    pub x: F,
    pub y: F,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
    #[derivative(Debug = "ignore")]
    _engine: PhantomData<ConstraintF>,
}

impl<P: TEModelParameters, ConstraintF: Field, F: FieldGadget<P::BaseField, ConstraintF>>
    AffineGadget<P, ConstraintF, F>
{
    pub fn new(x: F, y: F) -> Self {
        Self {
            x,
            y,
            _params: PhantomData,
            _engine: PhantomData,
        }
    }
}

impl<P, ConstraintF, F> PartialEq for AffineGadget<P, ConstraintF, F>
where
    P: TEModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<P, ConstraintF, F> Eq for AffineGadget<P, ConstraintF, F>
where
    P: TEModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
}

mod affine_impl {
    use super::*;
    use crate::Assignment;
    use algebra::{curves::AffineCurve, Field, PrimeField};
    use std::ops::Neg;

    impl<P, ConstraintF, F> GroupGadget<TEAffine<P>, ConstraintF> for AffineGadget<P, ConstraintF, F>
    where
        P: TEModelParameters,
        ConstraintF: Field,
        F: FieldGadget<P::BaseField, ConstraintF>,
    {
        type Value = TEAffine<P>;
        type Variable = (F::Variable, F::Variable);

        #[inline]
        fn get_value(&self) -> Option<Self::Value> {
            match (self.x.get_value(), self.y.get_value()) {
                (Some(x), Some(y)) => Some(TEAffine::new(x, y)),
                (..) => None,
            }
        }

        #[inline]
        fn get_variable(&self) -> Self::Variable {
            (self.x.get_variable(), self.y.get_variable())
        }

        #[inline]
        fn zero<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
            Ok(Self::new(
                F::zero(cs.ns(|| "zero"))?,
                F::one(cs.ns(|| "one"))?,
            ))
        }

        //TODO: Implement this using enforce_verdict
        #[inline]
        fn is_zero<CS: ConstraintSystem<ConstraintF>>(
            &self,
            _: CS,
        ) -> Result<Boolean, SynthesisError> {
            unimplemented!()
        }

        /// Optimized constraints for checking Edwards point addition from ZCash
        /// developers Daira Hopwood and Sean Bowe. Requires only 6 constraints
        /// compared to 7 for the straightforward version we had earlier.
        fn add<CS: ConstraintSystem<ConstraintF>>(
            &self,
            mut cs: CS,
            other: &Self,
        ) -> Result<Self, SynthesisError> {
            let a = P::COEFF_A;
            let d = P::COEFF_D;

            // Compute U = (x1 + y1) * (x2 + y2)
            let u1 = self
                .x
                .mul_by_constant(cs.ns(|| "-A * x1"), &a.neg())?
                .add(cs.ns(|| "-A * x1 + y1"), &self.y)?;
            let u2 = other.x.add(cs.ns(|| "x2 + y2"), &other.y)?;

            let u = u1.mul(cs.ns(|| "(-A * x1 + y1) * (x2 + y2)"), &u2)?;

            // Compute v0 = x1 * y2
            let v0 = other.y.mul(&mut cs.ns(|| "v0"), &self.x)?;

            // Compute v1 = x2 * y1
            let v1 = other.x.mul(&mut cs.ns(|| "v1"), &self.y)?;

            // Compute C = d*v0*v1
            let v2 = v0
                .mul(cs.ns(|| "v0 * v1"), &v1)?
                .mul_by_constant(cs.ns(|| "D * v0 * v1"), &d)?;

            // Compute x3 = (v0 + v1) / (1 + v2)
            let x3 = F::alloc(&mut cs.ns(|| "x3"), || {
                let t0 = v0.get_value().get()? + &v1.get_value().get()?;
                let t1 = P::BaseField::one() + &v2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let one = P::BaseField::one();
            let v2_plus_one = v2.add_constant(cs.ns(|| "v2 + 1"), &one)?;
            let v0_plus_v1 = v0.add(cs.ns(|| "v0 + v1"), &v1)?;
            x3.mul_equals(cs.ns(|| "check x3"), &v2_plus_one, &v0_plus_v1)?;

            // Compute y3 = (U + a * v0 - v1) / (1 - v2)
            let y3 = F::alloc(&mut cs.ns(|| "y3"), || {
                let t0 =
                    u.get_value().get()? + &(a * &v0.get_value().get()?) - &v1.get_value().get()?;
                let t1 = P::BaseField::one() - &v2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let one_minus_v2 = v2
                .add_constant(cs.ns(|| "v2 - 1"), &(-one))?
                .negate(cs.ns(|| "1 - v2"))?;
            let a_v0 = v0.mul_by_constant(cs.ns(|| "a * v0"), &a)?;
            let u_plus_a_v0_minus_v1 = u
                .add(cs.ns(|| "u + a * v0"), &a_v0)?
                .sub(cs.ns(|| "u + a * v0 - v1"), &v1)?;

            y3.mul_equals(cs.ns(|| "check y3"), &one_minus_v2, &u_plus_a_v0_minus_v1)?;

            Ok(Self::new(x3, y3))
        }

        fn add_constant<CS: ConstraintSystem<ConstraintF>>(
            &self,
            mut cs: CS,
            other: &TEAffine<P>,
        ) -> Result<Self, SynthesisError> {
            let a = P::COEFF_A;
            let d = P::COEFF_D;
            let other_x = other.x;
            let other_y = other.y;

            // Compute U = (x1 + y1) * (x2 + y2)
            let u1 = self
                .x
                .mul_by_constant(cs.ns(|| "-A * x1"), &a.neg())?
                .add(cs.ns(|| "-A * x1 + y1"), &self.y)?;
            let u2 = other_x + &other_y;

            let u = u1.mul_by_constant(cs.ns(|| "(-A * x1 + y1) * (x2 + y2)"), &u2)?;

            // Compute v0 = x1 * y2
            let v0 = self.x.mul_by_constant(&mut cs.ns(|| "v0"), &other_y)?;

            // Compute v1 = x2 * y1
            let v1 = self.y.mul_by_constant(&mut cs.ns(|| "v1"), &other.x)?;

            // Compute C = d*v0*v1
            let v2 = v0
                .mul(cs.ns(|| "v0 * v1"), &v1)?
                .mul_by_constant(cs.ns(|| "D * v0 * v1"), &d)?;

            // Compute x3 = (v0 + v1) / (1 + v2)
            let x3 = F::alloc(&mut cs.ns(|| "x3"), || {
                let t0 = v0.get_value().get()? + &v1.get_value().get()?;
                let t1 = P::BaseField::one() + &v2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let one = P::BaseField::one();
            let v2_plus_one = v2.add_constant(cs.ns(|| "v2 + 1"), &one)?;
            let v0_plus_v1 = v0.add(cs.ns(|| "v0 + v1"), &v1)?;
            x3.mul_equals(cs.ns(|| "check x3"), &v2_plus_one, &v0_plus_v1)?;

            // Compute y3 = (U + a * v0 - v1) / (1 - v2)
            let y3 = F::alloc(&mut cs.ns(|| "y3"), || {
                let t0 =
                    u.get_value().get()? + &(a * &v0.get_value().get()?) - &v1.get_value().get()?;
                let t1 = P::BaseField::one() - &v2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let one_minus_v2 = v2
                .add_constant(cs.ns(|| "v2 - 1"), &(-one))?
                .negate(cs.ns(|| "1 - v2"))?;
            let a_v0 = v0.mul_by_constant(cs.ns(|| "a * v0"), &a)?;
            let u_plus_a_v0_minus_v1 = u
                .add(cs.ns(|| "u + a * v0"), &a_v0)?
                .sub(cs.ns(|| "u + a * v0 - v1"), &v1)?;

            y3.mul_equals(cs.ns(|| "check y3"), &one_minus_v2, &u_plus_a_v0_minus_v1)?;

            Ok(Self::new(x3, y3))
        }

        fn double_in_place<CS: ConstraintSystem<ConstraintF>>(
            &mut self,
            mut cs: CS,
        ) -> Result<(), SynthesisError> {
            let a = P::COEFF_A;

            // xy
            let xy = self.x.mul(cs.ns(|| "x * y"), &self.y)?;
            let x2 = self.x.square(cs.ns(|| "x * x"))?;
            let y2 = self.y.square(cs.ns(|| "y * y"))?;

            let a_x2 = x2.mul_by_constant(cs.ns(|| "a * x^2"), &a)?;

            // Compute x3 = (2xy) / (ax^2 + y^2)
            let x3 = F::alloc(&mut cs.ns(|| "x3"), || {
                let t0 = xy.get_value().get()?.double();
                let t1 = a * &x2.get_value().get()? + &y2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let a_x2_plus_y2 = a_x2.add(cs.ns(|| "v2 + 1"), &y2)?;
            let two_xy = xy.double(cs.ns(|| "2xy"))?;
            x3.mul_equals(cs.ns(|| "check x3"), &a_x2_plus_y2, &two_xy)?;

            // Compute y3 = (y^2 - ax^2) / (2 - ax^2 - y^2)
            let two = P::BaseField::one().double();
            let y3 = F::alloc(&mut cs.ns(|| "y3"), || {
                let a_x2 = a * &x2.get_value().get()?;
                let t0 = y2.get_value().get()? - &a_x2;
                let t1 = two - &a_x2 - &y2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;
            let y2_minus_a_x2 = y2.sub(cs.ns(|| "y^2 - ax^2"), &a_x2)?;
            let two_minus_ax2_minus_y2 = a_x2
                .add(cs.ns(|| "ax2 + y2"), &y2)?
                .negate(cs.ns(|| "-ax2 - y2"))?
                .add_constant(cs.ns(|| "2 -ax2 - y2"), &two)?;

            y3.mul_equals(
                cs.ns(|| "check y3"),
                &two_minus_ax2_minus_y2,
                &y2_minus_a_x2,
            )?;
            self.x = x3;
            self.y = y3;

            Ok(())
        }

        fn negate<CS: ConstraintSystem<ConstraintF>>(
            &self,
            mut cs: CS,
        ) -> Result<Self, SynthesisError> {
            Ok(Self::new(
                self.x.negate(cs.ns(|| "negate x"))?,
                self.y.clone(),
            ))
        }

        fn cost_of_add() -> usize {
            4 * F::cost_of_mul() + 2 * F::cost_of_mul_equals()
        }

        fn cost_of_double() -> usize {
            3 + 2 * F::cost_of_mul_equals()
        }
    }

    impl<P, ConstraintF, F> AllocGadget<TEAffine<P>, ConstraintF> for AffineGadget<P, ConstraintF, F>
    where
        P: TEModelParameters,
        ConstraintF: Field,
        F: FieldGadget<P::BaseField, ConstraintF>,
        Self: GroupGadget<TEAffine<P>, ConstraintF>,
    {
        fn alloc<FN, T, CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEAffine<P>>,
        {
            let (x, y) = match value_gen() {
                Ok(ge) => {
                    let ge = *ge.borrow();
                    (Ok(ge.x), Ok(ge.y))
                }
                _ => (
                    Err(SynthesisError::AssignmentMissing),
                    Err(SynthesisError::AssignmentMissing),
                ),
            };

            let d = P::COEFF_D;
            let a = P::COEFF_A;

            let x = F::alloc(&mut cs.ns(|| "x"), || x)?;
            let y = F::alloc(&mut cs.ns(|| "y"), || y)?;

            // Check that ax^2 + y^2 = 1 + dx^2y^2
            // We do this by checking that ax^2 - 1 = y^2 * (dx^2 - 1)
            let x2 = x.square(&mut cs.ns(|| "x^2"))?;
            let y2 = y.square(&mut cs.ns(|| "y^2"))?;

            let one = P::BaseField::one();
            let d_x2_minus_one = x2
                .mul_by_constant(cs.ns(|| "d * x^2"), &d)?
                .add_constant(cs.ns(|| "d * x^2 - 1"), &one.neg())?;

            let a_x2_minus_one = x2
                .mul_by_constant(cs.ns(|| "a * x^2"), &a)?
                .add_constant(cs.ns(|| "a * x^2 - 1"), &one.neg())?;

            d_x2_minus_one.mul_equals(cs.ns(|| "on curve check"), &y2, &a_x2_minus_one)?;
            Ok(Self::new(x, y))
        }

        #[inline]
        fn alloc_without_check<FN, T, CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEAffine<P>>,
        {
            let (x, y) = match value_gen() {
                Ok(ge) => {
                    let ge = *ge.borrow();
                    (Ok(ge.x), Ok(ge.y))
                }
                _ => (
                    Err(SynthesisError::AssignmentMissing),
                    Err(SynthesisError::AssignmentMissing),
                ),
            };

            let x = F::alloc(&mut cs.ns(|| "x"), || x)?;
            let y = F::alloc(&mut cs.ns(|| "y"), || y)?;

            Ok(Self::new(x, y))
        }

        #[inline]
        fn alloc_checked<FN, T, CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEAffine<P>>,
        {
            let alloc_and_prime_order_check = |mut cs: r1cs_core::Namespace<_, _>,
                                               value_gen: FN|
             -> Result<Self, SynthesisError> {
                let cofactor_weight = BitIterator::new(P::COFACTOR).filter(|b| *b).count();
                // If we multiply by r, we actually multiply by r - 2.
                let r_minus_1 = (-P::ScalarField::one()).into_repr();
                let r_weight = BitIterator::new(&r_minus_1).filter(|b| *b).count();

                // We pick the most efficient method of performing the prime order check:
                // If the cofactor has lower hamming weight than the scalar field's modulus,
                // we first multiply by the inverse of the cofactor, and then, after allocating,
                // multiply by the cofactor. This ensures the resulting point has no cofactors
                //
                // Else, we multiply by the scalar field's modulus and ensure that the result
                // is zero.
                if cofactor_weight < r_weight {
                    let ge = Self::alloc(cs.ns(|| "Alloc checked"), || {
                        value_gen().map(|ge| ge.borrow().mul_by_cofactor_inv())
                    })?;
                    let mut seen_one = false;
                    let mut result = Self::zero(cs.ns(|| "result"))?;
                    for (i, b) in BitIterator::new(P::COFACTOR).enumerate() {
                        let mut cs = cs.ns(|| format!("Iteration {}", i));

                        let old_seen_one = seen_one;
                        if seen_one {
                            result.double_in_place(cs.ns(|| "Double"))?;
                        } else {
                            seen_one = b;
                        }

                        if b {
                            result = if old_seen_one {
                                result.add(cs.ns(|| "Add"), &ge)?
                            } else {
                                ge.clone()
                            };
                        }
                    }
                    Ok(result)
                } else {
                    let ge = Self::alloc(cs.ns(|| "Alloc checked"), value_gen)?;
                    let mut seen_one = false;
                    let mut result = Self::zero(cs.ns(|| "result"))?;
                    // Returns bits in big-endian order
                    for (i, b) in BitIterator::new(r_minus_1).enumerate() {
                        let mut cs = cs.ns(|| format!("Iteration {}", i));

                        let old_seen_one = seen_one;
                        if seen_one {
                            result.double_in_place(cs.ns(|| "Double"))?;
                        } else {
                            seen_one = b;
                        }

                        if b {
                            result = if old_seen_one {
                                result.add(cs.ns(|| "Add"), &ge)?
                            } else {
                                ge.clone()
                            };
                        }
                    }
                    let neg_ge = ge.negate(cs.ns(|| "Negate ge"))?;
                    neg_ge.enforce_equal(cs.ns(|| "Check equals"), &result)?;
                    Ok(ge)
                }
            };

            let ge =
                alloc_and_prime_order_check(cs.ns(|| "alloc and prime order check"), value_gen)?;

            Ok(ge)
        }

        fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEAffine<P>>,
        {
            let (x, y) = match value_gen() {
                Ok(ge) => {
                    let ge = *ge.borrow();
                    (Ok(ge.x), Ok(ge.y))
                }
                _ => (
                    Err(SynthesisError::AssignmentMissing),
                    Err(SynthesisError::AssignmentMissing),
                ),
            };

            let d = P::COEFF_D;
            let a = P::COEFF_A;

            let x = F::alloc_input(&mut cs.ns(|| "x"), || x)?;
            let y = F::alloc_input(&mut cs.ns(|| "y"), || y)?;

            // Check that ax^2 + y^2 = 1 + dx^2y^2
            // We do this by checking that ax^2 - 1 = y^2 * (dx^2 - 1)
            let x2 = x.square(&mut cs.ns(|| "x^2"))?;
            let y2 = y.square(&mut cs.ns(|| "y^2"))?;

            let one = P::BaseField::one();
            let d_x2_minus_one = x2
                .mul_by_constant(cs.ns(|| "d * x^2"), &d)?
                .add_constant(cs.ns(|| "d * x^2 - 1"), &one.neg())?;

            let a_x2_minus_one = x2
                .mul_by_constant(cs.ns(|| "a * x^2"), &a)?
                .add_constant(cs.ns(|| "a * x^2 - 1"), &one.neg())?;

            d_x2_minus_one.mul_equals(cs.ns(|| "on curve check"), &y2, &a_x2_minus_one)?;
            Ok(Self::new(x, y))
        }
    }

    impl<P, ConstraintF, F> ConstantGadget<TEAffine<P>, ConstraintF> for AffineGadget<P, ConstraintF, F>
    where
        P: TEModelParameters,
        ConstraintF: Field,
        F: FieldGadget<P::BaseField, ConstraintF>,
        Self: GroupGadget<TEAffine<P>, ConstraintF>,
    {
        #[inline]
        fn from_value<CS: ConstraintSystem<ConstraintF>>(mut cs: CS, value: &TEAffine<P>) -> Self {
            let x = F::from_value(cs.ns(|| "hardcode x"), &value.x);
            let y = F::from_value(cs.ns(|| "hardcode y"), &value.y);

            Self::new(x, y)
        }

        #[inline]
        fn get_constant(&self) -> TEAffine<P> {
            let x = self.x.get_value().unwrap();
            let y = self.y.get_value().unwrap();

            TEAffine::<P>::new(x, y)
        }
    }
}

mod projective_impl {
    use super::*;
    use crate::Assignment;
    use algebra::{
        curves::twisted_edwards_extended::GroupProjective as TEProjective, AffineCurve, Field,
        PrimeField, ProjectiveCurve,
    };
    use std::ops::Neg;

    impl<P, ConstraintF, F> GroupGadget<TEProjective<P>, ConstraintF>
        for AffineGadget<P, ConstraintF, F>
    where
        P: TEModelParameters,
        ConstraintF: Field,
        F: FieldGadget<P::BaseField, ConstraintF>,
    {
        type Value = TEProjective<P>;
        type Variable = (F::Variable, F::Variable);

        #[inline]
        fn get_value(&self) -> Option<Self::Value> {
            match (self.x.get_value(), self.y.get_value()) {
                (Some(x), Some(y)) => Some(TEAffine::new(x, y).into()),
                (..) => None,
            }
        }

        #[inline]
        fn get_variable(&self) -> Self::Variable {
            (self.x.get_variable(), self.y.get_variable())
        }

        #[inline]
        fn zero<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
            Ok(Self::new(
                F::zero(cs.ns(|| "zero"))?,
                F::one(cs.ns(|| "one"))?,
            ))
        }

        //TODO: Implement this using enforce_verdict
        #[inline]
        fn is_zero<CS: ConstraintSystem<ConstraintF>>(
            &self,
            _: CS,
        ) -> Result<Boolean, SynthesisError> {
            unimplemented!()
        }

        /// Optimized constraints for checking Edwards point addition from ZCash
        /// developers Daira Hopwood and Sean Bowe. Requires only 6 constraints
        /// compared to 7 for the straightforward version we had earlier.
        fn add<CS: ConstraintSystem<ConstraintF>>(
            &self,
            mut cs: CS,
            other: &Self,
        ) -> Result<Self, SynthesisError> {
            let a = P::COEFF_A;
            let d = P::COEFF_D;

            // Compute U = (x1 + y1) * (x2 + y2)
            let u1 = self
                .x
                .mul_by_constant(cs.ns(|| "-A * x1"), &a.neg())?
                .add(cs.ns(|| "-A * x1 + y1"), &self.y)?;
            let u2 = other.x.add(cs.ns(|| "x2 + y2"), &other.y)?;

            let u = u1.mul(cs.ns(|| "(-A * x1 + y1) * (x2 + y2)"), &u2)?;

            // Compute v0 = x1 * y2
            let v0 = other.y.mul(&mut cs.ns(|| "v0"), &self.x)?;

            // Compute v1 = x2 * y1
            let v1 = other.x.mul(&mut cs.ns(|| "v1"), &self.y)?;

            // Compute C = d*v0*v1
            let v2 = v0
                .mul(cs.ns(|| "v0 * v1"), &v1)?
                .mul_by_constant(cs.ns(|| "D * v0 * v1"), &d)?;

            // Compute x3 = (v0 + v1) / (1 + v2)
            let x3 = F::alloc(&mut cs.ns(|| "x3"), || {
                let t0 = v0.get_value().get()? + &v1.get_value().get()?;
                let t1 = P::BaseField::one() + &v2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let one = P::BaseField::one();
            let v2_plus_one = v2.add_constant(cs.ns(|| "v2 + 1"), &one)?;
            let v0_plus_v1 = v0.add(cs.ns(|| "v0 + v1"), &v1)?;
            x3.mul_equals(cs.ns(|| "check x3"), &v2_plus_one, &v0_plus_v1)?;

            // Compute y3 = (U + a * v0 - v1) / (1 - v2)
            let y3 = F::alloc(&mut cs.ns(|| "y3"), || {
                let t0 =
                    u.get_value().get()? + &(a * &v0.get_value().get()?) - &v1.get_value().get()?;
                let t1 = P::BaseField::one() - &v2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let one_minus_v2 = v2
                .add_constant(cs.ns(|| "v2 - 1"), &(-one))?
                .negate(cs.ns(|| "1 - v2"))?;
            let a_v0 = v0.mul_by_constant(cs.ns(|| "a * v0"), &a)?;
            let u_plus_a_v0_minus_v1 = u
                .add(cs.ns(|| "u + a * v0"), &a_v0)?
                .sub(cs.ns(|| "u + a * v0 - v1"), &v1)?;

            y3.mul_equals(cs.ns(|| "check y3"), &one_minus_v2, &u_plus_a_v0_minus_v1)?;

            Ok(Self::new(x3, y3))
        }

        fn add_constant<CS: ConstraintSystem<ConstraintF>>(
            &self,
            mut cs: CS,
            other: &TEProjective<P>,
        ) -> Result<Self, SynthesisError> {
            let a = P::COEFF_A;
            let d = P::COEFF_D;
            let other = other.into_affine();
            let other_x = other.x;
            let other_y = other.y;

            // Compute U = (x1 + y1) * (x2 + y2)
            let u1 = self
                .x
                .mul_by_constant(cs.ns(|| "-A * x1"), &a.neg())?
                .add(cs.ns(|| "-A * x1 + y1"), &self.y)?;
            let u2 = other_x + &other_y;

            let u = u1.mul_by_constant(cs.ns(|| "(-A * x1 + y1) * (x2 + y2)"), &u2)?;

            // Compute v0 = x1 * y2
            let v0 = self.x.mul_by_constant(&mut cs.ns(|| "v0"), &other_y)?;

            // Compute v1 = x2 * y1
            let v1 = self.y.mul_by_constant(&mut cs.ns(|| "v1"), &other.x)?;

            // Compute C = d*v0*v1
            let v2 = v0
                .mul(cs.ns(|| "v0 * v1"), &v1)?
                .mul_by_constant(cs.ns(|| "D * v0 * v1"), &d)?;

            // Compute x3 = (v0 + v1) / (1 + v2)
            let x3 = F::alloc(&mut cs.ns(|| "x3"), || {
                let t0 = v0.get_value().get()? + &v1.get_value().get()?;
                let t1 = P::BaseField::one() + &v2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let one = P::BaseField::one();
            let v2_plus_one = v2.add_constant(cs.ns(|| "v2 + 1"), &one)?;
            let v0_plus_v1 = v0.add(cs.ns(|| "v0 + v1"), &v1)?;
            x3.mul_equals(cs.ns(|| "check x3"), &v2_plus_one, &v0_plus_v1)?;

            // Compute y3 = (U + a * v0 - v1) / (1 - v2)
            let y3 = F::alloc(&mut cs.ns(|| "y3"), || {
                let t0 =
                    u.get_value().get()? + &(a * &v0.get_value().get()?) - &v1.get_value().get()?;
                let t1 = P::BaseField::one() - &v2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let one_minus_v2 = v2
                .add_constant(cs.ns(|| "v2 - 1"), &(-one))?
                .negate(cs.ns(|| "1 - v2"))?;
            let a_v0 = v0.mul_by_constant(cs.ns(|| "a * v0"), &a)?;
            let u_plus_a_v0_minus_v1 = u
                .add(cs.ns(|| "u + a * v0"), &a_v0)?
                .sub(cs.ns(|| "u + a * v0 - v1"), &v1)?;

            y3.mul_equals(cs.ns(|| "check y3"), &one_minus_v2, &u_plus_a_v0_minus_v1)?;

            Ok(Self::new(x3, y3))
        }

        fn double_in_place<CS: ConstraintSystem<ConstraintF>>(
            &mut self,
            mut cs: CS,
        ) -> Result<(), SynthesisError> {
            let a = P::COEFF_A;

            // xy
            let xy = self.x.mul(cs.ns(|| "x * y"), &self.y)?;
            let x2 = self.x.square(cs.ns(|| "x * x"))?;
            let y2 = self.y.square(cs.ns(|| "y * y"))?;

            let a_x2 = x2.mul_by_constant(cs.ns(|| "a * x^2"), &a)?;

            // Compute x3 = (2xy) / (ax^2 + y^2)
            let x3 = F::alloc(&mut cs.ns(|| "x3"), || {
                let t0 = xy.get_value().get()?.double();
                let t1 = a * &x2.get_value().get()? + &y2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;

            let a_x2_plus_y2 = a_x2.add(cs.ns(|| "v2 + 1"), &y2)?;
            let two_xy = xy.double(cs.ns(|| "2xy"))?;
            x3.mul_equals(cs.ns(|| "check x3"), &a_x2_plus_y2, &two_xy)?;

            // Compute y3 = (y^2 - ax^2) / (2 - ax^2 - y^2)
            let two = P::BaseField::one().double();
            let y3 = F::alloc(&mut cs.ns(|| "y3"), || {
                let a_x2 = a * &x2.get_value().get()?;
                let t0 = y2.get_value().get()? - &a_x2;
                let t1 = two - &a_x2 - &y2.get_value().get()?;
                Ok(t0 * &t1.inverse().get()?)
            })?;
            let y2_minus_a_x2 = y2.sub(cs.ns(|| "y^2 - ax^2"), &a_x2)?;
            let two_minus_ax2_minus_y2 = a_x2
                .add(cs.ns(|| "ax2 + y2"), &y2)?
                .negate(cs.ns(|| "-ax2 - y2"))?
                .add_constant(cs.ns(|| "2 -ax2 - y2"), &two)?;

            y3.mul_equals(
                cs.ns(|| "check y3"),
                &two_minus_ax2_minus_y2,
                &y2_minus_a_x2,
            )?;
            self.x = x3;
            self.y = y3;

            Ok(())
        }

        fn negate<CS: ConstraintSystem<ConstraintF>>(
            &self,
            mut cs: CS,
        ) -> Result<Self, SynthesisError> {
            Ok(Self::new(
                self.x.negate(cs.ns(|| "negate x"))?,
                self.y.clone(),
            ))
        }

        fn precomputed_base_scalar_mul<'a, CS, I, B>(
            &mut self,
            mut cs: CS,
            scalar_bits_with_base_powers: I,
        ) -> Result<(), SynthesisError>
        where
            CS: ConstraintSystem<ConstraintF>,
            I: Iterator<Item = (B, &'a TEProjective<P>)>,
            B: Borrow<Boolean>,
        {
            let scalar_bits_with_base_powers: Vec<_> = scalar_bits_with_base_powers
                .map(|(bit, base)| (bit.borrow().clone(), base.clone()))
                .collect();
            let zero = TEProjective::zero();
            for (i, bits_base_powers) in scalar_bits_with_base_powers.chunks(2).enumerate() {
                let mut cs = cs.ns(|| format!("Chunk {}", i));
                if bits_base_powers.len() == 2 {
                    let bits = [bits_base_powers[0].0, bits_base_powers[1].0];
                    let base_powers = [bits_base_powers[0].1, bits_base_powers[1].1];

                    let mut table = [
                        zero,
                        base_powers[0],
                        base_powers[1],
                        base_powers[0] + &base_powers[1],
                    ];

                    TEProjective::batch_normalization(&mut table);
                    let x_s = [table[0].x, table[1].x, table[2].x, table[3].x];
                    let y_s = [table[0].y, table[1].y, table[2].y, table[3].y];

                    let x: F = F::two_bit_lookup(cs.ns(|| "Lookup x"), &bits, &x_s)?;
                    let y: F = F::two_bit_lookup(cs.ns(|| "Lookup y"), &bits, &y_s)?;
                    let adder: Self = Self::new(x, y);
                    *self = <Self as GroupGadget<TEProjective<P>, ConstraintF>>::add(
                        self,
                        &mut cs.ns(|| "Add"),
                        &adder,
                    )?;
                } else if bits_base_powers.len() == 1 {
                    let bit = bits_base_powers[0].0;
                    let base_power = bits_base_powers[0].1;
                    let new_encoded =
                        self.add_constant(&mut cs.ns(|| "Add base power"), &base_power)?;
                    *self = Self::conditionally_select(
                        &mut cs.ns(|| "Conditional Select"),
                        &bit,
                        &new_encoded,
                        &self,
                    )?;
                }
            }

            Ok(())
        }

        fn precomputed_base_3_bit_signed_digit_scalar_mul<'a, CS, I, J, B>(
            mut cs: CS,
            bases: &[B],
            scalars: &[J],
        ) -> Result<Self, SynthesisError>
        where
            CS: ConstraintSystem<ConstraintF>,
            I: Borrow<[Boolean]>,
            J: Borrow<[I]>,
            B: Borrow<[TEProjective<P>]>,
        {
            const CHUNK_SIZE: usize = 3;
            let mut edwards_result: Option<AffineGadget<P, ConstraintF, F>> = None;
            let mut result: Option<MontgomeryAffineGadget<P, ConstraintF, F>> = None;

            let mut process_segment_result =
                |mut cs: r1cs_core::Namespace<_, _>,
                 result: &MontgomeryAffineGadget<P, ConstraintF, F>|
                 -> Result<(), SynthesisError> {
                    let segment_result = result.into_edwards(cs.ns(|| "segment result"))?;
                    match edwards_result {
                        None => {
                            edwards_result = Some(segment_result);
                        }
                        Some(ref mut edwards_result) => {
                            *edwards_result = GroupGadget::<TEAffine<P>, ConstraintF>::add(
                                &segment_result,
                                cs.ns(|| "edwards addition"),
                                edwards_result,
                            )?;
                        }
                    }

                    Ok(())
                };

            // Compute ∏(h_i^{m_i}) for all i.
            for (segment_i, (segment_bits_chunks, segment_powers)) in
                scalars.into_iter().zip(bases.iter()).enumerate()
            {
                for (i, (bits, base_power)) in segment_bits_chunks
                    .borrow()
                    .into_iter()
                    .zip(segment_powers.borrow().iter())
                    .enumerate()
                {
                    let base_power = base_power.borrow();
                    let mut acc_power = *base_power;
                    let mut coords = vec![];
                    for _ in 0..4 {
                        coords.push(acc_power);
                        acc_power = acc_power + base_power;
                    }

                    let bits = bits.borrow().to_bits(
                        &mut cs.ns(|| format!("Convert Scalar {}, {} to bits", segment_i, i)),
                    )?;
                    if bits.len() != CHUNK_SIZE {
                        return Err(SynthesisError::Unsatisfiable);
                    }

                    let coords = coords
                        .iter()
                        .map(|p| {
                            let p = p.into_affine();
                            MontgomeryAffineGadget::<P, ConstraintF, F>::from_edwards_to_coords(&p)
                                .unwrap()
                        })
                        .collect::<Vec<_>>();

                    let x_coeffs = coords.iter().map(|p| p.0).collect::<Vec<_>>();
                    let y_coeffs = coords.iter().map(|p| p.1).collect::<Vec<_>>();

                    let precomp = Boolean::and(
                        cs.ns(|| format!("precomp in window {}, {}", segment_i, i)),
                        &bits[0],
                        &bits[1],
                    )?;

                    let x = F::two_bit_lookup_lc(
                        cs.ns(|| format!("x in window {}, {}", segment_i, i)),
                        &precomp,
                        &[bits[0], bits[1]],
                        &x_coeffs,
                    )?;

                    let y = F::three_bit_cond_neg_lookup(
                        cs.ns(|| format!("y lookup in window {}, {}", segment_i, i)),
                        &bits,
                        &precomp,
                        &y_coeffs,
                    )?;

                    let tmp = MontgomeryAffineGadget::new(x, y);

                    match result {
                        None => {
                            result = Some(tmp);
                        }
                        Some(ref mut result) => {
                            *result = tmp.add(
                                cs.ns(|| format!("addition of window {}, {}", segment_i, i)),
                                result,
                            )?;
                        }
                    }
                }

                process_segment_result(
                    cs.ns(|| format!("window {}", segment_i)),
                    &result.unwrap(),
                )?;
                result = None;
            }
            if result.is_some() {
                process_segment_result(cs.ns(|| "leftover"), &result.unwrap())?;
            }
            Ok(edwards_result.unwrap())
        }

        fn cost_of_add() -> usize {
            4 * F::cost_of_mul() + 2 * F::cost_of_mul_equals()
        }

        fn cost_of_double() -> usize {
            3 + 2 * F::cost_of_mul_equals()
        }
    }

    impl<P, ConstraintF, F> AllocGadget<TEProjective<P>, ConstraintF>
        for AffineGadget<P, ConstraintF, F>
    where
        P: TEModelParameters,
        ConstraintF: Field,
        F: FieldGadget<P::BaseField, ConstraintF>,
        Self: GroupGadget<TEProjective<P>, ConstraintF>,
    {
        fn alloc<FN, T, CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEProjective<P>>,
        {
            let (x, y) = match value_gen() {
                Ok(ge) => {
                    let ge = ge.borrow().into_affine();
                    (Ok(ge.x), Ok(ge.y))
                }
                _ => (
                    Err(SynthesisError::AssignmentMissing),
                    Err(SynthesisError::AssignmentMissing),
                ),
            };

            let d = P::COEFF_D;
            let a = P::COEFF_A;

            let x = F::alloc(&mut cs.ns(|| "x"), || x)?;
            let y = F::alloc(&mut cs.ns(|| "y"), || y)?;

            // Check that ax^2 + y^2 = 1 + dx^2y^2
            // We do this by checking that ax^2 - 1 = y^2 * (dx^2 - 1)
            let x2 = x.square(&mut cs.ns(|| "x^2"))?;
            let y2 = y.square(&mut cs.ns(|| "y^2"))?;

            let one = P::BaseField::one();
            let d_x2_minus_one = x2
                .mul_by_constant(cs.ns(|| "d * x^2"), &d)?
                .add_constant(cs.ns(|| "d * x^2 - 1"), &one.neg())?;

            let a_x2_minus_one = x2
                .mul_by_constant(cs.ns(|| "a * x^2"), &a)?
                .add_constant(cs.ns(|| "a * x^2 - 1"), &one.neg())?;

            d_x2_minus_one.mul_equals(cs.ns(|| "on curve check"), &y2, &a_x2_minus_one)?;
            Ok(Self::new(x, y))
        }

        #[inline]
        fn alloc_without_check<FN, T, CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEProjective<P>>,
        {
            let (x, y) = match value_gen() {
                Ok(ge) => {
                    let ge = ge.borrow().into_affine();
                    (Ok(ge.x), Ok(ge.y))
                }
                _ => (
                    Err(SynthesisError::AssignmentMissing),
                    Err(SynthesisError::AssignmentMissing),
                ),
            };

            let x = F::alloc(&mut cs.ns(|| "x"), || x)?;
            let y = F::alloc(&mut cs.ns(|| "y"), || y)?;

            Ok(Self::new(x, y))
        }

        #[inline]
        fn alloc_checked<FN, T, CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEProjective<P>>,
        {
            let alloc_and_prime_order_check = |mut cs: r1cs_core::Namespace<_, _>,
                                               value_gen: FN|
             -> Result<Self, SynthesisError> {
                let cofactor_weight = BitIterator::new(P::COFACTOR).filter(|b| *b).count();
                // If we multiply by r, we actually multiply by r - 2.
                let r_minus_1 = (-P::ScalarField::one()).into_repr();
                let r_weight = BitIterator::new(&r_minus_1).filter(|b| *b).count();

                // We pick the most efficient method of performing the prime order check:
                // If the cofactor has lower hamming weight than the scalar field's modulus,
                // we first multiply by the inverse of the cofactor, and then, after allocating,
                // multiply by the cofactor. This ensures the resulting point has no cofactors
                //
                // Else, we multiply by the scalar field's modulus and ensure that the result
                // is zero.
                if cofactor_weight < r_weight {
                    let ge = Self::alloc(cs.ns(|| "Alloc checked"), || {
                        value_gen().map(|ge| {
                            ge.borrow()
                                .into_affine()
                                .mul_by_cofactor_inv()
                                .into_projective()
                        })
                    })?;
                    let mut seen_one = false;
                    let mut result = Self::zero(cs.ns(|| "result"))?;
                    for (i, b) in BitIterator::new(P::COFACTOR).enumerate() {
                        let mut cs = cs.ns(|| format!("Iteration {}", i));

                        let old_seen_one = seen_one;
                        if seen_one {
                            result.double_in_place(cs.ns(|| "Double"))?;
                        } else {
                            seen_one = b;
                        }

                        if b {
                            result = if old_seen_one {
                                result.add(cs.ns(|| "Add"), &ge)?
                            } else {
                                ge.clone()
                            };
                        }
                    }
                    Ok(result)
                } else {
                    let ge = Self::alloc(cs.ns(|| "Alloc checked"), value_gen)?;
                    let mut seen_one = false;
                    let mut result = Self::zero(cs.ns(|| "result"))?;
                    // Returns bits in big-endian order
                    for (i, b) in BitIterator::new(r_minus_1).enumerate() {
                        let mut cs = cs.ns(|| format!("Iteration {}", i));

                        let old_seen_one = seen_one;
                        if seen_one {
                            result.double_in_place(cs.ns(|| "Double"))?;
                        } else {
                            seen_one = b;
                        }

                        if b {
                            result = if old_seen_one {
                                result.add(cs.ns(|| "Add"), &ge)?
                            } else {
                                ge.clone()
                            };
                        }
                    }
                    let neg_ge = ge.negate(cs.ns(|| "Negate ge"))?;
                    neg_ge.enforce_equal(cs.ns(|| "Check equals"), &result)?;
                    Ok(ge)
                }
            };

            let ge =
                alloc_and_prime_order_check(cs.ns(|| "alloc and prime order check"), value_gen)?;

            Ok(ge)
        }

        fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEProjective<P>>,
        {
            let (x, y) = match value_gen() {
                Ok(ge) => {
                    let ge = ge.borrow().into_affine();
                    (Ok(ge.x), Ok(ge.y))
                }
                _ => (
                    Err(SynthesisError::AssignmentMissing),
                    Err(SynthesisError::AssignmentMissing),
                ),
            };

            let d = P::COEFF_D;
            let a = P::COEFF_A;

            let x = F::alloc_input(&mut cs.ns(|| "x"), || x)?;
            let y = F::alloc_input(&mut cs.ns(|| "y"), || y)?;

            // Check that ax^2 + y^2 = 1 + dx^2y^2
            // We do this by checking that ax^2 - 1 = y^2 * (dx^2 - 1)
            let x2 = x.square(&mut cs.ns(|| "x^2"))?;
            let y2 = y.square(&mut cs.ns(|| "y^2"))?;

            let one = P::BaseField::one();
            let d_x2_minus_one = x2
                .mul_by_constant(cs.ns(|| "d * x^2"), &d)?
                .add_constant(cs.ns(|| "d * x^2 - 1"), &one.neg())?;

            let a_x2_minus_one = x2
                .mul_by_constant(cs.ns(|| "a * x^2"), &a)?
                .add_constant(cs.ns(|| "a * x^2 - 1"), &one.neg())?;

            d_x2_minus_one.mul_equals(cs.ns(|| "on curve check"), &y2, &a_x2_minus_one)?;
            Ok(Self::new(x, y))
        }
    }

    impl<P, ConstraintF, F> ConstantGadget<TEProjective<P>, ConstraintF>
        for AffineGadget<P, ConstraintF, F>
    where
        P: TEModelParameters,
        ConstraintF: Field,
        F: FieldGadget<P::BaseField, ConstraintF>,
        Self: GroupGadget<TEProjective<P>, ConstraintF>,
    {
        #[inline]
        fn from_value<CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            value: &TEProjective<P>,
        ) -> Self {
            let value = value.into_affine();
            let x = F::from_value(cs.ns(|| "hardcode x"), &value.x);
            let y = F::from_value(cs.ns(|| "hardcode y"), &value.y);

            Self::new(x, y)
        }

        #[inline]
        fn get_constant(&self) -> TEProjective<P> {
            let value_proj =
                TEAffine::<P>::new(self.x.get_value().unwrap(), self.y.get_value().unwrap())
                    .into_projective();
            let x = value_proj.x;
            let y = value_proj.y;
            let t = value_proj.t;
            let z = value_proj.z;
            TEProjective::<P>::new(x, y, t, z)
        }
    }
}

impl<P, ConstraintF, F> CondSelectGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: TEModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError> {
        let x = F::conditionally_select(&mut cs.ns(|| "x"), cond, &first.x, &second.x)?;
        let y = F::conditionally_select(&mut cs.ns(|| "y"), cond, &first.y, &second.y)?;

        Ok(Self::new(x, y))
    }

    fn cost() -> usize {
        2 * <F as CondSelectGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF, F> EqGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: TEModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        let b0 = self.x.is_eq(cs.ns(|| "x"), &other.x)?;
        let b1 = self.y.is_eq(cs.ns(|| "y"), &other.y)?;
        Boolean::and(cs.ns(|| "x AND y"), &b0, &b1)
    }

    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.x
            .conditional_enforce_equal(cs.ns(|| "x"), &other.x, should_enforce)?;
        self.y
            .conditional_enforce_equal(cs.ns(|| "y"), &other.y, should_enforce)?;
        Ok(())
    }

    #[inline]
    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        let is_equal = self.is_eq(cs.ns(|| "is_eq(self, other)"), other)?;
        Boolean::and(
            cs.ns(|| "is_equal AND should_enforce"),
            &is_equal,
            should_enforce,
        )?
        .enforce_equal(
            cs.ns(|| "is_equal AND should_enforce == false"),
            &Boolean::Constant(false),
        )
    }
}

impl<P, ConstraintF, F> ToBitsGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: TEModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut x_bits = self.x.to_bits(cs.ns(|| "X Coordinate To Bits"))?;
        let y_bits = self.y.to_bits(cs.ns(|| "Y Coordinate To Bits"))?;
        x_bits.extend_from_slice(&y_bits);
        Ok(x_bits)
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut x_bits = self.x.to_bits_strict(cs.ns(|| "X Coordinate To Bits"))?;
        let y_bits = self.y.to_bits_strict(cs.ns(|| "Y Coordinate To Bits"))?;
        x_bits.extend_from_slice(&y_bits);

        Ok(x_bits)
    }
}

impl<P, ConstraintF, F> ToBytesGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: TEModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x_bytes = self.x.to_bytes(cs.ns(|| "x"))?;
        let y_bytes = self.y.to_bytes(cs.ns(|| "y"))?;
        x_bytes.extend_from_slice(&y_bytes);
        Ok(x_bytes)
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x_bytes = self.x.to_bytes_strict(cs.ns(|| "x"))?;
        let y_bytes = self.y.to_bytes_strict(cs.ns(|| "y"))?;
        x_bytes.extend_from_slice(&y_bytes);

        Ok(x_bytes)
    }
}

#[cfg(test)]
#[allow(dead_code)]
pub(crate) fn test<ConstraintF, P, GG>()
where
    ConstraintF: Field,
    P: TEModelParameters,
    GG: GroupGadget<TEAffine<P>, ConstraintF, Value = TEAffine<P>>,
{
    use crate::{
        boolean::AllocatedBit, groups::test::group_test, prelude::*,
        test_constraint_system::TestConstraintSystem,
    };
    use algebra::{Group, PrimeField, UniformRand};
    use rand::{thread_rng, Rng};

    group_test::<ConstraintF, TEAffine<P>, GG>();

    let mut cs = TestConstraintSystem::new();

    let a: TEAffine<P> = UniformRand::rand(&mut thread_rng());
    let gadget_a = GG::alloc(&mut cs.ns(|| "a"), || Ok(a)).unwrap();
    // Check mul_bits
    let scalar: <TEAffine<P> as Group>::ScalarField = UniformRand::rand(&mut thread_rng());
    let native_result = a.mul(&scalar);

    let mut scalar: Vec<bool> = BitIterator::new(scalar.into_repr()).collect();
    // Get the scalar bits into little-endian form.
    scalar.reverse();
    let input = Vec::<Boolean>::alloc(cs.ns(|| "Input"), || Ok(scalar)).unwrap();
    let zero = GG::zero(cs.ns(|| "zero")).unwrap();
    let result = gadget_a
        .mul_bits(cs.ns(|| "mul_bits"), &zero, input.iter())
        .unwrap();
    let gadget_value = result.get_value().expect("Gadget_result failed");
    assert_eq!(native_result, gadget_value);

    assert!(cs.is_satisfied());

    // Test the cost of allocation, conditional selection, and point addition.
    let mut cs = TestConstraintSystem::new();

    let bit = AllocatedBit::alloc(&mut cs.ns(|| "bool"), || Ok(true))
        .unwrap()
        .into();

    let mut rng = thread_rng();
    let a: TEAffine<P> = rng.gen();
    let b: TEAffine<P> = rng.gen();
    let gadget_a = GG::alloc(&mut cs.ns(|| "a"), || Ok(a)).unwrap();
    let gadget_b = GG::alloc(&mut cs.ns(|| "b"), || Ok(b)).unwrap();
    let alloc_cost = cs.num_constraints();
    let _ =
        GG::conditionally_select(&mut cs.ns(|| "cond_select"), &bit, &gadget_a, &gadget_b).unwrap();
    let cond_select_cost = cs.num_constraints() - alloc_cost;

    let _ = gadget_a.add(&mut cs.ns(|| "ab"), &gadget_b).unwrap();
    let add_cost = cs.num_constraints() - cond_select_cost - alloc_cost;
    assert_eq!(cond_select_cost, <GG as CondSelectGadget<_>>::cost());
    assert_eq!(add_cost, GG::cost_of_add());
    assert!(cs.is_satisfied());
}
