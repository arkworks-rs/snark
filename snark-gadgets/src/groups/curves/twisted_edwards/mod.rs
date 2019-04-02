use algebra::{
    curves::{twisted_edwards_extended::GroupAffine as TEAffine, TEModelParameters},
    BitIterator, PairingEngine,
};

use snark::{ConstraintSystem, SynthesisError};

use crate::{
    boolean::Boolean,
    fields::FieldGadget,
    groups::GroupGadget,
    uint8::UInt8,
    utils::{
        AllocGadget, CondSelectGadget, ConditionalEqGadget, EqGadget, NEqGadget, ToBitsGadget,
        ToBytesGadget,
    },
};

use std::{borrow::Borrow, marker::PhantomData};

pub mod edwards_bls12;
pub mod edwards_sw6;
pub mod jubjub;
#[cfg(test)]
mod test;

#[derive(Derivative)]
#[derivative(Debug, Clone)]
#[must_use]
pub struct AffineGadget<P: TEModelParameters, E: PairingEngine, F: FieldGadget<P::BaseField, E>> {
    pub x:   F,
    pub y:   F,
    _params: PhantomData<P>,
    _engine: PhantomData<E>,
}

impl<P: TEModelParameters, E: PairingEngine, F: FieldGadget<P::BaseField, E>>
    AffineGadget<P, E, F>
{
    pub fn new(x: F, y: F) -> Self {
        Self {
            x,
            y,
            _params: PhantomData,
            _engine: PhantomData,
        }
    }

    pub fn alloc_without_check<FN, CS: ConstraintSystem<E>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<TEAffine<P>, SynthesisError>,
    {
        let (x, y) = match value_gen() {
            Ok(fe) => (Ok(fe.x), Ok(fe.y)),
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let x = F::alloc(&mut cs.ns(|| "x"), || x)?;
        let y = F::alloc(&mut cs.ns(|| "y"), || y)?;

        Ok(Self::new(x, y))
    }
}

impl<P, E, F> PartialEq for AffineGadget<P, E, F>
where
    P: TEModelParameters,
    E: PairingEngine,
    F: FieldGadget<P::BaseField, E>,
{
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<P, E, F> Eq for AffineGadget<P, E, F>
where
    P: TEModelParameters,
    E: PairingEngine,
    F: FieldGadget<P::BaseField, E>,
{
}

mod affine_impl {
    use super::*;
    use crate::Assignment;
    use algebra::{curves::AffineCurve, Field, PrimeField};
    use std::ops::Neg;

    impl<P, E, F> GroupGadget<TEAffine<P>, E> for AffineGadget<P, E, F>
    where
        P: TEModelParameters,
        E: PairingEngine,
        F: FieldGadget<P::BaseField, E>,
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
        fn zero<CS: ConstraintSystem<E>>(mut cs: CS) -> Result<Self, SynthesisError> {
            Ok(Self::new(
                F::zero(cs.ns(|| "zero"))?,
                F::one(cs.ns(|| "one"))?,
            ))
        }

        /// Optimized constraints for checking Edwards point addition from ZCash
        /// developers Daira Hopwood and Sean Bowe. Requires only 6 constraints
        /// compared to 7 for the straightforward version we had earlier.
        fn add<CS: ConstraintSystem<E>>(
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

        fn add_constant<CS: ConstraintSystem<E>>(
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

        fn double_in_place<CS: ConstraintSystem<E>>(
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

        fn negate<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Self, SynthesisError> {
            Ok(Self::new(
                self.x.negate(cs.ns(|| "negate x"))?,
                self.y.clone(),
            ))
        }

        fn cost_of_add() -> usize {
            4 + 2 * F::cost_of_mul()
        }

        fn cost_of_double() -> usize {
            4 + F::cost_of_mul()
        }
    }

    impl<P, E, F> AllocGadget<TEAffine<P>, E> for AffineGadget<P, E, F>
    where
        P: TEModelParameters,
        E: PairingEngine,
        F: FieldGadget<P::BaseField, E>,
        Self: GroupGadget<TEAffine<P>, E>,
    {
        fn alloc<FN, T, CS: ConstraintSystem<E>>(
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
                },
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

        fn alloc_checked<FN, T, CS: ConstraintSystem<E>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEAffine<P>>,
        {
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
        }

        fn alloc_input<FN, T, CS: ConstraintSystem<E>>(
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
                },
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
}

mod projective_impl {
    use super::*;
    use crate::Assignment;
    use algebra::{
        curves::twisted_edwards_extended::GroupProjective as TEProjective, AffineCurve, Field,
        PrimeField, ProjectiveCurve,
    };
    use std::ops::Neg;

    impl<P, E, F> GroupGadget<TEProjective<P>, E> for AffineGadget<P, E, F>
    where
        P: TEModelParameters,
        E: PairingEngine,
        F: FieldGadget<P::BaseField, E>,
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
        fn zero<CS: ConstraintSystem<E>>(mut cs: CS) -> Result<Self, SynthesisError> {
            Ok(Self::new(
                F::zero(cs.ns(|| "zero"))?,
                F::one(cs.ns(|| "one"))?,
            ))
        }

        /// Optimized constraints for checking Edwards point addition from ZCash
        /// developers Daira Hopwood and Sean Bowe. Requires only 6 constraints
        /// compared to 7 for the straightforward version we had earlier.
        fn add<CS: ConstraintSystem<E>>(
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

        fn add_constant<CS: ConstraintSystem<E>>(
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

        fn double_in_place<CS: ConstraintSystem<E>>(
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

        fn negate<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Self, SynthesisError> {
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
            CS: ConstraintSystem<E>,
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
                    *self = <Self as GroupGadget<TEProjective<P>, E>>::add(
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

        fn cost_of_add() -> usize {
            4 + 2 * F::cost_of_mul()
        }

        fn cost_of_double() -> usize {
            4 + F::cost_of_mul()
        }
    }

    impl<P, E, F> AllocGadget<TEProjective<P>, E> for AffineGadget<P, E, F>
    where
        P: TEModelParameters,
        E: PairingEngine,
        F: FieldGadget<P::BaseField, E>,
        Self: GroupGadget<TEProjective<P>, E>,
    {
        fn alloc<FN, T, CS: ConstraintSystem<E>>(
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
                },
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

        fn alloc_checked<FN, T, CS: ConstraintSystem<E>>(
            mut cs: CS,
            value_gen: FN,
        ) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<TEProjective<P>>,
        {
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
        }

        fn alloc_input<FN, T, CS: ConstraintSystem<E>>(
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
                },
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
}

impl<P, E, F> CondSelectGadget<E> for AffineGadget<P, E, F>
where
    P: TEModelParameters,
    E: PairingEngine,
    F: FieldGadget<P::BaseField, E>,
{
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<E>>(
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
        2 * <F as CondSelectGadget<E>>::cost()
    }
}

impl<P, E, F> EqGadget<E> for AffineGadget<P, E, F>
where
    P: TEModelParameters,
    E: PairingEngine,
    F: FieldGadget<P::BaseField, E>,
{
}

impl<P, E, F> ConditionalEqGadget<E> for AffineGadget<P, E, F>
where
    P: TEModelParameters,
    E: PairingEngine,
    F: FieldGadget<P::BaseField, E>,
{
    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.x.conditional_enforce_equal(
            &mut cs.ns(|| "X Coordinate Conditional Equality"),
            &other.x,
            condition,
        )?;
        self.y.conditional_enforce_equal(
            &mut cs.ns(|| "Y Coordinate Conditional Equality"),
            &other.y,
            condition,
        )?;
        Ok(())
    }

    fn cost() -> usize {
        2 * <F as ConditionalEqGadget<E>>::cost()
    }
}

impl<P, E, F> NEqGadget<E> for AffineGadget<P, E, F>
where
    P: TEModelParameters,
    E: PairingEngine,
    F: FieldGadget<P::BaseField, E>,
{
    #[inline]
    fn enforce_not_equal<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        self.x
            .enforce_not_equal(&mut cs.ns(|| "X Coordinate Inequality"), &other.x)?;
        self.y
            .enforce_not_equal(&mut cs.ns(|| "Y Coordinate Inequality"), &other.y)?;
        Ok(())
    }

    fn cost() -> usize {
        2 * <F as NEqGadget<E>>::cost()
    }
}

impl<P, E, F> ToBitsGadget<E> for AffineGadget<P, E, F>
where
    P: TEModelParameters,
    E: PairingEngine,
    F: FieldGadget<P::BaseField, E>,
{
    fn to_bits<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        let mut x_bits = self.x.to_bits(cs.ns(|| "X Coordinate To Bits"))?;
        let y_bits = self.y.to_bits(cs.ns(|| "Y Coordinate To Bits"))?;
        x_bits.extend_from_slice(&y_bits);
        Ok(x_bits)
    }

    fn to_bits_strict<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut x_bits = self.x.to_bits_strict(cs.ns(|| "X Coordinate To Bits"))?;
        let y_bits = self.y.to_bits_strict(cs.ns(|| "Y Coordinate To Bits"))?;
        x_bits.extend_from_slice(&y_bits);

        Ok(x_bits)
    }
}

impl<P, E, F> ToBytesGadget<E> for AffineGadget<P, E, F>
where
    P: TEModelParameters,
    E: PairingEngine,
    F: FieldGadget<P::BaseField, E>,
{
    fn to_bytes<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x_bytes = self.x.to_bytes(cs.ns(|| "x"))?;
        let y_bytes = self.y.to_bytes(cs.ns(|| "y"))?;
        x_bytes.extend_from_slice(&y_bytes);
        Ok(x_bytes)
    }

    fn to_bytes_strict<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x_bytes = self.x.to_bytes_strict(cs.ns(|| "x"))?;
        let y_bytes = self.y.to_bytes_strict(cs.ns(|| "y"))?;
        x_bytes.extend_from_slice(&y_bytes);

        Ok(x_bytes)
    }
}
