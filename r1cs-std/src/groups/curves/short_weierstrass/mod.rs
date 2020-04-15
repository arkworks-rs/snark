use algebra::{
    curves::{
        short_weierstrass_jacobian::{GroupAffine as SWAffine, GroupProjective as SWProjective},
        SWModelParameters,
    },
    AffineCurve, BitIterator, Field, One, PrimeField, ProjectiveCurve, Zero,
};
use core::{borrow::Borrow, marker::PhantomData, ops::Neg};
use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{prelude::*, Assignment, Vec};

pub mod bls12;
pub mod mnt4;
pub mod mnt6;

#[derive(Derivative)]
#[derivative(Debug, Clone)]
#[must_use]
pub struct AffineGadget<
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
> {
    pub x: F,
    pub y: F,
    pub infinity: Boolean,
    _params: PhantomData<P>,
    _engine: PhantomData<ConstraintF>,
}

impl<P: SWModelParameters, ConstraintF: Field, F: FieldGadget<P::BaseField, ConstraintF>>
    AffineGadget<P, ConstraintF, F>
{
    pub fn new(x: F, y: F, infinity: Boolean) -> Self {
        Self {
            x,
            y,
            infinity,
            _params: PhantomData,
            _engine: PhantomData,
        }
    }

    pub fn alloc_without_check<FN, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<SWProjective<P>, SynthesisError>,
    {
        let (x, y, infinity) = match value_gen() {
            Ok(ge) => {
                let ge = ge.into_affine();
                (Ok(ge.x), Ok(ge.y), Ok(ge.infinity))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let x = F::alloc(&mut cs.ns(|| "x"), || x)?;
        let y = F::alloc(&mut cs.ns(|| "y"), || y)?;
        let infinity = Boolean::alloc(&mut cs.ns(|| "infinity"), || infinity)?;

        Ok(Self::new(x, y, infinity))
    }
}

impl<P, ConstraintF, F> PartialEq for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<P, ConstraintF, F> Eq for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
}

impl<P, ConstraintF, F> GroupGadget<SWProjective<P>, ConstraintF>
    for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: PrimeField,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    type Value = SWProjective<P>;
    type Variable = (F::Variable, F::Variable);

    #[inline]
    fn get_value(&self) -> Option<Self::Value> {
        match (
            self.x.get_value(),
            self.y.get_value(),
            self.infinity.get_value(),
        ) {
            (Some(x), Some(y), Some(infinity)) => {
                Some(SWAffine::new(x, y, infinity).into_projective())
            }
            (None, None, None) => None,
            _ => unreachable!(),
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
            Boolean::Constant(true),
        ))
    }

    #[inline]
    /// Incomplete addition: neither `self` nor `other` can be the neutral
    /// element.
    fn add<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        // lambda = (B.y - A.y)/(B.x - A.x)
        // C.x = lambda^2 - A.x - B.x
        // C.y = lambda(A.x - C.x) - A.y
        //
        // Special cases:
        //
        // doubling: if B.y = A.y and B.x = A.x then lambda is unbound and
        // C = (lambda^2, lambda^3)
        //
        // addition of negative point: if B.y = -A.y and B.x = A.x then no
        // lambda can satisfy the first equation unless B.y - A.y = 0. But
        // then this reduces to doubling.
        //
        // So we need to check that A.x - B.x != 0, which can be done by
        // enforcing I * (B.x - A.x) = 1
        // This is done below when we calculate inv (by F::inverse)

        let x2_minus_x1 = other.x.sub(cs.ns(|| "x2 - x1"), &self.x)?;
        let y2_minus_y1 = other.y.sub(cs.ns(|| "y2 - y1"), &self.y)?;

        let inv = x2_minus_x1.inverse(cs.ns(|| "compute inv"))?;

        let lambda = F::alloc(cs.ns(|| "lambda"), || {
            Ok(y2_minus_y1.get_value().get()? * &inv.get_value().get()?)
        })?;

        let x_3 = F::alloc(&mut cs.ns(|| "x_3"), || {
            let lambda_val = lambda.get_value().get()?;
            let x1 = self.x.get_value().get()?;
            let x2 = other.x.get_value().get()?;
            Ok((lambda_val.square() - &x1) - &x2)
        })?;

        let y_3 = F::alloc(&mut cs.ns(|| "y_3"), || {
            let lambda_val = lambda.get_value().get()?;
            let x_1 = self.x.get_value().get()?;
            let y_1 = self.y.get_value().get()?;
            let x_3 = x_3.get_value().get()?;
            Ok(lambda_val * &(x_1 - &x_3) - &y_1)
        })?;

        // Check lambda
        lambda.mul_equals(cs.ns(|| "check lambda"), &x2_minus_x1, &y2_minus_y1)?;

        // Check x3
        let x3_plus_x1_plus_x2 = x_3
            .add(cs.ns(|| "x3 + x1"), &self.x)?
            .add(cs.ns(|| "x3 + x1 + x2"), &other.x)?;
        lambda.mul_equals(cs.ns(|| "check x3"), &lambda, &x3_plus_x1_plus_x2)?;

        // Check y3
        let y3_plus_y1 = y_3.add(cs.ns(|| "y3 + y1"), &self.y)?;
        let x1_minus_x3 = self.x.sub(cs.ns(|| "x1 - x3"), &x_3)?;

        lambda.mul_equals(cs.ns(|| ""), &x1_minus_x3, &y3_plus_y1)?;

        Ok(Self::new(x_3, y_3, Boolean::Constant(false)))
    }

    /// Incomplete addition: neither `self` nor `other` can be the neutral
    /// element.
    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &SWProjective<P>,
    ) -> Result<Self, SynthesisError> {
        // lambda = (B.y - A.y)/(B.x - A.x)
        // C.x = lambda^2 - A.x - B.x
        // C.y = lambda(A.x - C.x) - A.y
        //
        // Special cases:
        //
        // doubling: if B.y = A.y and B.x = A.x then lambda is unbound and
        // C = (lambda^2, lambda^3)
        //
        // addition of negative point: if B.y = -A.y and B.x = A.x then no
        // lambda can satisfy the first equation unless B.y - A.y = 0. But
        // then this reduces to doubling.
        //
        // So we need to check that A.x - B.x != 0, which can be done by
        // enforcing I * (B.x - A.x) = 1
        if other.is_zero() {
            return Err(SynthesisError::AssignmentMissing);
        }
        let other = other.into_affine();
        let other_x = other.x;
        let other_y = other.y;

        let x2_minus_x1 = self
            .x
            .sub_constant(cs.ns(|| "x2 - x1"), &other_x)?
            .negate(cs.ns(|| "neg1"))?;
        let y2_minus_y1 = self
            .y
            .sub_constant(cs.ns(|| "y2 - y1"), &other_y)?
            .negate(cs.ns(|| "neg2"))?;

        let inv = x2_minus_x1.inverse(cs.ns(|| "compute inv"))?;

        let lambda = F::alloc(cs.ns(|| "lambda"), || {
            Ok(y2_minus_y1.get_value().get()? * &inv.get_value().get()?)
        })?;

        let x_3 = F::alloc(&mut cs.ns(|| "x_3"), || {
            let lambda_val = lambda.get_value().get()?;
            let x1 = self.x.get_value().get()?;
            let x2 = other_x;
            Ok((lambda_val.square() - &x1) - &x2)
        })?;

        let y_3 = F::alloc(&mut cs.ns(|| "y_3"), || {
            let lambda_val = lambda.get_value().get()?;
            let x_1 = self.x.get_value().get()?;
            let y_1 = self.y.get_value().get()?;
            let x_3 = x_3.get_value().get()?;
            Ok(lambda_val * &(x_1 - &x_3) - &y_1)
        })?;

        // Check lambda
        lambda.mul_equals(cs.ns(|| "check lambda"), &x2_minus_x1, &y2_minus_y1)?;

        // Check x3
        let x3_plus_x1_plus_x2 = x_3
            .add(cs.ns(|| "x3 + x1"), &self.x)?
            .add_constant(cs.ns(|| "x3 + x1 + x2"), &other_x)?;
        lambda.mul_equals(cs.ns(|| "check x3"), &lambda, &x3_plus_x1_plus_x2)?;

        // Check y3
        let y3_plus_y1 = y_3.add(cs.ns(|| "y3 + y1"), &self.y)?;
        let x1_minus_x3 = self.x.sub(cs.ns(|| "x1 - x3"), &x_3)?;

        lambda.mul_equals(cs.ns(|| ""), &x1_minus_x3, &y3_plus_y1)?;

        Ok(Self::new(x_3, y_3, Boolean::Constant(false)))
    }

    #[inline]
    fn double_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
    ) -> Result<(), SynthesisError> {
        let a = P::COEFF_A;
        let x_squared = self.x.square(cs.ns(|| "x^2"))?;

        let one = P::BaseField::one();
        let two = one.double();
        let three = two + &one;

        let three_x_squared = x_squared.mul_by_constant(cs.ns(|| "3 * x^2"), &three)?;
        let three_x_squared_plus_a = three_x_squared.add_constant(cs.ns(|| "3 * x^2 + a"), &a)?;

        let two_y = self.y.double(cs.ns(|| "2y"))?;

        let lambda = F::alloc(cs.ns(|| "lambda"), || {
            let y_doubled_inv = two_y.get_value().get()?.inverse().get()?;
            Ok(three_x_squared_plus_a.get_value().get()? * &y_doubled_inv)
        })?;

        // Check lambda
        lambda.mul_equals(cs.ns(|| "check lambda"), &two_y, &three_x_squared_plus_a)?;

        let x = lambda
            .square(cs.ns(|| "lambda^2"))?
            .sub(cs.ns(|| "lambda^2 - x"), &self.x)?
            .sub(cs.ns(|| "lambda^2 - 2x"), &self.x)?;

        let y = self
            .x
            .sub(cs.ns(|| "x - self.x"), &x)?
            .mul(cs.ns(|| "times lambda"), &lambda)?
            .sub(cs.ns(|| "plus self.y"), &self.y)?;

        *self = Self::new(x, y, Boolean::Constant(false));
        Ok(())
    }

    fn negate<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
        Ok(Self::new(
            self.x.clone(),
            self.y.negate(cs.ns(|| "negate y"))?,
            self.infinity,
        ))
    }

    fn cost_of_add() -> usize {
        3 * F::cost_of_mul_equals() + F::cost_of_inv()
    }

    fn cost_of_double() -> usize {
        4 * F::cost_of_mul()
    }
}

impl<P, ConstraintF, F> CondSelectGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: PrimeField,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        true_value: &Self,
        false_value: &Self,
    ) -> Result<Self, SynthesisError> {
        let x = F::conditionally_select(&mut cs.ns(|| "x"), cond, &true_value.x, &false_value.x)?;
        let y = F::conditionally_select(&mut cs.ns(|| "y"), cond, &true_value.y, &false_value.y)?;
        let infinity = Boolean::conditionally_select(
            &mut cs.ns(|| "infinity"),
            cond,
            &true_value.infinity,
            &false_value.infinity,
        )?;

        Ok(Self::new(x, y, infinity))
    }

    fn cost() -> usize {
        2 * <F as CondSelectGadget<ConstraintF>>::cost()
            + <Boolean as CondSelectGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF, F> EqGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
}

impl<P, ConstraintF, F> ConditionalEqGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
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
        self.infinity.conditional_enforce_equal(
            &mut cs.ns(|| "Infinity Conditional Equality"),
            &other.infinity,
            condition,
        )?;
        Ok(())
    }

    fn cost() -> usize {
        2 * <F as ConditionalEqGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF, F> NEqGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    #[inline]
    fn enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
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
        2 * <F as NEqGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF, F> AllocGadget<SWProjective<P>, ConstraintF>
    for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: PrimeField,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    #[inline]
    fn alloc_constant<T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<SWProjective<P>>,
    {
        let p = t.borrow().into_affine();
        Ok(Self {
            x: F::alloc_constant(cs.ns(|| "x"), &p.x)?,
            y: F::alloc_constant(cs.ns(|| "y"), &p.y)?,
            infinity: Boolean::constant(p.infinity),
            _params: PhantomData,
            _engine: PhantomData,
        })
    }

    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SWProjective<P>>,
    {
        let (x, y, infinity) = match value_gen() {
            Ok(ge) => {
                let ge = ge.borrow().into_affine();
                (Ok(ge.x), Ok(ge.y), Ok(ge.infinity))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        // Perform on-curve check.
        let b = P::COEFF_B;
        let a = P::COEFF_A;

        let x = F::alloc(&mut cs.ns(|| "x"), || x)?;
        let y = F::alloc(&mut cs.ns(|| "y"), || y)?;
        let infinity = Boolean::alloc(&mut cs.ns(|| "infinity"), || infinity)?;

        // Check that y^2 = x^3 + ax +b
        // We do this by checking that y^2 - b = x * (x^2 +a)
        let x2 = x.square(&mut cs.ns(|| "x^2"))?;
        let y2 = y.square(&mut cs.ns(|| "y^2"))?;

        let x2_plus_a = x2.add_constant(cs.ns(|| "x^2 + a"), &a)?;
        let y2_minus_b = y2.add_constant(cs.ns(|| "y^2 - b"), &b.neg())?;

        x2_plus_a.mul_equals(cs.ns(|| "on curve check"), &x, &y2_minus_b)?;

        Ok(Self::new(x, y, infinity))
    }

    #[inline]
    fn alloc_checked<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SWProjective<P>>,
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

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SWProjective<P>>,
    {
        // When allocating the input we assume that the verifier has performed
        // any on curve checks already.
        let (x, y, infinity) = match value_gen() {
            Ok(ge) => {
                let ge = ge.borrow().into_affine();
                (Ok(ge.x), Ok(ge.y), Ok(ge.infinity))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let x = F::alloc_input(&mut cs.ns(|| "x"), || x)?;
        let y = F::alloc_input(&mut cs.ns(|| "y"), || y)?;
        let infinity = Boolean::alloc_input(&mut cs.ns(|| "infinity"), || infinity)?;

        Ok(Self::new(x, y, infinity))
    }
}

impl<P, ConstraintF, F> ToBitsGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut x_bits = self.x.to_bits(&mut cs.ns(|| "X Coordinate To Bits"))?;
        let y_bits = self.y.to_bits(&mut cs.ns(|| "Y Coordinate To Bits"))?;
        x_bits.extend_from_slice(&y_bits);
        x_bits.push(self.infinity);
        Ok(x_bits)
    }

    fn to_non_unique_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut x_bits = self
            .x
            .to_non_unique_bits(&mut cs.ns(|| "X Coordinate To Bits"))?;
        let y_bits = self
            .y
            .to_non_unique_bits(&mut cs.ns(|| "Y Coordinate To Bits"))?;
        x_bits.extend_from_slice(&y_bits);
        x_bits.push(self.infinity);

        Ok(x_bits)
    }
}

impl<P, ConstraintF, F> ToBytesGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x_bytes = self.x.to_bytes(&mut cs.ns(|| "X Coordinate To Bytes"))?;
        let y_bytes = self.y.to_bytes(&mut cs.ns(|| "Y Coordinate To Bytes"))?;
        let inf_bytes = self.infinity.to_bytes(&mut cs.ns(|| "Infinity to Bytes"))?;
        x_bytes.extend_from_slice(&y_bytes);
        x_bytes.extend_from_slice(&inf_bytes);
        Ok(x_bytes)
    }

    fn to_non_unique_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x_bytes = self
            .x
            .to_non_unique_bytes(&mut cs.ns(|| "X Coordinate To Bytes"))?;
        let y_bytes = self
            .y
            .to_non_unique_bytes(&mut cs.ns(|| "Y Coordinate To Bytes"))?;
        let inf_bytes = self
            .infinity
            .to_non_unique_bytes(&mut cs.ns(|| "Infinity to Bytes"))?;
        x_bytes.extend_from_slice(&y_bytes);
        x_bytes.extend_from_slice(&inf_bytes);

        Ok(x_bytes)
    }
}

#[cfg(test)]
#[allow(dead_code)]
pub(crate) fn test<ConstraintF, P, GG>()
where
    ConstraintF: Field,
    P: SWModelParameters,
    GG: GroupGadget<SWProjective<P>, ConstraintF, Value = SWProjective<P>>,
{
    use crate::{boolean::AllocatedBit, prelude::*, test_constraint_system::TestConstraintSystem};
    use algebra::{test_rng, Group, UniformRand};
    use rand::Rng;

    // Incomplete addition doesn't allow us to call the group_test.
    // group_test::<ConstraintF, SWProjective<P>, GG>();

    let mut rng = test_rng();

    let mut cs = TestConstraintSystem::<ConstraintF>::new();

    let a = SWProjective::<P>::rand(&mut rng);
    let b = SWProjective::<P>::rand(&mut rng);
    let a_affine = a.into_affine();
    let b_affine = b.into_affine();
    let mut gadget_a = GG::alloc(&mut cs.ns(|| "a"), || Ok(a)).unwrap();
    let gadget_b = GG::alloc_checked(&mut cs.ns(|| "b"), || Ok(b)).unwrap();
    assert_eq!(gadget_a.get_value().unwrap().x, a_affine.x);
    assert_eq!(gadget_a.get_value().unwrap().y, a_affine.y);
    assert_eq!(gadget_b.get_value().unwrap().x, b_affine.x);
    assert_eq!(gadget_b.get_value().unwrap().y, b_affine.y);

    // Check addition
    let ab = a + &b;
    let ab_affine = ab.into_affine();
    let gadget_ab = gadget_a.add(&mut cs.ns(|| "ab"), &gadget_b).unwrap();
    let gadget_ba = gadget_b.add(&mut cs.ns(|| "ba"), &gadget_a).unwrap();
    gadget_ba
        .enforce_equal(&mut cs.ns(|| "b + a == a + b?"), &gadget_ab)
        .unwrap();

    let ab_val = gadget_ab
        .get_value()
        .expect("Doubling should be successful")
        .into_affine();
    assert_eq!(ab_val, ab_affine, "Result of addition is unequal");

    // Check doubling
    let aa = Group::double(&a);
    let aa_affine = aa.into_affine();
    gadget_a.double_in_place(&mut cs.ns(|| "2a")).unwrap();
    let aa_val = gadget_a
        .get_value()
        .expect("Doubling should be successful")
        .into_affine();
    assert_eq!(
        aa_val, aa_affine,
        "Gadget and native values are unequal after double."
    );

    // Check mul_bits
    let scalar = P::ScalarField::rand(&mut rng);
    let native_result = aa.into_affine().mul(scalar) + &b;
    let native_result = native_result.into_affine();

    let mut scalar: Vec<bool> = BitIterator::new(scalar.into_repr()).collect();
    // Get the scalar bits into little-endian form.
    scalar.reverse();
    let input = Vec::<Boolean>::alloc(cs.ns(|| "Input"), || Ok(scalar)).unwrap();
    let result = gadget_a
        .mul_bits(cs.ns(|| "mul_bits"), &gadget_b, input.iter())
        .unwrap();
    let result_val = result.get_value().unwrap().into_affine();
    assert_eq!(
        result_val, native_result,
        "gadget & native values are diff. after scalar mul"
    );

    if !cs.is_satisfied() {
        println!("{:?}", cs.which_is_unsatisfied().unwrap());
    }

    assert!(cs.is_satisfied());

    // Constraint cost etc.
    let mut cs = TestConstraintSystem::<ConstraintF>::new();

    let bit = AllocatedBit::alloc(&mut cs.ns(|| "bool"), || Ok(true))
        .unwrap()
        .into();

    let mut rng = test_rng();
    let a: SWProjective<P> = rng.gen();
    let b: SWProjective<P> = rng.gen();
    let gadget_a = GG::alloc(&mut cs.ns(|| "a"), || Ok(a)).unwrap();
    let gadget_b = GG::alloc(&mut cs.ns(|| "b"), || Ok(b)).unwrap();
    let alloc_cost = cs.num_constraints();
    let _ =
        GG::conditionally_select(&mut cs.ns(|| "cond_select"), &bit, &gadget_a, &gadget_b).unwrap();
    let cond_select_cost = cs.num_constraints() - alloc_cost;

    let _ = gadget_a.add(&mut cs.ns(|| "ab"), &gadget_b).unwrap();
    let add_cost = cs.num_constraints() - cond_select_cost - alloc_cost;

    assert!(cs.is_satisfied());
    assert_eq!(
        cond_select_cost,
        <GG as CondSelectGadget<ConstraintF>>::cost()
    );
    assert_eq!(add_cost, GG::cost_of_add());
}
