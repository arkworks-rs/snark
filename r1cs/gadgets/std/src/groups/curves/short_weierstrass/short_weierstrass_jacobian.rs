use algebra::{
    curves::short_weierstrass_jacobian::{
        GroupAffine as SWAffine, GroupProjective as SWProjective,
    },
    AffineCurve, BitIterator, Field, PrimeField, ProjectiveCurve, SWModelParameters,
};
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::{borrow::Borrow, marker::PhantomData, ops::Neg};

use crate::{prelude::*, Assignment};

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

impl<P, ConstraintF, F> AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
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

    #[inline]
    /// Incomplete addition: neither `self` nor `other` can be the neutral
    /// element.
    /// If `safe` is set, enforce in the circuit exceptional cases not occurring.
    fn add_internal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        safe: bool,
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
        let x2_minus_x1 = other.x.sub(cs.ns(|| "x2 - x1"), &self.x)?;
        let y2_minus_y1 = other.y.sub(cs.ns(|| "y2 - y1"), &self.y)?;

        let lambda = if safe {
            // Check that A.x - B.x != 0, which can be done by
            // enforcing I * (B.x - A.x) = 1
            // This is done below when we calculate inv (by F::inverse)
            let inv = x2_minus_x1.inverse(cs.ns(|| "compute inv"))?;
            F::alloc(cs.ns(|| "lambda"), || {
                Ok(y2_minus_y1.get_value().get()? * &inv.get_value().get()?)
            })
        } else {
            F::alloc(cs.ns(|| "lambda"), || {
                Ok(y2_minus_y1.get_value().get()?
                    * &x2_minus_x1.get_value().get()?.inverse().get()?)
            })
        }?;

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

    #[inline]
    /// Incomplete, unsafe, addition: neither `self` nor `other` can be the neutral
    /// element.
    pub fn add_unsafe<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        self.add_internal(cs, other, false)
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
    ConstraintF: Field,
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
            Boolean::constant(true),
        ))
    }

    #[inline]
    fn is_zero<CS: ConstraintSystem<ConstraintF>>(&self, _: CS) -> Result<Boolean, SynthesisError> {
        Ok(self.infinity)
    }

    #[inline]
    /// Incomplete, safe, addition: neither `self` nor `other` can be the neutral
    /// element.
    fn add<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        self.add_internal(cs, other, true)
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

        // Allocate fresh x and y as a temporary workaround to reduce the R1CS density.
        let x = F::alloc(cs.ns(|| "new x"), || {
            let lambda_val = lambda.get_value().get()?;
            let x_val = self.x.get_value().get()?;
            Ok((lambda_val * &lambda_val) - &x_val - &x_val)
        })?;

        // lambda * lambda = new_x + 2_old_x
        let new_x_plus_two_x = self
            .x
            .add(cs.ns(|| "2old_x"), &self.x)?
            .add(cs.ns(|| "new_x + 2old_x"), &x)?;
        lambda.mul_equals(cs.ns(|| "check new x"), &lambda, &new_x_plus_two_x)?;

        let y = F::alloc(cs.ns(|| "new y"), || {
            let lambda_val = lambda.get_value().get()?;
            let x_val = self.x.get_value().get()?;
            let y_val = self.y.get_value().get()?;
            let new_x_val = x.get_value().get()?;
            Ok(((x_val - &new_x_val) * &lambda_val) - &y_val)
        })?;

        //lambda * (old_x - new_x) = new_y + old_y
        let old_x_minus_new_x = self.x.sub(cs.ns(|| "old_x - new_x"), &x)?;
        let old_y_plus_new_y = self.y.add(cs.ns(|| "old_y + new_y"), &y)?;
        lambda.mul_equals(
            cs.ns(|| "check new y"),
            &old_x_minus_new_x,
            &old_y_plus_new_y,
        )?;

        *self = Self::new(x, y, Boolean::constant(false));
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

    ///This will take [(4 + 1) * ceil(len(bits)/2)] constraints to put the x lookup constraint
    ///into the addition formula. See coda/src/lib/snarky_curves/snarky_curves.ml "scale_known"
    ///Note: `self` must be different from `result` due to SW incomplete addition.
    #[inline]
    fn mul_bits_fixed_base<'a, CS: ConstraintSystem<ConstraintF>>(
        base: &'a SWProjective<P>,
        mut cs: CS,
        result: &Self,
        bits: &[Boolean],
    ) -> Result<Self, SynthesisError> {
        let mut to_sub = SWProjective::<P>::zero();

        let mut t = base.clone();
        let sigma = base.clone();
        let mut result = result.clone();

        let mut bit_vec = Vec::new();
        bit_vec.extend_from_slice(bits);
        //Simply add padding. This should be safe, since the padding bit will be part of the
        //circuit. (It is also done elsewhere).
        if bits.len() % 2 != 0 {
            bit_vec.push(Boolean::constant(false))
        }

        for (i, bits) in bit_vec.chunks(2).enumerate() {
            let ti = t.clone();
            let two_ti = ti.double();
            let mut table = [sigma, sigma + &ti, sigma + &two_ti, sigma + &ti + &two_ti];

            //Compute constants
            SWProjective::batch_normalization(&mut table);
            let x_coords = [table[0].x, table[1].x, table[2].x, table[3].x];
            let y_coords = [table[0].y, table[1].y, table[2].y, table[3].y];
            let precomp = Boolean::and(cs.ns(|| format!("b0 AND b1_{}", i)), &bits[0], &bits[1])?;

            //Lookup x and y
            let x = F::two_bit_lookup_lc(
                cs.ns(|| format!("Lookup x_{}", i)),
                &precomp,
                &[bits[0], bits[1]],
                &x_coords,
            )?;
            let y = F::two_bit_lookup_lc(
                cs.ns(|| format!("Lookup y_{}", i)),
                &precomp,
                &[bits[0], bits[1]],
                &y_coords,
            )?;

            //Perform addition
            let adder: Self = Self::new(x, y, Boolean::constant(false));
            result = result.add(cs.ns(|| format!("Add_{}", i)), &adder)?;
            t = t.double().double();
            to_sub += &sigma;
        }
        result = result.sub_constant(cs.ns(|| "result - sigma*n_div_2"), &to_sub)?;
        Ok(result)
    }

    /// Useful in context when you have some signed representation of the scalar's digits, like
    /// in BH hash. I decided here to keep the same logic as TE implementation  for future extensibility:
    /// in fact there is no actual difference between "outer" and "inner" sums since they all
    /// are SW unsafe additions. The code could be simplified, but nothing changes from a number
    /// of constraints point of view.
    fn precomputed_base_3_bit_signed_digit_scalar_mul<'a, CS, I, J, B>(
        mut cs: CS,
        bases: &[B],
        scalars: &[J],
    ) -> Result<Self, SynthesisError>
    where
        CS: ConstraintSystem<ConstraintF>,
        I: Borrow<[Boolean]>,
        J: Borrow<[I]>,
        B: Borrow<[SWProjective<P>]>,
    {
        const CHUNK_SIZE: usize = 3;
        let mut sw_result: Option<AffineGadget<P, ConstraintF, F>> = None;
        let mut result: Option<AffineGadget<P, ConstraintF, F>> = None;
        let mut process_segment_result = |mut cs: r1cs_core::Namespace<_, _>,
                                          result: &AffineGadget<P, ConstraintF, F>|
         -> Result<(), SynthesisError> {
            let segment_result = result.clone();
            match sw_result {
                None => {
                    sw_result = Some(segment_result);
                }
                Some(ref mut sw_result) => {
                    *sw_result =
                        segment_result.add_unsafe(cs.ns(|| "sw outer addition"), sw_result)?;
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
                let coords = coords.iter().map(|p| p.into_affine()).collect::<Vec<_>>();
                let x_coeffs = coords.iter().map(|p| p.x).collect::<Vec<_>>();
                let y_coeffs = coords.iter().map(|p| p.y).collect::<Vec<_>>();
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
                let tmp = Self::new(x, y, Boolean::constant(false));
                match result {
                    None => {
                        result = Some(tmp);
                    }
                    Some(ref mut result) => {
                        *result = tmp.add_unsafe(
                            cs.ns(|| format!("addition of window {}, {}", segment_i, i)),
                            result,
                        )?;
                    }
                }
            }
            process_segment_result(cs.ns(|| format!("window {}", segment_i)), &result.unwrap())?;
            result = None;
        }
        if result.is_some() {
            process_segment_result(cs.ns(|| "leftover"), &result.unwrap())?;
        }
        Ok(sw_result.unwrap())
    }

    fn cost_of_add() -> usize {
        3 * F::cost_of_mul_equals() + F::cost_of_inv()
    }

    fn cost_of_double() -> usize {
        3 * F::cost_of_mul() + F::cost_of_mul_equals()
    }
}

impl<P, ConstraintF, F> CondSelectGadget<ConstraintF> for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
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
        let infinity = Boolean::conditionally_select(
            &mut cs.ns(|| "infinity"),
            cond,
            &first.infinity,
            &second.infinity,
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
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        let b0 = self.x.is_eq(cs.ns(|| "x"), &other.x)?;
        let b1 = self.y.is_eq(cs.ns(|| "y"), &other.y)?;
        let coordinates_equal = Boolean::and(cs.ns(|| "x AND y"), &b0, &b1)?;
        let both_are_zero = Boolean::and(
            cs.ns(|| "self.infinity AND other.infinity"),
            &self.infinity,
            &other.infinity,
        )?;
        Boolean::or(
            cs.ns(|| "coordinates_equal OR both_are_zero"),
            &coordinates_equal,
            &both_are_zero,
        )
    }

    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.is_eq(cs.ns(|| "is_eq(self, other)"), &other)?
            .conditional_enforce_equal(
                cs.ns(|| "enforce condition"),
                &Boolean::constant(true),
                &should_enforce,
            )?;
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

impl<P, ConstraintF, F> AllocGadget<SWProjective<P>, ConstraintF>
    for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
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
    fn alloc_without_check<FN, T, CS: ConstraintSystem<ConstraintF>>(
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

        let x = F::alloc(&mut cs.ns(|| "x"), || x)?;
        let y = F::alloc(&mut cs.ns(|| "y"), || y)?;
        let infinity = Boolean::alloc(&mut cs.ns(|| "infinity"), || infinity)?;

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
        let alloc_and_prime_order_check =
            |mut cs: r1cs_core::Namespace<_, _>, value_gen: FN| -> Result<Self, SynthesisError> {
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
        let ge = alloc_and_prime_order_check(cs.ns(|| "alloc and prime order check"), value_gen)?;

        Ok(ge)
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

        let b = P::COEFF_B;
        let a = P::COEFF_A;

        let x = F::alloc_input(&mut cs.ns(|| "x"), || x)?;
        let y = F::alloc_input(&mut cs.ns(|| "y"), || y)?;
        let infinity = Boolean::alloc_input(&mut cs.ns(|| "infinity"), || infinity)?;

        // Check that y^2 = x^3 + ax +b
        // We do this by checking that y^2 - b = x * (x^2 +a)
        let x2 = x.square(&mut cs.ns(|| "x^2"))?;
        let y2 = y.square(&mut cs.ns(|| "y^2"))?;

        let x2_plus_a = x2.add_constant(cs.ns(|| "x^2 + a"), &a)?;
        let y2_minus_b = y2.add_constant(cs.ns(|| "y^2 - b"), &b.neg())?;

        x2_plus_a.mul_equals(cs.ns(|| "on curve check"), &x, &y2_minus_b)?;

        Ok(Self::new(x, y, infinity))
    }
}

impl<P, ConstraintF, F> ConstantGadget<SWProjective<P>, ConstraintF>
    for AffineGadget<P, ConstraintF, F>
where
    P: SWModelParameters,
    ConstraintF: Field,
    F: FieldGadget<P::BaseField, ConstraintF>,
{
    fn from_value<CS: ConstraintSystem<ConstraintF>>(mut cs: CS, value: &SWProjective<P>) -> Self {
        let value = value.into_affine();
        let x = F::from_value(cs.ns(|| "hardcode x"), &value.x);
        let y = F::from_value(cs.ns(|| "hardcode y"), &value.y);
        let infinity = Boolean::constant(value.infinity);

        Self::new(x, y, infinity)
    }

    fn get_constant(&self) -> SWProjective<P> {
        let value_proj = SWAffine::<P>::new(
            self.x.get_value().unwrap(),
            self.y.get_value().unwrap(),
            self.infinity.get_value().unwrap(),
        )
        .into_projective();
        let x = value_proj.x;
        let y = value_proj.y;
        let z = value_proj.z;
        SWProjective::<P>::new(x, y, z)
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

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut x_bits = self
            .x
            .to_bits_strict(&mut cs.ns(|| "X Coordinate To Bits"))?;
        let y_bits = self
            .y
            .to_bits_strict(&mut cs.ns(|| "Y Coordinate To Bits"))?;
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

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x_bytes = self
            .x
            .to_bytes_strict(&mut cs.ns(|| "X Coordinate To Bytes"))?;
        let y_bytes = self
            .y
            .to_bytes_strict(&mut cs.ns(|| "Y Coordinate To Bytes"))?;
        let inf_bytes = self.infinity.to_bytes(&mut cs.ns(|| "Infinity to Bytes"))?;
        x_bytes.extend_from_slice(&y_bytes);
        x_bytes.extend_from_slice(&inf_bytes);

        Ok(x_bytes)
    }
}

#[derive(Derivative)]
#[derivative(Debug, Clone)]
#[must_use]
pub struct CompressAffinePointGadget<ConstraintF: PrimeField> {
    pub x: FpGadget<ConstraintF>,
    pub y: FpGadget<ConstraintF>,
    pub infinity: Boolean,
    _engine: PhantomData<ConstraintF>,
}

impl<ConstraintF> CompressAffinePointGadget<ConstraintF>
where
    ConstraintF: PrimeField,
{
    pub fn new(x: FpGadget<ConstraintF>, y: FpGadget<ConstraintF>, infinity: Boolean) -> Self {
        Self {
            x,
            y,
            infinity,
            _engine: PhantomData,
        }
    }
}

use crate::fields::fp::FpGadget;
use crate::ToCompressedBitsGadget;
impl<ConstraintF> ToCompressedBitsGadget<ConstraintF> for CompressAffinePointGadget<ConstraintF>
where
    ConstraintF: PrimeField,
{
    /// Enforce compression of a point through serialization of the x coordinate and storing
    /// a sign bit for the y coordinate.
    fn to_compressed<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        //Enforce x_coordinate to bytes
        let mut compressed_bits = self.x.to_bits_strict(cs.ns(|| "x_to_bits_strict"))?;
        compressed_bits.push(self.infinity);
        let is_odd = self.y.is_odd(cs.ns(|| "y parity"))?;
        compressed_bits.push(is_odd);
        Ok(compressed_bits)
    }
}
