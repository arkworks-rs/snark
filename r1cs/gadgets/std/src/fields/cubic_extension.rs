use algebra::{CubicExtField, CubicExtParameters, Field, PrimeField, SquareRootField};
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::{borrow::Borrow, marker::PhantomData};

use crate::{fields::FieldGadget, prelude::*, Assignment};

pub trait CubicExtParametersGadget<ConstraintF: PrimeField>:
    CubicExtParameters<BasePrimeField = ConstraintF>
{
    type BaseFieldGadget: FieldGadget<Self::BaseField, ConstraintF>;

    /// Multiply a BaseFieldGadget by quadratic nonresidue.
    fn mul_base_field_gadget_by_nonresidue<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        fe: &Self::BaseFieldGadget,
    ) -> Result<Self::BaseFieldGadget, SynthesisError>;

    /// Multiply a BaseFieldGadget by the Frobenius Coefficient at given power
    fn mul_base_field_gadget_by_frobenius_coeff<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        c1: &mut Self::BaseFieldGadget,
        c2: &mut Self::BaseFieldGadget,
        power: usize,
    ) -> Result<(), SynthesisError>;
}

#[derive(Derivative)]
#[derivative(Debug(
    bound = "P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField"
))]
#[must_use]
pub struct CubicExtFieldGadget<
    P: CubicExtParametersGadget<ConstraintF>,
    ConstraintF: PrimeField + SquareRootField,
> {
    pub c0: P::BaseFieldGadget,
    pub c1: P::BaseFieldGadget,
    pub c2: P::BaseFieldGadget,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    CubicExtFieldGadget<P, ConstraintF>
{
    #[inline]
    pub fn new(c0: P::BaseFieldGadget, c1: P::BaseFieldGadget, c2: P::BaseFieldGadget) -> Self {
        Self {
            c0,
            c1,
            c2,
            _params: PhantomData,
        }
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    FieldGadget<CubicExtField<P>, ConstraintF> for CubicExtFieldGadget<P, ConstraintF>
{
    type Variable = (
        <P::BaseFieldGadget as FieldGadget<P::BaseField, ConstraintF>>::Variable,
        <P::BaseFieldGadget as FieldGadget<P::BaseField, ConstraintF>>::Variable,
        <P::BaseFieldGadget as FieldGadget<P::BaseField, ConstraintF>>::Variable,
    );

    #[inline]
    fn get_value(&self) -> Option<CubicExtField<P>> {
        match (
            self.c0.get_value(),
            self.c1.get_value(),
            self.c2.get_value(),
        ) {
            (Some(c0), Some(c1), Some(c2)) => Some(CubicExtField::<P>::new(c0, c1, c2)),
            (..) => None,
        }
    }

    #[inline]
    fn get_variable(&self) -> Self::Variable {
        (
            self.c0.get_variable(),
            self.c1.get_variable(),
            self.c2.get_variable(),
        )
    }

    #[inline]
    fn zero<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = P::BaseFieldGadget::zero(cs.ns(|| "c0"))?;
        let c1 = P::BaseFieldGadget::zero(cs.ns(|| "c1"))?;
        let c2 = P::BaseFieldGadget::zero(cs.ns(|| "c2"))?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn one<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = P::BaseFieldGadget::one(cs.ns(|| "c0"))?;
        let c1 = P::BaseFieldGadget::zero(cs.ns(|| "c1"))?;
        let c2 = P::BaseFieldGadget::zero(cs.ns(|| "c2"))?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn conditionally_add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        bit: &Boolean,
        coeff: CubicExtField<P>,
    ) -> Result<Self, SynthesisError> {
        let c0 = self
            .c0
            .conditionally_add_constant(cs.ns(|| "c0"), bit, coeff.c0)?;
        let c1 = self
            .c1
            .conditionally_add_constant(cs.ns(|| "c1"), bit, coeff.c1)?;
        let c2 = self
            .c2
            .conditionally_add_constant(cs.ns(|| "c2"), bit, coeff.c2)?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn add<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = self.c0.add(&mut cs.ns(|| "add c0"), &other.c0)?;
        let c1 = self.c1.add(&mut cs.ns(|| "add c1"), &other.c1)?;
        let c2 = self.c2.add(&mut cs.ns(|| "add c2"), &other.c2)?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn sub<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = self.c0.sub(&mut cs.ns(|| "sub c0"), &other.c0)?;
        let c1 = self.c1.sub(&mut cs.ns(|| "sub c1"), &other.c1)?;
        let c2 = self.c2.sub(&mut cs.ns(|| "sub c2"), &other.c2)?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn negate<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
        let c0 = self.c0.negate(&mut cs.ns(|| "negate c0"))?;
        let c1 = self.c1.negate(&mut cs.ns(|| "negate c1"))?;
        let c2 = self.c2.negate(&mut cs.ns(|| "negate c2"))?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn negate_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.negate_in_place(&mut cs.ns(|| "negate c0"))?;
        self.c1.negate_in_place(&mut cs.ns(|| "negate c1"))?;
        self.c2.negate_in_place(&mut cs.ns(|| "negate c2"))?;
        Ok(self)
    }

    /// Use the Toom-Cook-3x method to compute multiplication.
    #[inline]
    fn mul<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        // Uses Toom-Cool-3x multiplication
        //
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        //    Devegili, OhEigeartaigh, Scott, Dahab

        // v0 = a(0)b(0)   = a0 * b0
        let v0 = self.c0.mul(&mut cs.ns(|| "Calc v0"), &other.c0)?;

        // v1 = a(1)b(1)   = (a0 + a1 + a2)(b0 + b1 + b2)
        let v1 = {
            let mut v1_cs = cs.ns(|| "compute v1");
            let a0_plus_a1_plus_a2 = self
                .c0
                .add(v1_cs.ns(|| "a0 + a1"), &self.c1)?
                .add(v1_cs.ns(|| "a0 + a1 + a2"), &self.c2)?;
            let b0_plus_b1_plus_b2 = other
                .c0
                .add(v1_cs.ns(|| "b0 + b1"), &other.c1)?
                .add(v1_cs.ns(|| "b0 + b1 + b2"), &other.c2)?;

            a0_plus_a1_plus_a2.mul(
                v1_cs.ns(|| "(a0 + a1 + a2)(b0 + b1 + b2)"),
                &b0_plus_b1_plus_b2,
            )?
        };

        // v2 = a(−1)b(−1) = (a0 − a1 + a2)(b0 − b1 + b2)
        let v2 = {
            let mut v2_cs = cs.ns(|| "compute v2");

            let a0_minus_a1_plus_a2 = self
                .c0
                .sub(v2_cs.ns(|| "a0 - a1"), &self.c1)?
                .add(v2_cs.ns(|| "a0 - a1 + a2"), &self.c2)?;

            let b0_minus_b1_plus_b2 = other
                .c0
                .sub(v2_cs.ns(|| "b0 - b1"), &other.c1)?
                .add(v2_cs.ns(|| "b0 - b1 + b2"), &other.c2)?;

            a0_minus_a1_plus_a2.mul(
                v2_cs.ns(|| "(a0 - a1 + a2)(b0 - b1 + b2)"),
                &b0_minus_b1_plus_b2,
            )?
        };

        // v3 = a(2)b(2)   = (a0 + 2a1 + 4a2)(b0 + 2b1 + 4b2)
        let v3 = {
            let v3_cs = &mut cs.ns(|| "compute v3");

            let a1_double = self.c1.double(v3_cs.ns(|| "2 * a1"))?;
            let a2_quad = self
                .c2
                .double(v3_cs.ns(|| "2 * a2"))?
                .double(v3_cs.ns(|| "4 * a2"))?;

            let a0_plus_2_a1_plus_4_a2 = self
                .c0
                .add(v3_cs.ns(|| "a0 + 2a1"), &a1_double)?
                .add(v3_cs.ns(|| "a0 + 2a1 + 4a2"), &a2_quad)?;

            let b1_double = other.c1.double(v3_cs.ns(|| "2 * b1"))?;
            let b2_quad = other
                .c2
                .double(v3_cs.ns(|| "2 * b2"))?
                .double(v3_cs.ns(|| "4 * b2"))?;
            let b0_plus_2_b1_plus_4_b2 = other
                .c0
                .add(v3_cs.ns(|| "b0 + 2b1"), &b1_double)?
                .add(v3_cs.ns(|| "b0 + 2b1 + 4b2"), &b2_quad)?;

            a0_plus_2_a1_plus_4_a2.mul(
                v3_cs.ns(|| "(a0 + 2a1 + 4a2)(b0 + 2b1 + 4b2)"),
                &b0_plus_2_b1_plus_4_b2,
            )?
        };

        // v4 = a(∞)b(∞)   = a2 * b2
        let v4 = self.c2.mul(cs.ns(|| "v2: a2 * b2"), &other.c2)?;

        let two = P::BaseField::one().double();
        let six = two.double() + &two;
        let mut two_and_six = [two, six];
        algebra::fields::batch_inversion(&mut two_and_six);
        let (two_inverse, six_inverse) = (two_and_six[0], two_and_six[1]);

        let half_v0 = v0.mul_by_constant(cs.ns(|| "half_v0"), &two_inverse)?;
        let half_v1 = v1.mul_by_constant(cs.ns(|| "half_v1"), &two_inverse)?;
        let one_sixth_v2 = v2.mul_by_constant(cs.ns(|| "v2_by_six"), &six_inverse)?;
        let one_sixth_v3 = v3.mul_by_constant(cs.ns(|| "v3_by_six"), &six_inverse)?;
        let two_v4 = v4.double(cs.ns(|| "2 * v4"))?;

        // c0 = v0 + β((1/2)v0 − (1/2)v1 − (1/6)v2 + (1/6)v3 − 2v4)
        let c0 = {
            let c0_cs = &mut cs.ns(|| "c0");

            // No constraints, only get a linear combination back.
            let temp = half_v0
                .sub(c0_cs.ns(|| "sub1"), &half_v1)?
                .sub(c0_cs.ns(|| "sub2"), &one_sixth_v2)?
                .add(c0_cs.ns(|| "add3"), &one_sixth_v3)?
                .sub(c0_cs.ns(|| "sub4"), &two_v4)?;
            let non_residue_times_inner =
                temp.mul_by_constant(&mut c0_cs.ns(|| "mul5"), &P::NONRESIDUE)?;
            v0.add(c0_cs.ns(|| "add6"), &non_residue_times_inner)?
        };

        // −(1/2)v0 + v1 − (1/3)v2 − (1/6)v3 + 2v4 + βv4
        let c1 = {
            let c1_cs = &mut cs.ns(|| "c1");
            let one_third_v2 = one_sixth_v2.double(&mut c1_cs.ns(|| "v2_by_3"))?;
            let non_residue_v4 =
                v4.mul_by_constant(&mut c1_cs.ns(|| "mul_by_beta"), &P::NONRESIDUE)?;

            let result = half_v0
                .negate(c1_cs.ns(|| "neg1"))?
                .add(c1_cs.ns(|| "add2"), &v1)?
                .sub(c1_cs.ns(|| "sub3"), &one_third_v2)?
                .sub(c1_cs.ns(|| "sub4"), &one_sixth_v3)?
                .add(c1_cs.ns(|| "sub5"), &two_v4)?
                .add(c1_cs.ns(|| "sub6"), &non_residue_v4)?;
            result
        };

        // -v0 + (1/2)v1 + (1/2)v2 −v4
        let c2 = {
            let c2_cs = &mut cs.ns(|| "c2");
            let half_v2 = v2.mul_by_constant(&mut c2_cs.ns(|| "mul1"), &two_inverse)?;
            let result = half_v1
                .add(c2_cs.ns(|| "add1"), &half_v2)?
                .sub(c2_cs.ns(|| "sub1"), &v4)?
                .sub(c2_cs.ns(|| "sub2"), &v0)?;
            result
        };
        Ok(Self::new(c0, c1, c2))
    }

    /// Use the Chung-Hasan asymmetric squaring formula.
    ///
    /// (Devegili OhEig Scott Dahab --- Multiplication and Squaring on
    /// Abstract Pairing-Friendly
    /// Fields.pdf; Section 4 (CH-SQR2))
    #[inline]
    fn square<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
        let a = self.c0.clone();
        let b = self.c1.clone();
        let c = self.c2.clone();

        let s0 = a.square(cs.ns(|| "s0"))?;
        let ab = a.mul(cs.ns(|| "ab"), &b)?;
        let s1 = ab.double(cs.ns(|| "s1"))?;
        let s2 = a
            .sub(cs.ns(|| "a-b"), &b)?
            .add(cs.ns(|| "plus c"), &c)?
            .square(cs.ns(|| "s2"))?;
        let s3 = b.mul(cs.ns(|| "bc"), &c)?.double(cs.ns(|| "s3"))?;
        let s4 = c.square(cs.ns(|| "s4"))?;

        let c0 = P::mul_base_field_gadget_by_nonresidue(cs.ns(|| "c0 part 1"), &s3)?
            .add(cs.ns(|| "c0"), &s0)?;

        let c1 = P::mul_base_field_gadget_by_nonresidue(cs.ns(|| "c1 part 1"), &s4)?
            .add(cs.ns(|| "c1"), &s1)?;

        let c2 = s1
            .add(cs.ns(|| "c2 part1"), &s2)?
            .add(cs.ns(|| "c2 part2"), &s3)?
            .sub(cs.ns(|| "c2 part3"), &s0)?
            .sub(cs.ns(|| "c2 part4"), &s4)?;

        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn mul_equals<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        // Uses Toom-Cook-3x multiplication
        //
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        //    Devegili, OhEigeartaigh, Scott, Dahab

        // The prover chooses lambda_1 and lambda_2.
        // We can explicitly enforce lambda_1 without wasting constraints as it turns out
        // that a(∞)b(∞) = lambda_1 at X = ∞.
        let lambda_1 = self.c2.mul(cs.ns(|| "lambda_1 <=> check 5"), &other.c2)?;

        let lambda_2 = P::BaseFieldGadget::alloc(cs.ns(|| "lambda_2"), || {
            let a1b2 = self.c1.get_value().get()? * &other.c2.get_value().get()?;
            let a2b1 = self.c2.get_value().get()? * &other.c1.get_value().get()?;
            Ok(a1b2 + &a2b1)
        })?;

        let one = P::BaseField::one();

        // a0b0 = c0 - β*lambda_2 at X = 0
        {
            let c0_plus_nr_lambda_2 = lambda_2
                .mul_by_constant(cs.ns(|| "nr * lambda_2"), &P::NONRESIDUE)?
                .negate(cs.ns(|| "-(nr * lambda_2)"))?
                .add(cs.ns(|| "c0 - nr * lambda_2"), &result.c0)?;

            self.c0
                .mul_equals(cs.ns(|| "check 1"), &other.c0, &c0_plus_nr_lambda_2)?;
        }

        //(a0 + a1 + a2)(b0 + b1 + b2) = (c0 + c1 + c2) + (lambda_1 + lambda_2)*(1 - β) at X = 1
        {
            let a0_plus_a1_plus_a2 = self
                .c0
                .add(cs.ns(|| "a0 + a1"), &self.c1)?
                .add(cs.ns(|| "a0 + a1 + a2"), &self.c2)?;
            let b0_plus_b1_plus_b2 = other
                .c0
                .add(cs.ns(|| "b0 + b1"), &other.c1)?
                .add(cs.ns(|| "b0 + b1 + b2"), &other.c2)?;
            let c0_plus_c1_plus_c2 = result
                .c0
                .add(cs.ns(|| "c0 + c1"), &result.c1)?
                .add(cs.ns(|| "c0 + c1 + c2"), &result.c2)?;
            let lambda_1_plus_lambda_2_times_one_minus_nr = lambda_1
                .add(cs.ns(|| "lambda_1 + lambda_2"), &lambda_2)?
                .mul_by_constant(
                    cs.ns(|| "(lambda_1 + lambda_2)*(1 - nr)"),
                    &(one - &P::NONRESIDUE),
                )?;

            let to_check = c0_plus_c1_plus_c2.add(
                cs.ns(|| "c0 + c1 + c2 + (lambda_1 + lambda_2)*(1 - nr)"),
                &lambda_1_plus_lambda_2_times_one_minus_nr,
            )?;
            a0_plus_a1_plus_a2.mul_equals(cs.ns(|| "check 2"), &b0_plus_b1_plus_b2, &to_check)?;
        }

        //(a0 - a1 + a2)(b0 - b1 + b2) = (c0 - c1 + c2) + (lambda_1 - lambda_2)*(1 + β) at X = -1
        {
            let a0_minus_a1_plus_a2 = self
                .c0
                .sub(cs.ns(|| "a0 - a1"), &self.c1)?
                .add(cs.ns(|| "a0 - a1 + a2"), &self.c2)?;
            let b0_minus_b1_plus_b2 = other
                .c0
                .sub(cs.ns(|| "b0 - b1"), &other.c1)?
                .add(cs.ns(|| "b0 - b1 + b2"), &other.c2)?;
            let c0_minus_c1_plus_c2 = result
                .c0
                .sub(cs.ns(|| "c0 - c1"), &result.c1)?
                .add(cs.ns(|| "c0 - c1 + c2"), &result.c2)?;
            let lambda_1_minus_lambda_2_times_one_plus_nr = lambda_1
                .sub(cs.ns(|| "lambda_1 - lambda_2"), &lambda_2)?
                .mul_by_constant(
                    cs.ns(|| "(lambda_1 - lambda_2)*(1 + nr)"),
                    &(one + &P::NONRESIDUE),
                )?;

            let to_check = c0_minus_c1_plus_c2.add(
                cs.ns(|| "c0 - c1 + c2 + (lambda_1 - lambda_2)*(1 + nr)"),
                &lambda_1_minus_lambda_2_times_one_plus_nr,
            )?;
            a0_minus_a1_plus_a2.mul_equals(cs.ns(|| "check 3"), &b0_minus_b1_plus_b2, &to_check)?;
        }

        // (a0 + 2a1 + 4a2)(b0 + 2b1 + 4b2) = (c0 + 2c1 + 4c2) + (2lambda_1 + lambda_2)(8 - β) at X = 2
        {
            let a0_plus_2_a1_plus_4_a2 = {
                let a1_double = self.c1.double(cs.ns(|| "2 * a1"))?;
                let a2_quad = self
                    .c2
                    .double(cs.ns(|| "2 * a2"))?
                    .double(cs.ns(|| "4 * a2"))?;

                self.c0
                    .add(cs.ns(|| "a0 + 2a1"), &a1_double)?
                    .add(cs.ns(|| "a0 + 2a1 + 4a2"), &a2_quad)?
            };

            let b0_plus_2_b1_plus_4_b2 = {
                let b1_double = other.c1.double(cs.ns(|| "2 * b1"))?;
                let b2_quad = other
                    .c2
                    .double(cs.ns(|| "2 * b2"))?
                    .double(cs.ns(|| "4 * b2"))?;

                other
                    .c0
                    .add(cs.ns(|| "b0 + 2b1"), &b1_double)?
                    .add(cs.ns(|| "b0 + 2b1 + 4b2"), &b2_quad)?
            };

            let c0_plus_2_c1_plus_4_c2 = {
                let c1_double = result.c1.double(cs.ns(|| "2 * c1"))?;
                let c2_quad = result
                    .c2
                    .double(cs.ns(|| "2 * c2"))?
                    .double(cs.ns(|| "4 * c2"))?;

                result
                    .c0
                    .add(cs.ns(|| "c0 + 2c1"), &c1_double)?
                    .add(cs.ns(|| "c0 + 2c1 + 4c2"), &c2_quad)?
            };

            let eight = one.double().double().double();
            let two_lambda_1_plus_lambda_2_times_eight_minus_nr = lambda_1
                .double(cs.ns(|| "2*lambda_1"))?
                .add(cs.ns(|| "2*lambda_1 + lambda_2"), &lambda_2)?
                .mul_by_constant(
                    cs.ns(|| "(2*lambda_1 + lambda_2)*(8 - nr)"),
                    &(eight - &P::NONRESIDUE),
                )?;

            let to_check = c0_plus_2_c1_plus_4_c2.add(
                cs.ns(|| "(c0 + 2c1 + 4c2) + (2*lambda_1 + lambda_2)*(8 - nr)"),
                &two_lambda_1_plus_lambda_2_times_eight_minus_nr,
            )?;
            a0_plus_2_a1_plus_4_a2.mul_equals(
                cs.ns(|| "check 4"),
                &b0_plus_2_b1_plus_4_b2,
                &to_check,
            )?;
        }

        Ok(())
    }

    #[inline]
    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &CubicExtField<P>,
    ) -> Result<Self, SynthesisError> {
        let c0 = self.c0.add_constant(cs.ns(|| "c0"), &other.c0)?;
        let c1 = self.c1.add_constant(cs.ns(|| "c1"), &other.c1)?;
        let c2 = self.c2.add_constant(cs.ns(|| "c2"), &other.c2)?;

        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn add_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        other: &CubicExtField<P>,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.add_constant_in_place(cs.ns(|| "c0"), &other.c0)?;
        self.c1.add_constant_in_place(cs.ns(|| "c1"), &other.c1)?;
        self.c2.add_constant_in_place(cs.ns(|| "c2"), &other.c2)?;
        Ok(self)
    }

    #[inline]
    fn mul_by_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &CubicExtField<P>,
    ) -> Result<Self, SynthesisError> {
        // Naive Fp3 multiplication

        // c0 = b0*a0 + β*b2*a1 + β*b1*a2,
        let c0 = {
            let a0_b0 = self.c0.mul_by_constant(cs.ns(|| "a0 * b0"), &other.c0)?;
            let a1_b2_nr = self
                .c1
                .mul_by_constant(cs.ns(|| "a1 * b2 * nr"), &(other.c2 * &P::NONRESIDUE))?;
            let a2_b1_nr = self
                .c2
                .mul_by_constant(cs.ns(|| "a2 * b1 * nr"), &(other.c1 * &P::NONRESIDUE))?;

            a0_b0
                .add(cs.ns(|| "a0 * b0 + a1 * b2 * nr"), &a1_b2_nr)?
                .add(cs.ns(|| "a0 * b0 + a1 * b2 * nr + a2 * b1 * nr"), &a2_b1_nr)
        }?;

        // c1 = b1*a0 + b0*a1 + β*b2*a2,
        let c1 = {
            let a0_b1 = self.c0.mul_by_constant(cs.ns(|| "a0 * b1"), &other.c1)?;
            let a1_b0 = self.c1.mul_by_constant(cs.ns(|| "a1 * b0"), &other.c0)?;
            let a2_b2_nr = self
                .c2
                .mul_by_constant(cs.ns(|| "a2 * b2 * nr"), &(other.c2 * &P::NONRESIDUE))?;

            a0_b1
                .add(cs.ns(|| "a0 * b1 + a1 * b0"), &a1_b0)?
                .add(cs.ns(|| "a0 * b1 + a1 * b0 + a2 * b2 * nr"), &a2_b2_nr)
        }?;

        // c2 = b2*a0 + b1*a1 + b0*a2.
        let c2 = {
            let a0_b2 = self.c0.mul_by_constant(cs.ns(|| "a0 * b2"), &other.c2)?;
            let a1_b1 = self.c1.mul_by_constant(cs.ns(|| "a1 * b1"), &other.c1)?;
            let a2_b0 = self.c2.mul_by_constant(cs.ns(|| "a2 * b0"), &other.c0)?;

            a0_b2
                .add(cs.ns(|| "a0 * b2 + a1 * b1"), &a1_b1)?
                .add(cs.ns(|| "a0 * b2 + a1 * b1 + a2 * b0"), &a2_b0)
        }?;

        Ok(Self::new(c0, c1, c2))
    }

    fn frobenius_map<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        power: usize,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.frobenius_map_in_place(cs, power)?;
        Ok(result)
    }

    fn frobenius_map_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        power: usize,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.frobenius_map_in_place(&mut cs.ns(|| "c0"), power)?;
        self.c1.frobenius_map_in_place(&mut cs.ns(|| "c1"), power)?;
        self.c2.frobenius_map_in_place(&mut cs.ns(|| "c2"), power)?;

        P::mul_base_field_gadget_by_frobenius_coeff(
            &mut cs.ns(|| "c1 and c2 powers"),
            &mut self.c1,
            &mut self.c2,
            power,
        )?;

        Ok(self)
    }

    fn cost_of_mul() -> usize {
        5
    }

    fn cost_of_mul_equals() -> usize {
        5
    }

    fn cost_of_inv() -> usize {
        Self::cost_of_mul_equals()
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField> PartialEq
    for CubicExtFieldGadget<P, ConstraintF>
{
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1 && self.c2 == other.c2
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField> Eq
    for CubicExtFieldGadget<P, ConstraintF>
{
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    EqGadget<ConstraintF> for CubicExtFieldGadget<P, ConstraintF>
{
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        let b0 = self.c0.is_eq(cs.ns(|| "c0"), &other.c0)?;
        let b1 = self.c1.is_eq(cs.ns(|| "c1"), &other.c1)?;
        let b2 = self.c2.is_eq(cs.ns(|| "c2"), &other.c2)?;
        let temp = Boolean::and(cs.ns(|| "b0 AND b1"), &b0, &b1)?;
        Boolean::and(cs.ns(|| "b0 AND b1 AND b2"), &temp, &b2)
    }

    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.c0
            .conditional_enforce_equal(cs.ns(|| "c0"), &other.c0, should_enforce)?;
        self.c1
            .conditional_enforce_equal(cs.ns(|| "c1"), &other.c1, should_enforce)?;
        self.c2
            .conditional_enforce_equal(cs.ns(|| "c2"), &other.c2, should_enforce)?;
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

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ToBitsGadget<ConstraintF> for CubicExtFieldGadget<P, ConstraintF>
{
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits(&mut cs)?;
        let mut c1 = self.c1.to_bits(&mut cs)?;
        let mut c2 = self.c2.to_bits(cs)?;

        c0.append(&mut c1);
        c0.append(&mut c2);

        Ok(c0)
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits_strict(&mut cs)?;
        let mut c1 = self.c1.to_bits_strict(&mut cs)?;
        let mut c2 = self.c2.to_bits_strict(cs)?;

        c0.append(&mut c1);
        c0.append(&mut c2);

        Ok(c0)
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ToBytesGadget<ConstraintF> for CubicExtFieldGadget<P, ConstraintF>
{
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c0 = self.c0.to_bytes(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_bytes(cs.ns(|| "c1"))?;
        let mut c2 = self.c2.to_bytes(cs.ns(|| "c2"))?;

        c0.append(&mut c1);
        c0.append(&mut c2);

        Ok(c0)
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c0 = self.c0.to_bytes_strict(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_bytes_strict(cs.ns(|| "c1"))?;
        let mut c2 = self.c2.to_bytes_strict(cs.ns(|| "c2"))?;

        c0.append(&mut c1);
        c0.append(&mut c2);

        Ok(c0)
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField> Clone
    for CubicExtFieldGadget<P, ConstraintF>
{
    fn clone(&self) -> Self {
        Self::new(self.c0.clone(), self.c1.clone(), self.c2.clone())
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    CondSelectGadget<ConstraintF> for CubicExtFieldGadget<P, ConstraintF>
{
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = P::BaseFieldGadget::conditionally_select(
            &mut cs.ns(|| "c0"),
            cond,
            &first.c0,
            &second.c0,
        )?;
        let c1 = P::BaseFieldGadget::conditionally_select(
            &mut cs.ns(|| "c1"),
            cond,
            &first.c1,
            &second.c1,
        )?;
        let c2 = P::BaseFieldGadget::conditionally_select(
            &mut cs.ns(|| "c2"),
            cond,
            &first.c2,
            &second.c2,
        )?;

        Ok(Self::new(c0, c1, c2))
    }

    fn cost() -> usize {
        3 * <P::BaseFieldGadget as CondSelectGadget<ConstraintF>>::cost()
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    TwoBitLookupGadget<ConstraintF> for CubicExtFieldGadget<P, ConstraintF>
{
    type TableConstant = CubicExtField<P>;
    fn two_bit_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c2s = c.iter().map(|f| f.c2).collect::<Vec<_>>();
        let c0 = P::BaseFieldGadget::two_bit_lookup(cs.ns(|| "Lookup c0"), b, &c0s)?;
        let c1 = P::BaseFieldGadget::two_bit_lookup(cs.ns(|| "Lookup c1"), b, &c1s)?;
        let c2 = P::BaseFieldGadget::two_bit_lookup(cs.ns(|| "Lookup c2"), b, &c2s)?;
        Ok(Self::new(c0, c1, c2))
    }

    fn two_bit_lookup_lc<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        precomp: &Boolean,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c2s = c.iter().map(|f| f.c2).collect::<Vec<_>>();
        let c0 = P::BaseFieldGadget::two_bit_lookup_lc(cs.ns(|| "Lookup c0"), precomp, b, &c0s)?;
        let c1 = P::BaseFieldGadget::two_bit_lookup_lc(cs.ns(|| "Lookup c1"), precomp, b, &c1s)?;
        let c2 = P::BaseFieldGadget::two_bit_lookup_lc(cs.ns(|| "Lookup c2"), precomp, b, &c2s)?;
        Ok(Self::new(c0, c1, c2))
    }

    fn cost() -> usize {
        3 * <P::BaseFieldGadget as TwoBitLookupGadget<ConstraintF>>::cost()
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ThreeBitCondNegLookupGadget<ConstraintF> for CubicExtFieldGadget<P, ConstraintF>
{
    type TableConstant = CubicExtField<P>;

    fn three_bit_cond_neg_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        b0b1: &Boolean,
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c2s = c.iter().map(|f| f.c2).collect::<Vec<_>>();
        let c0 =
            P::BaseFieldGadget::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c0"), b, b0b1, &c0s)?;
        let c1 =
            P::BaseFieldGadget::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c1"), b, b0b1, &c1s)?;
        let c2 =
            P::BaseFieldGadget::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c2"), b, b0b1, &c2s)?;
        Ok(Self::new(c0, c1, c2))
    }

    fn cost() -> usize {
        3 * <P::BaseFieldGadget as ThreeBitCondNegLookupGadget<ConstraintF>>::cost()
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    AllocGadget<CubicExtField<P>, ConstraintF> for CubicExtFieldGadget<P, ConstraintF>
{
    #[inline]
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<CubicExtField<P>>,
    {
        let (c0, c1, c2) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1), Ok(fe.c2))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = P::BaseFieldGadget::alloc(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = P::BaseFieldGadget::alloc(&mut cs.ns(|| "c1"), || c1)?;
        let c2 = P::BaseFieldGadget::alloc(&mut cs.ns(|| "c2"), || c2)?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<CubicExtField<P>>,
    {
        let (c0, c1, c2) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1), Ok(fe.c2))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = P::BaseFieldGadget::alloc_input(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = P::BaseFieldGadget::alloc_input(&mut cs.ns(|| "c1"), || c1)?;
        let c2 = P::BaseFieldGadget::alloc_input(&mut cs.ns(|| "c2"), || c2)?;
        Ok(Self::new(c0, c1, c2))
    }
}

impl<P: CubicExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ConstantGadget<CubicExtField<P>, ConstraintF> for CubicExtFieldGadget<P, ConstraintF>
{
    #[inline]
    fn from_value<CS: ConstraintSystem<ConstraintF>>(mut cs: CS, value: &CubicExtField<P>) -> Self {
        let c0 = P::BaseFieldGadget::from_value(&mut cs.ns(|| "c0"), &value.c0);
        let c1 = P::BaseFieldGadget::from_value(&mut cs.ns(|| "c1"), &value.c1);
        let c2 = P::BaseFieldGadget::from_value(&mut cs.ns(|| "c2"), &value.c2);
        Self::new(c0, c1, c2)
    }

    #[inline]
    fn get_constant(&self) -> CubicExtField<P> {
        self.get_value().unwrap()
    }
}
