use algebra::{
    fields::fp3::{Fp3, Fp3Parameters},
    PrimeField, SquareRootField,
};
use core::{borrow::Borrow, marker::PhantomData};
use r1cs_core::{ConstraintSystem, ConstraintVar, SynthesisError};

use crate::{fields::fp::FpGadget, prelude::*, Vec};

#[derive(Derivative)]
#[derivative(Debug(
    bound = "P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField"
))]
#[must_use]
pub struct Fp3Gadget<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
{
    pub c0:  FpGadget<ConstraintF>,
    pub c1:  FpGadget<ConstraintF>,
    pub c2:  FpGadget<ConstraintF>,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    Fp3Gadget<P, ConstraintF>
{
    #[inline]
    pub fn new(
        c0: FpGadget<ConstraintF>,
        c1: FpGadget<ConstraintF>,
        c2: FpGadget<ConstraintF>,
    ) -> Self {
        Self {
            c0,
            c1,
            c2,
            _params: PhantomData,
        }
    }

    /// Multiply a FpGadget by quadratic nonresidue P::NONRESIDUE.
    #[inline]
    pub fn mul_fp_gadget_by_nonresidue<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        fe: &FpGadget<ConstraintF>,
    ) -> Result<FpGadget<ConstraintF>, SynthesisError> {
        fe.mul_by_constant(cs, &P::NONRESIDUE)
    }

    /// Multiply a Fp3Gadget by an element of fp.
    #[inline]
    pub fn mul_by_fp_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        fe: &P::Fp,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.mul_by_constant_in_place(cs.ns(|| "c0"), fe)?;
        self.c1.mul_by_constant_in_place(cs.ns(|| "c1"), fe)?;
        self.c2.mul_by_constant_in_place(cs.ns(|| "c2"), fe)?;
        Ok(self)
    }

    /// Multiply a Fp3Gadget by an element of fp.
    #[inline]
    pub fn mul_by_fp_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        fe: &P::Fp,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.mul_by_fp_constant_in_place(cs, fe)?;
        Ok(result)
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    FieldGadget<Fp3<P>, ConstraintF> for Fp3Gadget<P, ConstraintF>
{
    type Variable = (
        ConstraintVar<ConstraintF>,
        ConstraintVar<ConstraintF>,
        ConstraintVar<ConstraintF>,
    );

    #[inline]
    fn get_value(&self) -> Option<Fp3<P>> {
        match (
            self.c0.get_value(),
            self.c1.get_value(),
            self.c2.get_value(),
        ) {
            (Some(c0), Some(c1), Some(c2)) => Some(Fp3::new(c0, c1, c2)),
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
        let c0 = FpGadget::<ConstraintF>::zero(cs.ns(|| "c0"))?;
        let c1 = FpGadget::<ConstraintF>::zero(cs.ns(|| "c1"))?;
        let c2 = FpGadget::<ConstraintF>::zero(cs.ns(|| "c2"))?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn one<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = FpGadget::<ConstraintF>::one(cs.ns(|| "c0"))?;
        let c1 = FpGadget::<ConstraintF>::zero(cs.ns(|| "c1"))?;
        let c2 = FpGadget::<ConstraintF>::zero(cs.ns(|| "c2"))?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn conditionally_add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        bit: &Boolean,
        coeff: Fp3<P>,
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
        // Uses Toom-Cook-3x multiplication from
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

        let two = P::Fp::one().double();
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

    #[inline]
    fn mul_equals<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        // Karatsuba multiplication for Fp3:
        //     v0 = A.c0 * B.c0
        //     v1 = A.c1 * B.c1
        //     v2 = A.c2 * B.c2
        //     result.c0 = v0 + β((a1 + a2)(b1 + b2) − v1 − v2)
        //     result.c1 = (a0 + a1)(b0 + b1) − v0 − v1 + βv2
        //     result.c2 = (a0 + a2)(b0 + b2) − v0 + v1 − v2,
        // We enforce this with six constraints:
        //
        //     v0 = A.c0 * B.c0
        //     v1 = A.c1 * B.c1
        //     v2 = A.c2 * B.c2
        //
        //     result.c0 - v0 + \beta*(v1 + v2) = β(a1 + a2)(b1 + b2))
        //     result.c1 + v0 + v1 - βv2        = (a0 + a1)(b0 + b1)
        //     result.c2 + v0 - v1 + v2         = (a0 + a2)(b0 + b2)
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        // Devegili, OhEigeartaigh, Scott, Dahab
        //
        // This implementation adapted from
        // https://github.com/ZencashOfficial/ginger-lib/blob/development/r1cs/gadgets/std/src/fields/fp3.rs
        let v0 = self.c0.mul(cs.ns(|| "v0 = a0 * b0"), &other.c0)?;
        let v1 = self.c1.mul(cs.ns(|| "v1 = a1 * b1"), &other.c1)?;
        let v2 = self.c2.mul(cs.ns(|| "v2 = a2 * b2"), &other.c2)?;

        // Check c0
        let nr_a1_plus_a2 = self
            .c1
            .add(cs.ns(|| "a1 + a2"), &self.c2)?
            .mul_by_constant(cs.ns(|| "nr*(a1 + a2)"), &P::NONRESIDUE)?;
        let b1_plus_b2 = other.c1.add(cs.ns(|| "b1 + b2"), &other.c2)?;
        let nr_v1 = v1.mul_by_constant(cs.ns(|| "nr * v1"), &P::NONRESIDUE)?;
        let nr_v2 = v2.mul_by_constant(cs.ns(|| "nr * v2"), &P::NONRESIDUE)?;
        let to_check = result
            .c0
            .sub(cs.ns(|| "c0 - v0"), &v0)?
            .add(cs.ns(|| "c0 - v0 + nr * v1"), &nr_v1)?
            .add(cs.ns(|| "c0 - v0 + nr * v1 + nr * v2"), &nr_v2)?;
        nr_a1_plus_a2.mul_equals(cs.ns(|| "check c0"), &b1_plus_b2, &to_check)?;

        // Check c1
        let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;
        let b0_plus_b1 = other.c0.add(cs.ns(|| "b0 + b1"), &other.c1)?;
        let to_check = result
            .c1
            .sub(cs.ns(|| "c1 - nr * v2"), &nr_v2)?
            .add(cs.ns(|| "c1 - nr * v2 + v0"), &v0)?
            .add(cs.ns(|| "c1 - nr * v2 + v0 + v1"), &v1)?;
        a0_plus_a1.mul_equals(cs.ns(|| "check c1"), &b0_plus_b1, &to_check)?;

        // Check c2
        let a0_plus_a2 = self.c0.add(cs.ns(|| "a0 + a2"), &self.c2)?;
        let b0_plus_b2 = other.c0.add(cs.ns(|| "b0 + b2"), &other.c2)?;
        let to_check = result
            .c2
            .add(cs.ns(|| "c2 + v0"), &v0)?
            .sub(cs.ns(|| "c2 + v0 - v1"), &v1)?
            .add(cs.ns(|| "c2 + v0 - v1 + v2"), &v2)?;
        a0_plus_a2.mul_equals(cs.ns(|| "check c2"), &b0_plus_b2, &to_check)?;
        Ok(())
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

        let c0 = Self::mul_fp_gadget_by_nonresidue(cs.ns(|| "c0 part 1"), &s3)?
            .add(cs.ns(|| "c0"), &s0)?;

        let c1 = Self::mul_fp_gadget_by_nonresidue(cs.ns(|| "c1 part 1"), &s4)?
            .add(cs.ns(|| "c1"), &s1)?;

        let c2 = s1
            .add(cs.ns(|| "c2 part1"), &s2)?
            .add(cs.ns(|| "c2 part2"), &s3)?
            .sub(cs.ns(|| "c2 part3"), &s0)?
            .sub(cs.ns(|| "c2 part4"), &s4)?;

        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Fp3<P>,
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
        other: &Fp3<P>,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.add_constant_in_place(cs.ns(|| "c0"), &other.c0)?;
        self.c1.add_constant_in_place(cs.ns(|| "c1"), &other.c1)?;
        self.c2.add_constant_in_place(cs.ns(|| "c2"), &other.c2)?;
        Ok(self)
    }

    /// Use the Toom-Cook-3x method to compute multiplication.
    #[inline]
    fn mul_by_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Fp3<P>,
    ) -> Result<Self, SynthesisError> {
        // Uses Toom-Cook-3x multiplication from
        //
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        //    Devegili, OhEigeartaigh, Scott, Dahab

        // v0 = a(0)b(0)   = a0 * b0
        let v0 = self.c0.mul_by_constant(cs.ns(|| "v0"), &other.c0)?;

        // v1 = a(1)b(1)   = (a0 + a1 + a2)(b0 + b1 + b2)
        let v1 = {
            let mut v1_cs = cs.ns(|| "v1");
            let mut a0_plus_a1_plus_a2 = self
                .c0
                .add(v1_cs.ns(|| "a0 + a1"), &self.c1)?
                .add(v1_cs.ns(|| "a0 + a1 + a2"), &self.c2)?;
            let b0_plus_b1_plus_b2 = other.c0 + &other.c1 + &other.c2;

            a0_plus_a1_plus_a2.mul_by_constant_in_place(
                v1_cs.ns(|| "(a0 + a1 + a2)*(b0 + b1 + b2)"),
                &b0_plus_b1_plus_b2,
            )?;
            a0_plus_a1_plus_a2
        };

        // v2 = a(−1)b(−1) = (a0 − a1 + a2)(b0 − b1 + b2)
        let mut v2 = {
            let mut v2_cs = cs.ns(|| "v2");
            let mut a0_minus_a1_plus_a2 = self
                .c0
                .sub(v2_cs.ns(|| "sub1"), &self.c1)?
                .add(v2_cs.ns(|| "add2"), &self.c2)?;
            let b0_minus_b1_plus_b2 = other.c0 - &other.c1 + &other.c2;
            a0_minus_a1_plus_a2.mul_by_constant_in_place(
                v2_cs.ns(|| "(a0 - a1 + a2)*(b0 - b1 + b2)"),
                &b0_minus_b1_plus_b2,
            )?;
            a0_minus_a1_plus_a2
        };

        // v3 = a(2)b(2)   = (a0 + 2a1 + 4a2)(b0 + 2b1 + 4b2)
        let mut v3 = {
            let mut v3_cs = cs.ns(|| "v3");
            let a1_double = self.c1.double(v3_cs.ns(|| "2a1"))?;
            let a2_quad = self
                .c2
                .double(v3_cs.ns(|| "2a2"))?
                .double(v3_cs.ns(|| "4a2"))?;
            let mut a0_plus_2_a1_plus_4_a2 = self
                .c0
                .add(v3_cs.ns(|| "a0 + 2a1"), &a1_double)?
                .add(v3_cs.ns(|| "a0 + 2a1 + 4a2"), &a2_quad)?;

            let b1_double = other.c1.double();
            let b2_quad = other.c2.double().double();
            let b0_plus_2_b1_plus_4_b2 = other.c0 + &b1_double + &b2_quad;

            a0_plus_2_a1_plus_4_a2.mul_by_constant_in_place(
                v3_cs.ns(|| "(a0 + 2a1 + 4a2)*(b0 + 2b1 + 4b2)"),
                &b0_plus_2_b1_plus_4_b2,
            )?;
            a0_plus_2_a1_plus_4_a2
        };

        // v4 = a(∞)b(∞)   = a2 * b2
        let v4 = self.c2.mul_by_constant(cs.ns(|| "v4"), &other.c2)?;

        let two = P::Fp::one().double();
        let six = two.double() + &two;
        let mut two_and_six = [two, six];
        algebra::fields::batch_inversion(&mut two_and_six);
        let (two_inverse, six_inverse) = (two_and_six[0], two_and_six[1]);

        let mut half_v0 = v0.mul_by_constant(cs.ns(|| "half_v0"), &two_inverse)?;
        let half_v1 = v1.mul_by_constant(cs.ns(|| "half_v1"), &two_inverse)?;
        let mut one_sixth_v2 = v2.mul_by_constant(cs.ns(|| "v2_by_6"), &six_inverse)?;
        let one_sixth_v3 = v3.mul_by_constant_in_place(cs.ns(|| "v3_by_6"), &six_inverse)?;
        let two_v4 = v4.double(cs.ns(|| "2v4"))?;

        // c0 = v0 + β((1/2)v0 − (1/2)v1 − (1/6)v2 + (1/6)v3 − 2v4)
        let c0 = {
            let mut c0_cs = cs.ns(|| "c0");

            // No constraints, only get a linear combination back.
            let mut inner = half_v0
                .sub(c0_cs.ns(|| "sub1"), &half_v1)?
                .sub(c0_cs.ns(|| "sub2"), &one_sixth_v2)?
                .add(c0_cs.ns(|| "add3"), &one_sixth_v3)?
                .sub(c0_cs.ns(|| "sub4"), &two_v4)?;
            let non_residue_times_inner =
                inner.mul_by_constant_in_place(&mut c0_cs, &P::NONRESIDUE)?;
            v0.add(c0_cs.ns(|| "add5"), non_residue_times_inner)?
        };

        // −(1/2)v0 + v1 − (1/3)v2 − (1/6)v3 + 2v4 + βv4
        let c1 = {
            let mut c1_cs = cs.ns(|| "c1");
            let one_third_v2 = one_sixth_v2.double_in_place(c1_cs.ns(|| "double1"))?;
            let non_residue_v4 =
                v4.mul_by_constant(c1_cs.ns(|| "mul_by_const1"), &P::NONRESIDUE)?;

            half_v0
                .negate_in_place(c1_cs.ns(|| "neg1"))?
                .add(c1_cs.ns(|| "add1"), &v1)?
                .sub(c1_cs.ns(|| "sub2"), one_third_v2)?
                .sub(c1_cs.ns(|| "sub3"), &one_sixth_v3)?
                .add(c1_cs.ns(|| "add4"), &two_v4)?
                .add(c1_cs.ns(|| "add5"), &non_residue_v4)?
        };

        // -v0 + (1/2)v1 + (1/2)v2 −v4
        let c2 = {
            let mut c2_cs = cs.ns(|| "c2");
            let half_v2 = v2.mul_by_constant_in_place(c2_cs.ns(|| "half_v2"), &two_inverse)?;
            half_v1
                .add(c2_cs.ns(|| "add1"), half_v2)?
                .sub(c2_cs.ns(|| "sub2"), &v4)?
                .sub(c2_cs.ns(|| "sub3"), &v0)?
        };

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
        self.c1.mul_by_constant_in_place(
            cs.ns(|| "c1_power"),
            &P::FROBENIUS_COEFF_FP3_C1[power % 3],
        )?;
        self.c2.mul_by_constant_in_place(
            cs.ns(|| "c2_power"),
            &P::FROBENIUS_COEFF_FP3_C2[power % 3],
        )?;

        Ok(self)
    }

    fn cost_of_mul() -> usize {
        5 * FpGadget::<ConstraintF>::cost_of_mul()
    }

    fn cost_of_mul_equals() -> usize {
        6 * FpGadget::<ConstraintF>::cost_of_mul()
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> PartialEq
    for Fp3Gadget<P, ConstraintF>
{
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1 && self.c2 == other.c2
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> Eq
    for Fp3Gadget<P, ConstraintF>
{
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    EqGadget<ConstraintF> for Fp3Gadget<P, ConstraintF>
{
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ConditionalEqGadget<ConstraintF> for Fp3Gadget<P, ConstraintF>
{
    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.c0
            .conditional_enforce_equal(&mut cs.ns(|| "c0"), &other.c0, condition)?;
        self.c1
            .conditional_enforce_equal(&mut cs.ns(|| "c1"), &other.c1, condition)?;
        self.c2
            .conditional_enforce_equal(&mut cs.ns(|| "c2"), &other.c2, condition)?;
        Ok(())
    }

    fn cost() -> usize {
        3 * <FpGadget<ConstraintF> as ConditionalEqGadget<ConstraintF>>::cost()
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    NEqGadget<ConstraintF> for Fp3Gadget<P, ConstraintF>
{
    #[inline]
    fn enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        self.c0.enforce_not_equal(&mut cs.ns(|| "c0"), &other.c0)?;
        self.c1.enforce_not_equal(&mut cs.ns(|| "c1"), &other.c1)?;
        self.c2.enforce_not_equal(&mut cs.ns(|| "c2"), &other.c2)?;
        Ok(())
    }

    fn cost() -> usize {
        3 * <FpGadget<ConstraintF> as NEqGadget<ConstraintF>>::cost()
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ToBitsGadget<ConstraintF> for Fp3Gadget<P, ConstraintF>
{
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_bits(cs.ns(|| "c1"))?;
        let mut c2 = self.c2.to_bits(cs.ns(|| "c2"))?;

        c0.append(&mut c1);
        c0.append(&mut c2);

        Ok(c0)
    }

    fn to_non_unique_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_non_unique_bits(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_non_unique_bits(cs.ns(|| "c1"))?;
        let mut c2 = self.c2.to_non_unique_bits(cs.ns(|| "c2"))?;

        c0.append(&mut c1);
        c0.append(&mut c2);

        Ok(c0)
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ToBytesGadget<ConstraintF> for Fp3Gadget<P, ConstraintF>
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

    fn to_non_unique_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c0 = self.c0.to_non_unique_bytes(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_non_unique_bytes(cs.ns(|| "c1"))?;
        let mut c2 = self.c2.to_non_unique_bytes(cs.ns(|| "c2"))?;

        c0.append(&mut c1);
        c0.append(&mut c2);

        Ok(c0)
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField> Clone
    for Fp3Gadget<P, ConstraintF>
{
    fn clone(&self) -> Self {
        Self::new(self.c0.clone(), self.c1.clone(), self.c2.clone())
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    CondSelectGadget<ConstraintF> for Fp3Gadget<P, ConstraintF>
{
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = FpGadget::<ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c0"),
            cond,
            &first.c0,
            &second.c0,
        )?;
        let c1 = FpGadget::<ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c1"),
            cond,
            &first.c1,
            &second.c1,
        )?;
        let c2 = FpGadget::<ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c2"),
            cond,
            &first.c2,
            &second.c2,
        )?;

        Ok(Self::new(c0, c1, c2))
    }

    fn cost() -> usize {
        3 * <FpGadget<ConstraintF> as CondSelectGadget<ConstraintF>>::cost()
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    TwoBitLookupGadget<ConstraintF> for Fp3Gadget<P, ConstraintF>
{
    type TableConstant = Fp3<P>;
    fn two_bit_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c2s = c.iter().map(|f| f.c2).collect::<Vec<_>>();
        let c0 = FpGadget::<ConstraintF>::two_bit_lookup(cs.ns(|| "Lookup c0"), b, &c0s)?;
        let c1 = FpGadget::<ConstraintF>::two_bit_lookup(cs.ns(|| "Lookup c1"), b, &c1s)?;
        let c2 = FpGadget::<ConstraintF>::two_bit_lookup(cs.ns(|| "Lookup c2"), b, &c2s)?;
        Ok(Self::new(c0, c1, c2))
    }

    fn cost() -> usize {
        3 * <FpGadget<ConstraintF> as TwoBitLookupGadget<ConstraintF>>::cost()
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ThreeBitCondNegLookupGadget<ConstraintF> for Fp3Gadget<P, ConstraintF>
{
    type TableConstant = Fp3<P>;

    fn three_bit_cond_neg_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        b0b1: &Boolean,
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c2s = c.iter().map(|f| f.c2).collect::<Vec<_>>();
        let c0 = FpGadget::<ConstraintF>::three_bit_cond_neg_lookup(
            cs.ns(|| "Lookup c0"),
            b,
            b0b1,
            &c0s,
        )?;
        let c1 = FpGadget::<ConstraintF>::three_bit_cond_neg_lookup(
            cs.ns(|| "Lookup c1"),
            b,
            b0b1,
            &c1s,
        )?;
        let c2 = FpGadget::<ConstraintF>::three_bit_cond_neg_lookup(
            cs.ns(|| "Lookup c2"),
            b,
            b0b1,
            &c2s,
        )?;
        Ok(Self::new(c0, c1, c2))
    }

    fn cost() -> usize {
        3 * <FpGadget<ConstraintF> as ThreeBitCondNegLookupGadget<ConstraintF>>::cost()
    }
}

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    AllocGadget<Fp3<P>, ConstraintF> for Fp3Gadget<P, ConstraintF>
{
    #[inline]
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Fp3<P>>,
    {
        let (c0, c1, c2) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1), Ok(fe.c2))
            },
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = FpGadget::<ConstraintF>::alloc(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = FpGadget::<ConstraintF>::alloc(&mut cs.ns(|| "c1"), || c1)?;
        let c2 = FpGadget::<ConstraintF>::alloc(&mut cs.ns(|| "c2"), || c2)?;
        Ok(Self::new(c0, c1, c2))
    }

    #[inline]
    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Fp3<P>>,
    {
        let (c0, c1, c2) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1), Ok(fe.c2))
            },
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = FpGadget::<ConstraintF>::alloc_input(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = FpGadget::<ConstraintF>::alloc_input(&mut cs.ns(|| "c1"), || c1)?;
        let c2 = FpGadget::<ConstraintF>::alloc_input(&mut cs.ns(|| "c2"), || c2)?;
        Ok(Self::new(c0, c1, c2))
    }
}
