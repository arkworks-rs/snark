use algebra::{
    biginteger::arithmetic::find_wnaf, Field, PrimeField, QuadExtField, QuadExtParameters,
    SquareRootField,
};
use r1cs_core::{ConstraintSystem, SynthesisError};
use std::{borrow::Borrow, marker::PhantomData};

use crate::{fields::FieldGadget, prelude::*};

pub trait QuadExtParametersGadget<ConstraintF: PrimeField>:
    QuadExtParameters<BasePrimeField = ConstraintF>
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
        power: usize,
    ) -> Result<(), SynthesisError>;

    /// Compute the cyclotomic square of fe, which must be in the cyclotomic subgroup.
    fn cyclotomic_square_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        fe: &QuadExtFieldGadget<Self, ConstraintF>,
    ) -> Result<QuadExtFieldGadget<Self, ConstraintF>, SynthesisError>
    where
        ConstraintF: PrimeField + SquareRootField,
    {
        fe.square(cs)
    }
}

#[derive(Derivative)]
#[derivative(Debug(
    bound = "P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField"
))]
#[must_use]
pub struct QuadExtFieldGadget<
    P: QuadExtParametersGadget<ConstraintF>,
    ConstraintF: PrimeField + SquareRootField,
> {
    pub c0: P::BaseFieldGadget,
    pub c1: P::BaseFieldGadget,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    QuadExtFieldGadget<P, ConstraintF>
{
    pub fn new(c0: P::BaseFieldGadget, c1: P::BaseFieldGadget) -> Self {
        Self {
            c0,
            c1,
            _params: PhantomData,
        }
    }

    #[inline]
    pub fn unitary_inverse<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
        let new_c0 = self.c0.clone();
        let new_c1 = self.c1.clone().negate(cs.ns(|| "c1 negation"))?;
        Ok(Self::new(new_c0, new_c1))
    }

    #[inline]
    pub fn conjugate_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.c1.negate_in_place(cs)?;
        Ok(self)
    }

    #[inline]
    pub fn cyclotomic_exp<CS: ConstraintSystem<ConstraintF>, S: AsRef<[u64]>>(
        &self,
        mut cs: CS,
        exp: S,
    ) -> Result<Self, SynthesisError> {
        let mut res = Self::one(cs.ns(|| "one"))?;
        let self_inverse = self.unitary_inverse(cs.ns(|| "unitary inverse"))?;
        let mut found_nonzero = false;
        let naf = find_wnaf(exp.as_ref());

        for (j, &value) in naf.iter().rev().enumerate() {
            if found_nonzero {
                res = P::cyclotomic_square_gadget(cs.ns(|| format!("res_square_{:?}", j)), &res)?;
            }
            if value != 0 {
                found_nonzero = true;

                if value > 0 {
                    res.mul_in_place(cs.ns(|| format!("res_mul_{:?}", j)), self)?;
                } else {
                    res.mul_in_place(cs.ns(|| format!("res_mul_inverse_{:?}", j)), &self_inverse)?;
                }
            }
        }
        Ok(res)
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    FieldGadget<QuadExtField<P>, ConstraintF> for QuadExtFieldGadget<P, ConstraintF>
{
    type Variable = (
        <P::BaseFieldGadget as FieldGadget<P::BaseField, ConstraintF>>::Variable,
        <P::BaseFieldGadget as FieldGadget<P::BaseField, ConstraintF>>::Variable,
    );

    #[inline]
    fn get_value(&self) -> Option<QuadExtField<P>> {
        match (self.c0.get_value(), self.c1.get_value()) {
            (Some(c0), Some(c1)) => Some(QuadExtField::<P>::new(c0, c1)),
            (..) => None,
        }
    }

    #[inline]
    fn get_variable(&self) -> Self::Variable {
        (
            self.c0.get_variable().clone(),
            self.c1.get_variable().clone(),
        )
    }

    #[inline]
    fn zero<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = P::BaseFieldGadget::zero(cs.ns(|| "c0"))?;
        let c1 = P::BaseFieldGadget::zero(cs.ns(|| "c1"))?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn one<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = P::BaseFieldGadget::one(cs.ns(|| "c0"))?;
        let c1 = P::BaseFieldGadget::zero(cs.ns(|| "c1"))?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn conditionally_add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        bit: &Boolean,
        coeff: QuadExtField<P>,
    ) -> Result<Self, SynthesisError> {
        let c0 = self
            .c0
            .conditionally_add_constant(cs.ns(|| "c0"), bit, coeff.c0)?;
        let c1 = self
            .c1
            .conditionally_add_constant(cs.ns(|| "c1"), bit, coeff.c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn add<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = self.c0.add(&mut cs.ns(|| "add c0"), &other.c0)?;
        let c1 = self.c1.add(&mut cs.ns(|| "add c1"), &other.c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn double<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.double_in_place(cs)?;
        Ok(result)
    }

    #[inline]
    fn double_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.double_in_place(&mut cs.ns(|| "double c0"))?;
        self.c1.double_in_place(&mut cs.ns(|| "double c1"))?;
        Ok(self)
    }

    #[inline]
    fn sub<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = self.c0.sub(&mut cs.ns(|| "sub c0"), &other.c0)?;
        let c1 = self.c1.sub(&mut cs.ns(|| "sub c1"), &other.c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn negate<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.negate_in_place(cs)?;
        Ok(result)
    }

    #[inline]
    fn negate_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.negate_in_place(&mut cs.ns(|| "negate c0"))?;
        self.c1.negate_in_place(&mut cs.ns(|| "negate c1"))?;
        Ok(self)
    }

    #[inline]
    fn mul<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        // Karatsuba multiplication for Fp2:
        //     v0 = A.c0 * B.c0
        //     v1 = A.c1 * B.c1
        //     result.c0 = v0 + non_residue * v1
        //     result.c1 = (A.c0 + A.c1) * (B.c0 + B.c1) - v0 - v1
        // Enforced with 3 constraints:
        //     A.c1 * B.c1 = v1
        //     A.c0 * B.c0 = result.c0 - non_residue * v1
        //     (A.c0+A.c1)*(B.c0+B.c1) = result.c1 + result.c0 + (1 - non_residue) * v1
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        // Devegili, OhEigeartaigh, Scott, Dahab
        let mul_cs = &mut cs.ns(|| "mul");

        let v0 = self.c0.mul(mul_cs.ns(|| "v0"), &other.c0)?;
        let v1 = self.c1.mul(mul_cs.ns(|| "v1"), &other.c1)?;
        let c0 = {
            let non_residue_times_v1 =
                v1.mul_by_constant(mul_cs.ns(|| "non_residue * v0"), &P::NONRESIDUE)?;
            v0.add(mul_cs.ns(|| "v0 + beta * v1"), &non_residue_times_v1)?
        };
        let c1 = {
            let a0_plus_a1 = self.c0.add(mul_cs.ns(|| "a0 + a1"), &self.c1)?;
            let b0_plus_b1 = other.c0.add(mul_cs.ns(|| "b0 + b1"), &other.c1)?;
            let a0_plus_a1_times_b0_plus_b1 =
                a0_plus_a1.mul(&mut mul_cs.ns(|| "(a0 + a1) * (b0 + b1)"), &b0_plus_b1)?;
            a0_plus_a1_times_b0_plus_b1
                .sub(mul_cs.ns(|| "res - v0"), &v0)?
                .sub(mul_cs.ns(|| "res - v0 - v1"), &v1)?
        };
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn square<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Self, SynthesisError> {
        // From Libsnark/fp2_gadget.tcc
        // Complex multiplication for Fp2:
        //     v0 = A.c0 * A.c1
        //     result.c0 = (A.c0 + A.c1) * (A.c0 + non_residue * A.c1) - (1 +
        // non_residue) * v0     result.c1 = 2 * v0
        // Enforced with 2 constraints:
        //     (2*A.c0) * A.c1 = result.c1
        //     (A.c0 + A.c1) * (A.c0 + non_residue * A.c1) = result.c0 + result.c1 * (1
        // + non_residue)/2 Reference:
        //     "Multiplication and Squaring on Pairing-Friendly Fields"
        //     Devegili, OhEigeartaigh, Scott, Dahab

        let mut v0 = self.c0.mul(cs.ns(|| "v0"), &self.c1)?;
        let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;

        let non_residue_c1 = self
            .c1
            .mul_by_constant(cs.ns(|| "non_residue * a1"), &P::NONRESIDUE)?;
        let a0_plus_non_residue_c1 = self
            .c0
            .add(cs.ns(|| "a0 + non_residue * a1"), &non_residue_c1)?;
        let one_plus_non_residue_v0 = v0.mul_by_constant(
            cs.ns(|| "1 + non_residue * v0"),
            &(P::BaseField::one() + &P::NONRESIDUE),
        )?;

        let c0 = a0_plus_a1
            .mul(
                cs.ns(|| "(a0 + a1) * (a0 + non_residue * a1)"),
                &a0_plus_non_residue_c1,
            )?
            .sub(cs.ns(|| "- (1 + non_residue) v0"), &one_plus_non_residue_v0)?;

        v0.double_in_place(cs.ns(|| "2v0"))?;
        let c1 = v0;

        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn square_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        // From Libsnark/fp2_gadget.tcc
        // Complex multiplication for Fp2:
        //     v0 = A.c0 * A.c1
        //     result.c0 = (A.c0 + A.c1) * (A.c0 + non_residue * A.c1) - (1 +
        // non_residue) * v0     result.c1 = 2 * v0
        // Enforced with 2 constraints:
        //     (2*A.c0) * A.c1 = result.c1
        //     (A.c0 + A.c1) * (A.c0 + non_residue * A.c1) = result.c0 + result.c1 * (1
        // + non_residue)/2 Reference:
        //     "Multiplication and Squaring on Pairing-Friendly Fields"
        //     Devegili, OhEigeartaigh, Scott, Dahab

        let mut v0 = self.c0.mul(cs.ns(|| "v0"), &self.c1)?;
        let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;

        let _ = self
            .c1
            .mul_by_constant_in_place(cs.ns(|| "non_residue * a1"), &P::NONRESIDUE)?;
        let a0_plus_non_residue_c1 = self.c0.add(cs.ns(|| "a0 + non_residue * a1"), &self.c1)?;
        let one_plus_non_residue_v0 = v0.mul_by_constant(
            cs.ns(|| "1 + non_residue * v0"),
            &(P::BaseField::one() + &P::NONRESIDUE),
        )?;

        self.c0 = a0_plus_a1
            .mul(
                cs.ns(|| "(a0 + a1) * (a0 + non_residue * a1)"),
                &a0_plus_non_residue_c1,
            )?
            .sub(cs.ns(|| "- (1 + non_residue) v0"), &one_plus_non_residue_v0)?;

        v0.double_in_place(cs.ns(|| "2v0"))?;
        self.c1 = v0;

        Ok(self)
    }

    #[inline]
    fn mul_equals<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        // Karatsuba multiplication for Fp2:
        //     v0 = A.c0 * B.c0
        //     v1 = A.c1 * B.c1
        //     result.c0 = v0 + non_residue * v1
        //     result.c1 = (A.c0 + A.c1) * (B.c0 + B.c1) - v0 - v1
        // Enforced with 3 constraints:
        //     A.c1 * B.c1 = v1
        //     A.c0 * B.c0 = result.c0 - non_residue * v1
        //     (A.c0+A.c1)*(B.c0+B.c1) = result.c1 + result.c0 + (1 - non_residue) * v1
        // Reference:
        // "Multiplication and Squaring on Pairing-Friendly Fields"
        // Devegili, OhEigeartaigh, Scott, Dahab
        let mul_cs = &mut cs.ns(|| "mul");

        // Compute v1
        let mut v1 = self.c1.mul(mul_cs.ns(|| "v1"), &other.c1)?;

        // Perform second check
        let non_residue_times_v1 =
            v1.mul_by_constant(mul_cs.ns(|| "non_residue * v0"), &P::NONRESIDUE)?;
        let rhs = result
            .c0
            .sub(mul_cs.ns(|| "sub from result.c0"), &non_residue_times_v1)?;
        self.c0
            .mul_equals(mul_cs.ns(|| "second check"), &other.c0, &rhs)?;

        // Last check
        let a0_plus_a1 = self.c0.add(mul_cs.ns(|| "a0 + a1"), &self.c1)?;
        let b0_plus_b1 = other.c0.add(mul_cs.ns(|| "b0 + b1"), &other.c1)?;
        let one_minus_non_residue_v1 =
            v1.sub_in_place(mul_cs.ns(|| "sub from v1"), &non_residue_times_v1)?;

        let result_c1_plus_result_c0_plus_one_minus_non_residue_v1 = result
            .c1
            .add(mul_cs.ns(|| "c1 + c0"), &result.c0)?
            .add(mul_cs.ns(|| "rest of stuff"), one_minus_non_residue_v1)?;

        a0_plus_a1.mul_equals(
            mul_cs.ns(|| "third check"),
            &b0_plus_b1,
            &result_c1_plus_result_c0_plus_one_minus_non_residue_v1,
        )?;

        Ok(())
    }

    #[inline]
    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &QuadExtField<P>,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        let _ = result.add_constant_in_place(cs, other)?;
        Ok(result)
    }

    #[inline]
    fn add_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        other: &QuadExtField<P>,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.add_constant_in_place(cs.ns(|| "c0"), &other.c0)?;
        self.c1.add_constant_in_place(cs.ns(|| "c1"), &other.c1)?;
        Ok(self)
    }

    fn mul_by_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        fe: &QuadExtField<P>,
    ) -> Result<Self, SynthesisError> {
        // Karatsuba multiplication (see mul above).
        // Doesn't need any constraints; returns linear combinations of
        // `self`'s variables.
        //
        // (The operations below are guaranteed to return linear combinations)
        let (a0, a1) = (&self.c0, &self.c1);
        let (b0, b1) = (fe.c0, fe.c1);
        let mut v0 = a0.mul_by_constant(&mut cs.ns(|| "v0"), &b0)?;
        let beta_v1 = a1.mul_by_constant(&mut cs.ns(|| "v1"), &(b1 * &P::NONRESIDUE))?;

        v0.add_in_place(&mut cs.ns(|| "c0"), &beta_v1)?;
        let c0 = v0;

        let mut a0b1 = a0.mul_by_constant(&mut cs.ns(|| "a0b1"), &b1)?;
        let a1b0 = a1.mul_by_constant(&mut cs.ns(|| "a1b0"), &b0)?;
        a0b1.add_in_place(&mut cs.ns(|| "c1"), &a1b0)?;
        let c1 = a0b1;
        Ok(Self::new(c0, c1))
    }

    fn frobenius_map<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        power: usize,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        let _ = result.frobenius_map_in_place(cs, power)?;
        Ok(result)
    }

    fn frobenius_map_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        power: usize,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.frobenius_map_in_place(&mut cs.ns(|| "c0"), power)?;
        self.c1.frobenius_map_in_place(&mut cs.ns(|| "c1"), power)?;

        P::mul_base_field_gadget_by_frobenius_coeff(
            &mut cs.ns(|| "c1_power"),
            &mut self.c1,
            power,
        )?;
        Ok(self)
    }

    fn cost_of_mul() -> usize {
        3
    }

    fn cost_of_mul_equals() -> usize {
        3
    }

    fn cost_of_inv() -> usize {
        Self::cost_of_mul_equals()
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField> PartialEq
    for QuadExtFieldGadget<P, ConstraintF>
{
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField> Eq
    for QuadExtFieldGadget<P, ConstraintF>
{
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    EqGadget<ConstraintF> for QuadExtFieldGadget<P, ConstraintF>
{
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        let b0 = self.c0.is_eq(cs.ns(|| "c0"), &other.c0)?;
        let b1 = self.c1.is_eq(cs.ns(|| "c1"), &other.c1)?;
        Boolean::and(cs.ns(|| "b0 AND b1"), &b0, &b1)
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

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ToBitsGadget<ConstraintF> for QuadExtFieldGadget<P, ConstraintF>
{
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits(&mut cs)?;
        let mut c1 = self.c1.to_bits(cs)?;
        c0.append(&mut c1);
        Ok(c0)
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits_strict(&mut cs)?;
        let mut c1 = self.c1.to_bits_strict(cs)?;
        c0.append(&mut c1);
        Ok(c0)
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ToBytesGadget<ConstraintF> for QuadExtFieldGadget<P, ConstraintF>
{
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c0 = self.c0.to_bytes(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_bytes(cs.ns(|| "c1"))?;
        c0.append(&mut c1);
        Ok(c0)
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c0 = self.c0.to_bytes_strict(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_bytes_strict(cs.ns(|| "c1"))?;
        c0.append(&mut c1);
        Ok(c0)
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField> Clone
    for QuadExtFieldGadget<P, ConstraintF>
{
    fn clone(&self) -> Self {
        Self {
            c0: self.c0.clone(),
            c1: self.c1.clone(),
            _params: PhantomData,
        }
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    CondSelectGadget<ConstraintF> for QuadExtFieldGadget<P, ConstraintF>
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

        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    TwoBitLookupGadget<ConstraintF> for QuadExtFieldGadget<P, ConstraintF>
{
    type TableConstant = QuadExtField<P>;
    fn two_bit_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 = P::BaseFieldGadget::two_bit_lookup(cs.ns(|| "Lookup c0"), b, &c0s)?;
        let c1 = P::BaseFieldGadget::two_bit_lookup(cs.ns(|| "Lookup c1"), b, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn two_bit_lookup_lc<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        precomp: &Boolean,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 = P::BaseFieldGadget::two_bit_lookup_lc(cs.ns(|| "Lookup c0"), precomp, b, &c0s)?;
        let c1 = P::BaseFieldGadget::two_bit_lookup_lc(cs.ns(|| "Lookup c1"), precomp, b, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <P::BaseFieldGadget as TwoBitLookupGadget<ConstraintF>>::cost()
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ThreeBitCondNegLookupGadget<ConstraintF> for QuadExtFieldGadget<P, ConstraintF>
{
    type TableConstant = QuadExtField<P>;

    fn three_bit_cond_neg_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        b0b1: &Boolean,
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 =
            P::BaseFieldGadget::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c0"), b, b0b1, &c0s)?;
        let c1 =
            P::BaseFieldGadget::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c1"), b, b0b1, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <P::BaseFieldGadget as ThreeBitCondNegLookupGadget<ConstraintF>>::cost()
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    AllocGadget<QuadExtField<P>, ConstraintF> for QuadExtFieldGadget<P, ConstraintF>
{
    #[inline]
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<QuadExtField<P>>,
    {
        let (c0, c1) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1))
            }
            Err(_) => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = P::BaseFieldGadget::alloc(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = P::BaseFieldGadget::alloc(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<QuadExtField<P>>,
    {
        let (c0, c1) = match value_gen() {
            Ok(fe) => {
                let fe = *fe.borrow();
                (Ok(fe.c0), Ok(fe.c1))
            }
            Err(_) => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let c0 = P::BaseFieldGadget::alloc_input(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = P::BaseFieldGadget::alloc_input(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }
}

impl<P: QuadExtParametersGadget<ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    ConstantGadget<QuadExtField<P>, ConstraintF> for QuadExtFieldGadget<P, ConstraintF>
{
    #[inline]
    fn from_value<CS: ConstraintSystem<ConstraintF>>(mut cs: CS, value: &QuadExtField<P>) -> Self {
        let c0 = P::BaseFieldGadget::from_value(&mut cs.ns(|| "c0"), &value.c0);
        let c1 = P::BaseFieldGadget::from_value(&mut cs.ns(|| "c1"), &value.c1);
        Self::new(c0, c1)
    }

    #[inline]
    fn get_constant(&self) -> QuadExtField<P> {
        self.get_value().unwrap()
    }
}
