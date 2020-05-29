use algebra::{
    fields::{Fp2, Fp2Parameters},
    PrimeField,
};
use core::{borrow::Borrow, marker::PhantomData};
use r1cs_core::{ConstraintSystem, ConstraintVar, SynthesisError};

use crate::{fields::fp::FpGadget, prelude::*, Vec};

#[derive(Derivative)]
#[derivative(Debug(bound = "P: Fp2Parameters, ConstraintF: PrimeField"))]
#[must_use]
pub struct Fp2Gadget<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> {
    pub c0: FpGadget<ConstraintF>,
    pub c1: FpGadget<ConstraintF>,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> Fp2Gadget<P, ConstraintF> {
    pub fn new(c0: FpGadget<ConstraintF>, c1: FpGadget<ConstraintF>) -> Self {
        Self {
            c0,
            c1,
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

    /// Multiply a Fp2Gadget by an element of fp.
    #[inline]
    pub fn mul_by_fp_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        fe: &P::Fp,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.mul_by_constant_in_place(cs.ns(|| "c0"), fe)?;
        self.c1.mul_by_constant_in_place(cs.ns(|| "c1"), fe)?;
        Ok(self)
    }

    /// Multiply a Fp2Gadget by an element of fp.
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

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> FieldGadget<Fp2<P>, ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    type Variable = (ConstraintVar<ConstraintF>, ConstraintVar<ConstraintF>);

    #[inline]
    fn get_value(&self) -> Option<Fp2<P>> {
        match (self.c0.value, self.c1.value) {
            (Some(c0), Some(c1)) => Some(Fp2::new(c0, c1)),
            (..) => None,
        }
    }

    #[inline]
    fn get_variable(&self) -> Self::Variable {
        (self.c0.get_variable(), self.c1.get_variable())
    }

    #[inline]
    fn zero<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = FpGadget::zero(cs.ns(|| "c0"))?;
        let c1 = FpGadget::zero(cs.ns(|| "c1"))?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn one<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = FpGadget::one(cs.ns(|| "c0"))?;
        let c1 = FpGadget::zero(cs.ns(|| "c1"))?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn conditionally_add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        bit: &Boolean,
        coeff: Fp2<P>,
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
            &(P::Fp::one() + &P::NONRESIDUE),
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
            &(P::Fp::one() + &P::NONRESIDUE),
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
        cs: CS,
        power: usize,
    ) -> Result<&mut Self, SynthesisError> {
        self.c1
            .mul_by_constant_in_place(cs, &P::FROBENIUS_COEFF_FP2_C1[power % 2])?;
        Ok(self)
    }

    #[inline]
    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Fp2<P>,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        let _ = result.add_constant_in_place(cs, other)?;
        Ok(result)
    }

    #[inline]
    fn add_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        other: &Fp2<P>,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.add_constant_in_place(cs.ns(|| "c0"), &other.c0)?;
        self.c1.add_constant_in_place(cs.ns(|| "c1"), &other.c1)?;
        Ok(self)
    }

    fn mul_by_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        fe: &Fp2<P>,
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

    fn cost_of_mul() -> usize {
        3 * FpGadget::<ConstraintF>::cost_of_mul()
    }

    fn cost_of_mul_equals() -> usize {
        Self::cost_of_mul()
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> PartialEq
    for Fp2Gadget<P, ConstraintF>
{
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> Eq for Fp2Gadget<P, ConstraintF> {}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> EqGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> ConditionalEqGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
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
        Ok(())
    }

    fn cost() -> usize {
        2
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> NEqGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    #[inline]
    fn enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        self.c0.enforce_not_equal(&mut cs.ns(|| "c0"), &other.c0)?;
        self.c1.enforce_not_equal(&mut cs.ns(|| "c1"), &other.c1)?;
        Ok(())
    }

    fn cost() -> usize {
        2
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> ToBitsGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_bits(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_bits(cs.ns(|| "c1"))?;
        c0.append(&mut c1);
        Ok(c0)
    }

    fn to_non_unique_bits<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let mut c0 = self.c0.to_non_unique_bits(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_non_unique_bits(cs.ns(|| "c1"))?;
        c0.append(&mut c1);
        Ok(c0)
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> ToBytesGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
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

    fn to_non_unique_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c0 = self.c0.to_non_unique_bytes(cs.ns(|| "c0"))?;
        let mut c1 = self.c1.to_non_unique_bytes(cs.ns(|| "c1"))?;
        c0.append(&mut c1);
        Ok(c0)
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> Clone
    for Fp2Gadget<P, ConstraintF>
{
    fn clone(&self) -> Self {
        Self {
            c0: self.c0.clone(),
            c1: self.c1.clone(),
            _params: PhantomData,
        }
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> CondSelectGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        true_value: &Self,
        false_value: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = FpGadget::<ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c0"),
            cond,
            &true_value.c0,
            &false_value.c0,
        )?;
        let c1 = FpGadget::<ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c1"),
            cond,
            &true_value.c1,
            &false_value.c1,
        )?;

        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> TwoBitLookupGadget<ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    type TableConstant = Fp2<P>;
    fn two_bit_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 = FpGadget::two_bit_lookup(cs.ns(|| "Lookup c0"), b, &c0s)?;
        let c1 = FpGadget::two_bit_lookup(cs.ns(|| "Lookup c1"), b, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <FpGadget<ConstraintF> as TwoBitLookupGadget<ConstraintF>>::cost()
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField>
    ThreeBitCondNegLookupGadget<ConstraintF> for Fp2Gadget<P, ConstraintF>
{
    type TableConstant = Fp2<P>;

    fn three_bit_cond_neg_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        b0b1: &Boolean,
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 = FpGadget::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c0"), b, b0b1, &c0s)?;
        let c1 = FpGadget::three_bit_cond_neg_lookup(cs.ns(|| "Lookup c1"), b, b0b1, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <FpGadget<ConstraintF> as ThreeBitCondNegLookupGadget<ConstraintF>>::cost()
    }
}

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField> AllocGadget<Fp2<P>, ConstraintF>
    for Fp2Gadget<P, ConstraintF>
{
    #[inline]
    fn alloc_constant<T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<Fp2<P>>,
    {
        Self::zero(cs.ns(|| "zero"))?.add_constant(cs.ns(|| "add constant"), t.borrow())
    }

    #[inline]
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Fp2<P>>,
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

        let c0 = FpGadget::alloc(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = FpGadget::alloc(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Fp2<P>>,
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

        let c0 = FpGadget::alloc_input(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = FpGadget::alloc_input(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }
}
