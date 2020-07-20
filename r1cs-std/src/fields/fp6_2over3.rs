use algebra::{
    fields::{
        fp6_2over3::{Fp6, Fp6Parameters},
        Fp3, Fp3Parameters,
    },
    BigInteger, PrimeField, SquareRootField,
};
use core::{borrow::Borrow, marker::PhantomData};
use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{prelude::*, Vec};

type Fp3Gadget<P, ConstraintF> =
    super::fp3::Fp3Gadget<<P as Fp6Parameters>::Fp3Params, ConstraintF>;
type Fp3GadgetVariable<P, ConstraintF> = <Fp3Gadget<P, ConstraintF> as FieldGadget<
    Fp3<<P as Fp6Parameters>::Fp3Params>,
    ConstraintF,
>>::Variable;

#[derive(Derivative)]
#[derivative(Debug(bound = "ConstraintF: PrimeField + SquareRootField"))]
#[must_use]
pub struct Fp6Gadget<P, ConstraintF: PrimeField + SquareRootField>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    pub c0: Fp3Gadget<P, ConstraintF>,
    pub c1: Fp3Gadget<P, ConstraintF>,
    #[derivative(Debug = "ignore")]
    _params: PhantomData<P>,
}

impl<P, ConstraintF: PrimeField + SquareRootField> ToConstraintFieldGadget<ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    fn to_constraint_field<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<FpGadget<ConstraintF>>, SynthesisError> {
        let mut res = Vec::new();

        let mut c0_gadget = self.c0.to_constraint_field(&mut cs.ns(|| "c0"))?;
        let mut c1_gadget = self.c1.to_constraint_field(&mut cs.ns(|| "c1"))?;

        res.append(&mut c0_gadget);
        res.append(&mut c1_gadget);

        Ok(res)
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    pub fn new(c0: Fp3Gadget<P, ConstraintF>, c1: Fp3Gadget<P, ConstraintF>) -> Self {
        Self {
            c0,
            c1,
            _params: PhantomData,
        }
    }

    /// Multiply a Fp3Gadget by quadratic nonresidue P::NONRESIDUE.
    #[inline]
    pub fn mul_fp3_gadget_by_nonresidue<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        fe: &Fp3Gadget<P, ConstraintF>,
    ) -> Result<Fp3Gadget<P, ConstraintF>, SynthesisError> {
        let mut res = Fp3Gadget::<P, ConstraintF>::new(fe.c2.clone(), fe.c0.clone(), fe.c1.clone());
        res.c0.mul_by_constant_in_place(
            cs.ns(|| "res * non_residue"),
            &<P::Fp3Params as Fp3Parameters>::NONRESIDUE,
        )?;
        Ok(res)
    }

    pub fn unitary_inverse<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Self, SynthesisError> {
        Ok(Self::new(self.c0.clone(), self.c1.negate(cs)?))
    }

    #[inline]
    pub fn cyclotomic_exp<CS: ConstraintSystem<ConstraintF>, B: BigInteger>(
        &self,
        mut cs: CS,
        exponent: &B,
    ) -> Result<Self, SynthesisError> {
        let mut res = Self::one(cs.ns(|| "one"))?;
        let self_inverse = self.unitary_inverse(cs.ns(|| "unitary inverse"))?;

        let mut found_nonzero = false;
        let naf = exponent.find_wnaf();

        for (i, &value) in naf.iter().rev().enumerate() {
            if found_nonzero {
                res.square_in_place(cs.ns(|| format!("square {}", i)))?;
            }

            if value != 0 {
                found_nonzero = true;

                if value > 0 {
                    res.mul_in_place(cs.ns(|| format!("res *= self {}", i)), &self)?;
                } else {
                    res.mul_in_place(
                        cs.ns(|| format!("res *= self_inverse {}", i)),
                        &self_inverse,
                    )?;
                }
            }
        }

        Ok(res)
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> FieldGadget<Fp6<P>, ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    type Variable = (
        Fp3GadgetVariable<P, ConstraintF>,
        Fp3GadgetVariable<P, ConstraintF>,
    );

    #[inline]
    fn get_value(&self) -> Option<Fp6<P>> {
        match (self.c0.get_value(), self.c1.get_value()) {
            (Some(c0), Some(c1)) => Some(Fp6::new(c0, c1)),
            (..) => None,
        }
    }

    #[inline]
    fn get_variable(&self) -> Self::Variable {
        (self.c0.get_variable(), self.c1.get_variable())
    }

    #[inline]
    fn zero<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = Fp3Gadget::<P, ConstraintF>::zero(cs.ns(|| "c0"))?;
        let c1 = Fp3Gadget::<P, ConstraintF>::zero(cs.ns(|| "c1"))?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn one<CS: ConstraintSystem<ConstraintF>>(mut cs: CS) -> Result<Self, SynthesisError> {
        let c0 = Fp3Gadget::<P, ConstraintF>::one(cs.ns(|| "c0"))?;
        let c1 = Fp3Gadget::<P, ConstraintF>::zero(cs.ns(|| "c1"))?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn conditionally_add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        bit: &Boolean,
        coeff: Fp6<P>,
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
        // Karatsuba multiplication for Fp6:
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
                Self::mul_fp3_gadget_by_nonresidue(mul_cs.ns(|| "first mul_by_nr"), &v1)?;
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
        // From Libsnark/fp4_gadget.tcc
        // Complex multiplication for Fp6:
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

        let non_residue_c1 =
            Self::mul_fp3_gadget_by_nonresidue(cs.ns(|| "non_residue * a1"), &self.c1)?;
        let a0_plus_non_residue_c1 = self
            .c0
            .add(cs.ns(|| "a0 + non_residue * a1"), &non_residue_c1)?;
        let one_plus_non_residue_v0 =
            Self::mul_fp3_gadget_by_nonresidue(cs.ns(|| "non_residue * v0"), &v0)?
                .add(cs.ns(|| "plus v0"), &v0)?;

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

    fn mul_equals<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        // Karatsuba multiplication for Fp6:
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
            Self::mul_fp3_gadget_by_nonresidue(mul_cs.ns(|| "nr * v1"), &v1)?;
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
        mut cs: CS,
        power: usize,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0
            .frobenius_map_in_place(cs.ns(|| "frob_map1"), power)?;
        self.c1
            .frobenius_map_in_place(cs.ns(|| "frob_map2"), power)?;
        self.c1
            .mul_by_fp_constant_in_place(cs.ns(|| "mul"), &P::FROBENIUS_COEFF_FP6_C1[power % 6])?;
        Ok(self)
    }

    #[inline]
    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Fp6<P>,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        let _ = result.add_constant_in_place(cs, other)?;
        Ok(result)
    }

    #[inline]
    fn add_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        other: &Fp6<P>,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.add_constant_in_place(cs.ns(|| "c0"), &other.c0)?;
        self.c1.add_constant_in_place(cs.ns(|| "c1"), &other.c1)?;
        Ok(self)
    }

    fn mul_by_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        fe: &Fp6<P>,
    ) -> Result<Self, SynthesisError> {
        // Karatsuba multiplication (see mul above).
        // Doesn't need any constraints; returns linear combinations of
        // `self`'s variables.
        //
        // (The operations below are guaranteed to return linear combinations)
        let (a0, a1) = (&self.c0, &self.c1);
        let (b0, b1) = (fe.c0, fe.c1);
        let mut v0 = a0.mul_by_constant(&mut cs.ns(|| "v0"), &b0)?;
        let mut v1 = Self::mul_fp3_gadget_by_nonresidue(&mut cs.ns(|| "v1"), a1)?;
        let beta_v1 = v1.mul_by_constant_in_place(&mut cs.ns(|| "beta * v1"), &b1)?;

        v0.add_in_place(&mut cs.ns(|| "c0"), &beta_v1)?;
        let c0 = v0;

        let mut a0b1 = a0.mul_by_constant(&mut cs.ns(|| "a0b1"), &b1)?;
        let a1b0 = a1.mul_by_constant(&mut cs.ns(|| "a1b0"), &b0)?;
        a0b1.add_in_place(&mut cs.ns(|| "c1"), &a1b0)?;
        let c1 = a0b1;
        Ok(Self::new(c0, c1))
    }

    fn cost_of_mul() -> usize {
        2 * Fp3Gadget::<P, ConstraintF>::cost_of_mul()
    }

    fn cost_of_mul_equals() -> usize {
        Self::cost_of_mul()
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> PartialEq for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    fn eq(&self, other: &Self) -> bool {
        self.c0 == other.c0 && self.c1 == other.c1
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> Eq for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
}

impl<P, ConstraintF: PrimeField + SquareRootField> EqGadget<ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
}

impl<P, ConstraintF: PrimeField + SquareRootField> ConditionalEqGadget<ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
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
        2 * <Fp3Gadget<P, ConstraintF> as ConditionalEqGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> NEqGadget<ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
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
        2 * <Fp3Gadget<P, ConstraintF> as NEqGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> ToBitsGadget<ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
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

impl<P, ConstraintF: PrimeField + SquareRootField> ToBytesGadget<ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
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

impl<P, ConstraintF: PrimeField + SquareRootField> Clone for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    fn clone(&self) -> Self {
        Self {
            c0: self.c0.clone(),
            c1: self.c1.clone(),
            _params: PhantomData,
        }
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> CondSelectGadget<ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError> {
        let c0 = Fp3Gadget::<P, ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c0"),
            cond,
            &first.c0,
            &second.c0,
        )?;
        let c1 = Fp3Gadget::<P, ConstraintF>::conditionally_select(
            &mut cs.ns(|| "c1"),
            cond,
            &first.c1,
            &second.c1,
        )?;

        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <Fp3Gadget<P, ConstraintF> as CondSelectGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> TwoBitLookupGadget<ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    type TableConstant = Fp6<P>;
    fn two_bit_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 = Fp3Gadget::<P, ConstraintF>::two_bit_lookup(cs.ns(|| "Lookup c0"), b, &c0s)?;
        let c1 = Fp3Gadget::<P, ConstraintF>::two_bit_lookup(cs.ns(|| "Lookup c1"), b, &c1s)?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <Fp3Gadget<P, ConstraintF> as TwoBitLookupGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> ThreeBitCondNegLookupGadget<ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    type TableConstant = Fp6<P>;

    fn three_bit_cond_neg_lookup<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        b: &[Boolean],
        b0b1: &Boolean,
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        let c0s = c.iter().map(|f| f.c0).collect::<Vec<_>>();
        let c1s = c.iter().map(|f| f.c1).collect::<Vec<_>>();
        let c0 = Fp3Gadget::<P, ConstraintF>::three_bit_cond_neg_lookup(
            cs.ns(|| "Lookup c0"),
            b,
            b0b1,
            &c0s,
        )?;
        let c1 = Fp3Gadget::<P, ConstraintF>::three_bit_cond_neg_lookup(
            cs.ns(|| "Lookup c1"),
            b,
            b0b1,
            &c1s,
        )?;
        Ok(Self::new(c0, c1))
    }

    fn cost() -> usize {
        2 * <Fp3Gadget<P, ConstraintF> as ThreeBitCondNegLookupGadget<ConstraintF>>::cost()
    }
}

impl<P, ConstraintF: PrimeField + SquareRootField> AllocGadget<Fp6<P>, ConstraintF>
    for Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    #[inline]
    fn alloc_constant<T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<Fp6<P>>,
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
        T: Borrow<Fp6<P>>,
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

        let c0 = Fp3Gadget::<P, ConstraintF>::alloc(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = Fp3Gadget::<P, ConstraintF>::alloc(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }

    #[inline]
    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Fp6<P>>,
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

        let c0 = Fp3Gadget::<P, ConstraintF>::alloc_input(&mut cs.ns(|| "c0"), || c0)?;
        let c1 = Fp3Gadget::<P, ConstraintF>::alloc_input(&mut cs.ns(|| "c1"), || c1)?;
        Ok(Self::new(c0, c1))
    }
}
