use algebra::PrimeField;
use primitives::SBox;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{
    alloc::AllocGadget,
    bits::boolean::Boolean,
    eq::EqGadget,
    fields::{fp::FpGadget, FieldGadget},
    Assignment,
};
use std::marker::PhantomData;

pub trait SBoxGadget<ConstraintF: PrimeField, SB: SBox<Field = ConstraintF>> {
    /// Enforce S(x)
    fn apply<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        x: &mut FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError>;
}

pub struct InverseSBoxGadget<ConstraintF: PrimeField, SB: SBox<Field = ConstraintF>> {
    _field: PhantomData<ConstraintF>,
    _sbox: PhantomData<SB>,
}

impl<ConstraintF: PrimeField, SB: SBox<Field = ConstraintF>> SBoxGadget<ConstraintF, SB>
    for InverseSBoxGadget<ConstraintF, SB>
{
    // Enforce S(x) = X^-1 if X != 0 otherwise X
    fn apply<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        x: &mut FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let b = Boolean::alloc(cs.ns(|| "alloc b"), || {
            let x_val = x.get_value().get()?;
            if x_val == ConstraintF::zero() {
                Ok(false)
            } else {
                Ok(true)
            }
        })?;

        let y = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc y"), || {
            let x_val = x.get_value().get()?;
            if x_val == ConstraintF::zero() {
                Ok(x_val)
            } else {
                let inv = x_val.inverse().get()?;
                Ok(inv)
            }
        })?;

        cs.enforce(
            || "b = x * y",
            |lc| &x.variable + lc,
            |lc| &y.variable + lc,
            |lc| lc + &b.lc(CS::one(), ConstraintF::one()),
        );
        x.conditional_enforce_equal(cs.ns(|| "0 = (1 - b) * (x - y)"), &y, &b.not())?;

        *x = y;
        Ok(())
    }
}

pub struct QuinticSBoxGadget<ConstraintF: PrimeField, SB: SBox<Field = ConstraintF>> {
    _field: PhantomData<ConstraintF>,
    _sbox: PhantomData<SB>,
}

impl<ConstraintF: PrimeField, SB: SBox<Field = ConstraintF>> SBoxGadget<ConstraintF, SB>
    for QuinticSBoxGadget<ConstraintF, SB>
{
    // Enforce S(X) = X^5
    fn apply<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        x: &mut FpGadget<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        let x4 = x.square(cs.ns(|| "x^2"))?.square(cs.ns(|| "x^4"))?;
        x.mul_in_place(cs.ns(|| "x^5"), &x4)?;
        Ok(())
    }
}
