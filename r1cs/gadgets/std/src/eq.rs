use crate::prelude::*;
use algebra::{Field, PrimeField, FpParameters};
use r1cs_core::{ConstraintSystem, LinearCombination, SynthesisError, Variable};

/// Specifies how to generate constraints that check for equality for two variables of type `Self`.
pub trait EqGadget<ConstraintF: Field>: Eq {
    /// Output a `Boolean` value representing whether `self.value() == other.value()`.
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS, other: &Self) -> Result<Boolean, SynthesisError>;

    /// Output a `Boolean` value representing whether `self.value() != other.value()`.
    ///
    /// By default, this is defined as `self.is_eq(other)?.not()`.
    fn is_neq<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS, other: &Self) -> Result<Boolean, SynthesisError> {
        Ok(self.is_eq(cs, other)?.not())
    }

    /// If `should_enforce == true`, enforce that `self` and `other` are equal; else,
    /// enforce a vacuously true statement.
    ///
    /// A safe default implementation is provided that generates the following constraints:
    /// `self.is_eq(other)?.conditional_enforce_equal(&Boolean::TRUE, should_enforce)`.
    ///
    /// More efficient specialized implementation may be possible; implementors
    /// are encouraged to carefully analyze the efficiency and safety of these.
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.is_eq(cs.ns(|| "is_eq(self, other)"), &other)?
            .conditional_enforce_equal(cs.ns(|| "enforce condition"), &Boolean::constant(true), should_enforce)
    }

    /// Enforce that `self` and `other` are equal.
    ///
    /// A safe default implementation is provided that generates the following constraints:
    /// `self.conditional_enforce_equal(other, &Boolean::TRUE)`.
    ///
    /// More efficient specialized implementation may be possible; implementors
    /// are encouraged to carefully analyze the efficiency and safety of these.
    fn enforce_equal<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS, other: &Self) -> Result<(), SynthesisError> {
        self.conditional_enforce_equal(cs, other, &Boolean::constant(true))
    }

    /// If `should_enforce == true`, enforce that `self` and `other` are *not* equal; else,
    /// enforce a vacuously true statement.
    ///
    /// A safe default implementation is provided that generates the following constraints:
    /// `self.is_neq(other)?.conditional_enforce_equal(&Boolean::TRUE, should_enforce)`.
    ///
    /// More efficient specialized implementation may be possible; implementors
    /// are encouraged to carefully analyze the efficiency and safety of these.
    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.is_neq(cs.ns(|| "is_neq(self, other)"), &other)?
            .conditional_enforce_equal(cs.ns(|| "enforce condition"), &Boolean::constant(true), should_enforce)
    }

    /// Enforce that `self` and `other` are *not* equal.
    ///
    /// A safe default implementation is provided that generates the following constraints:
    /// `self.conditional_enforce_not_equal(other, &Boolean::TRUE)`.
    ///
    /// More efficient specialized implementation may be possible; implementors
    /// are encouraged to carefully analyze the efficiency and safety of these.
    fn enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS, other: &Self) -> Result<(), SynthesisError> {
        self.conditional_enforce_not_equal(cs, other, &Boolean::constant(true))
    }
}

impl<T: EqGadget<ConstraintF>, ConstraintF: Field> EqGadget<ConstraintF> for [T] {
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(&self, mut cs: CS, other: &Self) -> Result<Boolean, SynthesisError> {
        assert_eq!(self.len(), other.len());
        assert!(!self.is_empty());
        let mut results = Vec::with_capacity(self.len());
        for (i ,(a, b)) in self.iter().zip(other).enumerate() {
            results.push(a.is_eq(cs.ns(|| format!("is_eq_{}", i)), b)?);
        }
        Boolean::kary_and(cs.ns(|| "kary and"), &results)
    }

    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        assert_eq!(self.len(), other.len());
        for (i, (a, b)) in self.iter().zip(other).enumerate() {
            a.conditional_enforce_equal(
                cs.ns(|| format!("conditional_enforce_equal_{}", i)),
                b,
                condition
            )?;
        }
        Ok(())
    }

    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        assert_eq!(self.len(), other.len());
        let some_are_different = self.is_neq(cs.ns(|| "is_neq"), other)?;
        if some_are_different.get_value().is_some() && should_enforce.get_value().is_some() {
            assert!(some_are_different.get_value().unwrap());
            Ok(())
        } else {
            some_are_different.conditional_enforce_equal(
                cs.ns(|| "conditional_enforce_equal"),
                should_enforce,
                should_enforce,
            )?;
            Ok(())
        }
    }
}

pub struct MultiEq<ConstraintF: PrimeField, CS: ConstraintSystem<ConstraintF>> {
    cs: CS,
    ops: usize,
    bits_used: usize,
    lhs: LinearCombination<ConstraintF>,
    rhs: LinearCombination<ConstraintF>,
}

impl<ConstraintF: PrimeField, CS: ConstraintSystem<ConstraintF>> MultiEq<ConstraintF, CS> {
    pub fn new(cs: CS) -> Self {
        MultiEq {
            cs,
            ops: 0,
            bits_used: 0,
            lhs: LinearCombination::zero(),
            rhs: LinearCombination::zero(),
        }
    }

    fn accumulate(&mut self) {
        let ops = self.ops;
        let lhs = self.lhs.clone();
        let rhs = self.rhs.clone();
        self.cs.enforce(
            || format!("multieq {}", ops),
            |_| lhs,
            |lc| lc + CS::one(),
            |_| rhs,
        );
        self.lhs = LinearCombination::zero();
        self.rhs = LinearCombination::zero();
        self.bits_used = 0;
        self.ops += 1;
    }

    pub fn enforce_equal(
        &mut self,
        num_bits: usize,
        lhs: &LinearCombination<ConstraintF>,
        rhs: &LinearCombination<ConstraintF>,
    ) {
        // Check if we will exceed the capacity
        if (ConstraintF::Params::CAPACITY as usize) <= (self.bits_used + num_bits) {
            self.accumulate();
        }

        assert!((ConstraintF::Params::CAPACITY as usize) > (self.bits_used + num_bits));

        let frmstr = ConstraintF::from_str("2").unwrap_or_default();

        let coeff = frmstr.pow(&[self.bits_used as u64]);
        self.lhs = self.lhs.clone() + (coeff, lhs);
        self.rhs = self.rhs.clone() + (coeff, rhs);
        self.bits_used += num_bits;
    }
}

impl<ConstraintF: PrimeField, CS: ConstraintSystem<ConstraintF>> Drop for MultiEq<ConstraintF, CS> {
    fn drop(&mut self) {
        if self.bits_used > 0 {
            self.accumulate();
        }
    }
}

impl<ConstraintF: PrimeField, CS: ConstraintSystem<ConstraintF>> ConstraintSystem<ConstraintF>
for MultiEq<ConstraintF, CS>
{
    type Root = Self;

    fn one() -> Variable {
        CS::one()
    }

    fn alloc<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
        where
            F: FnOnce() -> Result<ConstraintF, SynthesisError>,
            A: FnOnce() -> AR,
            AR: Into<String>,
    {
        self.cs.alloc(annotation, f)
    }

    fn alloc_input<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
        where
            F: FnOnce() -> Result<ConstraintF, SynthesisError>,
            A: FnOnce() -> AR,
            AR: Into<String>,
    {
        self.cs.alloc_input(annotation, f)
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
        where
            A: FnOnce() -> AR,
            AR: Into<String>,
            LA: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
            LB: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
            LC: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
    {
        self.cs.enforce(annotation, a, b, c)
    }

    fn push_namespace<NR, N>(&mut self, name_fn: N)
        where
            NR: Into<String>,
            N: FnOnce() -> NR,
    {
        self.cs.get_root().push_namespace(name_fn)
    }

    fn pop_namespace(&mut self) {
        self.cs.get_root().pop_namespace()
    }

    fn get_root(&mut self) -> &mut Self::Root {
        self
    }

    fn num_constraints(&self) -> usize {
        self.cs.num_constraints()
    }
}