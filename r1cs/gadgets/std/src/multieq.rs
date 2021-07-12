use algebra::{PrimeField, FpParameters};

use r1cs_core::{ConstraintSystem, LinearCombination, SynthesisError, Variable};

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
        unimplemented!()
    }
}