use algebra::Field;
use r1cs_core::{
    lc, ConstraintSynthesizer, ConstraintSystemRef, LinearCombination, SynthesisError, Variable,
};
use std::marker::PhantomData;

pub struct Benchmark<F: Field> {
    num_constraints: usize,
    _engine: PhantomData<F>,
}

impl<F: Field> Benchmark<F> {
    pub fn new(num_constraints: usize) -> Self {
        Self {
            num_constraints,
            _engine: PhantomData,
        }
    }
}

impl<F: Field> ConstraintSynthesizer<F> for Benchmark<F> {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
        let mut assignments = Vec::new();
        let mut a_val = F::one();
        let mut a_var = cs.new_input_variable(|| Ok(a_val))?;
        assignments.push((a_val, a_var));

        let mut b_val = F::one();
        let mut b_var = cs.new_input_variable(|| Ok(b_val))?;
        assignments.push((a_val, a_var));

        for i in 0..self.num_constraints - 1 {
            if i % 2 != 0 {
                let c_val = a_val * &b_val;
                let c_var = cs.new_witness_variable(|| Ok(c_val))?;

                cs.enforce_constraint(lc!() + a_var, lc!() + b_var, lc!() + c_var)?;

                assignments.push((c_val, c_var));
                a_val = b_val;
                a_var = b_var;
                b_val = c_val;
                b_var = c_var;
            } else {
                let c_val = a_val + &b_val;
                let c_var = cs.new_witness_variable(|| Ok(c_val))?;

                cs.enforce_constraint(lc!() + a_var + b_var, lc!() + Variable::One, lc!() + c_var)?;

                assignments.push((c_val, c_var));
                a_val = b_val;
                a_var = b_var;
                b_val = c_val;
                b_var = c_var;
            }
        }

        let mut a_lc = LinearCombination::zero();
        let mut b_lc = LinearCombination::zero();
        let mut c_val = F::zero();

        for (val, var) in assignments {
            a_lc = a_lc + var;
            b_lc = b_lc + var;
            c_val = c_val + &val;
        }
        c_val = c_val.square();

        let c_var = cs.new_witness_variable(|| Ok(c_val))?;

        cs.enforce_constraint(lc!() + a_lc, lc!() + b_lc, lc!() + c_var)?;

        Ok(())
    }
}
