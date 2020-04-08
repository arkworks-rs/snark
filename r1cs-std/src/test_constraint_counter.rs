use crate::String;
use algebra::Field;
use r1cs_core::{ConstraintSystem, Index, LinearCombination, SynthesisError, Variable};

/// Constraint counter for testing purposes.
pub struct ConstraintCounter {
    pub num_inputs:      usize,
    pub num_aux:         usize,
    pub num_constraints: usize,
}

impl ConstraintCounter {
    pub fn new() -> Self {
        Self {
            num_aux:         0,
            num_inputs:      0,
            num_constraints: 0,
        }
    }

    pub fn num_constraints(&self) -> usize {
        self.num_constraints
    }
}

impl<ConstraintF: Field> ConstraintSystem<ConstraintF> for ConstraintCounter {
    type Root = Self;

    fn alloc<F, A, AR>(&mut self, _: A, _: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<ConstraintF, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        let var = Variable::new_unchecked(Index::Aux(self.num_aux));
        self.num_aux += 1;
        Ok(var)
    }

    fn alloc_input<F, A, AR>(&mut self, _: A, _: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<ConstraintF, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        let var = Variable::new_unchecked(Index::Input(self.num_inputs));
        self.num_inputs += 1;

        Ok(var)
    }

    fn enforce<A, AR, LA, LB, LC>(&mut self, _: A, _: LA, _: LB, _: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
        LB: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
        LC: FnOnce(LinearCombination<ConstraintF>) -> LinearCombination<ConstraintF>,
    {
        self.num_constraints += 1;
    }

    fn push_namespace<NR, N>(&mut self, _: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
    }

    fn pop_namespace(&mut self) {}

    fn get_root(&mut self) -> &mut Self::Root {
        self
    }

    fn num_constraints(&self) -> usize {
        self.num_constraints
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_constraint_system::TestConstraintSystem;
    use algebra::{bls12_381::Fq, Field};
    use r1cs_core::{ConstraintSynthesizer, SynthesisError};

    // circuit proving knowledge of a square root
    #[derive(Clone, Debug)]
    struct TestCircuit<F>(Option<F>);

    impl<F: Field> ConstraintSynthesizer<F> for TestCircuit<F> {
        fn generate_constraints<CS: ConstraintSystem<F>>(
            self,
            cs: &mut CS,
        ) -> Result<(), SynthesisError> {
            let x = cs.alloc(|| "x", || self.0.ok_or(SynthesisError::AssignmentMissing))?;
            // 1 input!
            let out = cs.alloc_input(
                || "square",
                || {
                    self.0
                        .map(|x| x.square())
                        .ok_or(SynthesisError::AssignmentMissing)
                },
            )?;
            cs.enforce(|| "x * x = x^2", |lc| lc + x, |lc| lc + x, |lc| lc + out);
            Ok(())
        }
    }

    #[test]
    fn test_constraints_counter() {
        let empty_circuit = TestCircuit::<Fq>(None);
        let populated_circuit = TestCircuit(Some(Fq::from(10u32)));

        let mut counter = ConstraintCounter::new();
        let mut cs = TestConstraintSystem::new();

        empty_circuit
            .clone()
            .generate_constraints(&mut counter)
            .unwrap();
        // an empty circuit cannot be used with TestConstraintSystem
        empty_circuit.generate_constraints(&mut cs).unwrap_err();
        populated_circuit.generate_constraints(&mut cs).unwrap();

        assert_eq!(counter.num_constraints(), cs.num_constraints())
    }
}
