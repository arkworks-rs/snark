use ark_ff::{Field, One, Zero};
use ark_std::{cfg_iter, vec::Vec};

use crate::{
    gr1cs::{
        predicate::{self, polynomial_constraint::SR1CS_PREDICATE_LABEL},
        ConstraintSystem, ConstraintSystemRef, LinearCombination, OptimizationGoal, SynthesisError,
        SynthesisMode, Variable, R1CS_PREDICATE_LABEL,
    },
    lc, lc_diff,
};
use ark_std::{collections::BTreeMap, ops::AddAssign};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Helper struct to enable easy conversion from R1CS to SR1CS
pub struct Sr1csAdapter<F: Field> {
    _phantom: core::marker::PhantomData<F>,
}

impl<F: Field> Sr1csAdapter<F> {
    /// Computes the inner product of `terms` with `assignment`.
    pub fn evaluate_constraint<'a, LHS, RHS, R>(
        terms: &'a [(LHS, usize)],
        assignment: &'a [RHS],
    ) -> R
    where
        LHS: One + Send + Sync + PartialEq,
        RHS: Send + Sync + core::ops::Mul<&'a LHS, Output = RHS> + Copy,
        R: Zero + Send + Sync + AddAssign<RHS> + core::iter::Sum,
    {
        // Need to wrap in a closure when using Rayon
        #[cfg(feature = "parallel")]
        let zero = || R::zero();
        #[cfg(not(feature = "parallel"))]
        let zero = R::zero();

        let res = cfg_iter!(terms).fold(zero, |mut sum, (coeff, index)| {
            let val = &assignment[*index];

            if coeff.is_one() {
                sum += *val;
            } else {
                sum += val.mul(coeff);
            }

            sum
        });

        // Need to explicitly call `.sum()` when using Rayon
        #[cfg(feature = "parallel")]
        return res.sum();
        #[cfg(not(feature = "parallel"))]
        return res;
    }

    fn add_to_variable_maps(
        constraint: &[(F, usize)],
        public_variables: &mut BTreeMap<usize, Variable>,
        witness_variables: &mut BTreeMap<usize, Variable>,
        num_public: usize,
        new_cs: &ConstraintSystemRef<F>,
    ) -> LinearCombination<F> {
        let lc = constraint
            .iter()
            .map(|(coeff, index)| {
                let var = if *index == 0 {
                    Variable::One
                } else if *index < num_public {
                    *public_variables
                        .entry(*index)
                        .or_insert_with(|| new_cs.new_witness_variable(|| Ok(F::ONE)).unwrap())
                } else {
                    *witness_variables
                        .entry(*index)
                        .or_insert_with(|| new_cs.new_witness_variable(|| Ok(F::ONE)).unwrap())
                };
                (*coeff, var)
            })
            .collect();
        LinearCombination(lc)
    }

    fn add_to_variable_maps_witness(
        constraint: &[(F, usize)],
        public_variables: &mut BTreeMap<usize, Variable>,
        witness_variables: &mut BTreeMap<usize, Variable>,
        num_public: usize,
        r1cs_assignment: &[F],
        new_cs: &ConstraintSystemRef<F>,
    ) -> (LinearCombination<F>, F) {
        let (lc, vals): (Vec<_>, Vec<_>) = constraint
            .iter()
            .map(|(coeff, index)| {
                let (var, val) = if *index == 0 {
                    (Variable::One, F::ONE)
                } else if *index < num_public {
                    let val = r1cs_assignment[*index];
                    let var = *public_variables
                        .entry(*index)
                        .or_insert_with(|| new_cs.new_witness_variable(|| Ok(val)).unwrap());
                    (var, val)
                } else {
                    let val = r1cs_assignment[*index];
                    let var = *witness_variables
                        .entry(*index)
                        .or_insert_with(|| new_cs.new_witness_variable(|| Ok(val)).unwrap());
                    (var, val)
                };
                ((*coeff, var), *coeff * val)
            })
            .unzip();
        let val = vals.into_iter().sum::<F>();
        (LinearCombination(lc), val)
    }

    /// Converts an R1CS constraint system to an SR1CS constraint system.
    pub fn r1cs_to_sr1cs(
        cs: &ConstraintSystemRef<F>,
    ) -> Result<ConstraintSystemRef<F>, SynthesisError> {
        assert_eq!(cs.num_predicates(), 1);
        let matrices = &cs.to_matrices().unwrap()[R1CS_PREDICATE_LABEL];
        let mut public_variables = BTreeMap::new();
        let mut witness_variables = BTreeMap::new();
        let num_public = cs.num_instance_variables();

        let new_cs = ConstraintSystem::new_ref();
        // new_cs.set_optimization_goal(OptimizationGoal::Constraints);
        new_cs.set_mode(SynthesisMode::Setup);
        new_cs.remove_predicate(R1CS_PREDICATE_LABEL);
        let _ = new_cs.register_predicate(
            SR1CS_PREDICATE_LABEL,
            predicate::PredicateConstraintSystem::new_sr1cs_predicate().unwrap(),
        );
        for ((a_i, b_i), c_i) in matrices[0].iter().zip(&matrices[1]).zip(&matrices[2]) {
            let a_i = Self::add_to_variable_maps(
                a_i,
                &mut public_variables,
                &mut witness_variables,
                num_public,
                &new_cs,
            );
            let b_i = Self::add_to_variable_maps(
                b_i,
                &mut public_variables,
                &mut witness_variables,
                num_public,
                &new_cs,
            );
            let mut c_i = Self::add_to_variable_maps(
                c_i,
                &mut public_variables,
                &mut witness_variables,
                num_public,
                &new_cs,
            );
            let square_variable = new_cs.new_witness_variable(|| Ok(F::ONE))?;

            let left_1 = || a_i.clone() + &b_i;
            c_i.0.iter_mut().for_each(|(c, _)| {
                c.double_in_place().double_in_place();
            });
            let right_1 = || c_i + square_variable;
            new_cs.enforce_sr1cs_constraint(left_1, right_1)?;

            let left_2 = || a_i - b_i;
            let right_2 = || lc![square_variable];
            new_cs.enforce_sr1cs_constraint(left_2, right_2)?;
        }

        for old_var in public_variables.values() {
            let new_var = new_cs.new_input_variable(|| Ok(F::ONE))?;
            let lc = || lc_diff![*old_var, new_var];
            new_cs.enforce_sr1cs_constraint(lc, || lc![])?;
        }
        Ok(new_cs)
    }

    /// Converts an R1CS constraint system to an SR1CS constraint system,
    /// while also converting the R1CS assignment to an equivalent SR1CS assignment.
    pub fn r1cs_to_sr1cs_with_assignment(
        cs: &mut ConstraintSystem<F>,
    ) -> Result<ConstraintSystemRef<F>, SynthesisError> {
        let matrices = &cs.to_matrices().unwrap()[R1CS_PREDICATE_LABEL];
        let mut public_variables = BTreeMap::new();
        let mut witness_variables = BTreeMap::new();
        let num_public = cs.num_instance_variables();

        let mut r1cs_assignment = cs.assignments.instance_assignment.clone();
        r1cs_assignment.extend_from_slice(&cs.assignments.witness_assignment);

        let new_cs = ConstraintSystem::new_ref();
        new_cs.remove_predicate(R1CS_PREDICATE_LABEL);
        let _ = new_cs.register_predicate(
            SR1CS_PREDICATE_LABEL,
            predicate::PredicateConstraintSystem::new_sr1cs_predicate().unwrap(),
        );
        new_cs.set_optimization_goal(OptimizationGoal::Constraints);
        cs.set_mode(SynthesisMode::Prove {
            construct_matrices: true,
            generate_lc_assignments: true,
        });
        for ((a_i, b_i), c_i) in matrices[0].iter().zip(&matrices[1]).zip(&matrices[2]) {
            let (a_i, a_val) = Self::add_to_variable_maps_witness(
                a_i,
                &mut public_variables,
                &mut witness_variables,
                num_public,
                &r1cs_assignment,
                &new_cs,
            );
            let (b_i, b_val) = Self::add_to_variable_maps_witness(
                b_i,
                &mut public_variables,
                &mut witness_variables,
                num_public,
                &r1cs_assignment,
                &new_cs,
            );
            let (mut c_i, _c_val) = Self::add_to_variable_maps_witness(
                c_i,
                &mut public_variables,
                &mut witness_variables,
                num_public,
                &r1cs_assignment,
                &new_cs,
            );
            let square_variable = new_cs.new_witness_variable(|| Ok((a_val - b_val).square()))?;

            let left_1 = || a_i.clone() + &b_i;
            c_i.0.iter_mut().for_each(|(c, _)| {
                c.double_in_place().double_in_place();
            });
            let right_1 = || c_i + square_variable;
            new_cs.enforce_sr1cs_constraint(left_1, right_1)?;

            let left_2 = || a_i - &b_i;
            let right_2 = || lc![square_variable];
            new_cs.enforce_sr1cs_constraint(left_2, right_2)?;
        }

        let mut public_variable_polys = Vec::new();
        for old_var in public_variables.values() {
            let value = new_cs
                .assigned_value(*old_var)
                .ok_or(SynthesisError::AssignmentMissing)?;
            let constraint_number = new_cs.num_constraints();
            let new_var = new_cs.new_input_variable(|| Ok(value))?;
            let lc = || lc_diff![*old_var, new_var];
            new_cs.enforce_sr1cs_constraint(lc, || lc![])?;
            public_variable_polys.push(constraint_number);
        }
        new_cs.finalize();
        Ok(new_cs)
    }
}

#[cfg(test)]
mod tests {
    use crate::gr1cs::ConstraintSynthesizer;
    use ark_ff::PrimeField;
    use ark_test_curves::bls12_381::Fr;

    use super::*;
    #[derive(Copy)]
    struct DummyCircuit<F: PrimeField> {
        pub a: Option<F>,
        pub b: Option<F>,
        pub num_variables: usize,
        pub num_constraints: usize,
    }

    impl<F: PrimeField> Clone for DummyCircuit<F> {
        fn clone(&self) -> Self {
            DummyCircuit {
                a: self.a.clone(),
                b: self.b.clone(),
                num_variables: self.num_variables.clone(),
                num_constraints: self.num_constraints.clone(),
            }
        }
    }

    impl<F: PrimeField> ConstraintSynthesizer<F> for DummyCircuit<F> {
        fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
            let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
            let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
            let c = cs.new_input_variable(|| {
                let a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
                let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;

                Ok(a * b)
            })?;

            for _ in 0..(self.num_variables - 3) {
                let _ =
                    cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
            }

            for _ in 0..self.num_constraints - 1 {
                cs.enforce_r1cs_constraint(|| lc![a], || lc![b], || lc![c])?;
            }

            cs.enforce_r1cs_constraint(|| lc![], || lc![], || lc![])?;

            Ok(())
        }
    }
    #[test]
    fn r1cs_to_sr1cs() {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let circuit = DummyCircuit {
            a: Some(3u64.into()),
            b: Some(5u64.into()),
            num_variables: 128,
            num_constraints: 128,
        };
        circuit.generate_constraints(cs.clone()).unwrap();
    }
}
