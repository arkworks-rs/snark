use crate::{
    gr1cs::{predicate::PredicateConstraintSystem, ConstraintSynthesizer, ConstraintSystemRef},
    lc, ns,
};
use ark_ff::Field;
use ark_std::{collections::BTreeMap, string::ToString, vec::Vec};

use super::{Label, Matrix, R1CS_PREDICATE_LABEL};

#[derive(Debug, Clone)]
pub struct Circuit1<F: Field> {
    pub x1: F,
    pub x2: F,
    pub x3: F,
    pub x4: F,
    pub x5: F,
    pub w1: F,
    pub w2: F,
    pub w3: F,
    pub w4: F,
    pub w5: F,
    pub w6: F,
    pub w7: F,
    pub w8: F,
}

impl<F: Field> Circuit1<F> {
    pub fn get_matrices() -> BTreeMap<Label, Vec<Matrix<F>>> {
        let mut map: BTreeMap<Label, Vec<Matrix<F>>> = BTreeMap::new();
        map.insert(
            R1CS_PREDICATE_LABEL.to_string(),
            vec![vec![], vec![], vec![]],
        );
        map.insert(
            "poly-predicate-A".to_string(),
            vec![
                vec![vec![(F::one(), 1)]],
                vec![vec![(F::one(), 2)]],
                vec![vec![(F::one(), 3)]],
                vec![vec![(F::one(), 9)]],
            ],
        );

        map.insert(
            "poly-predicate-B".to_string(),
            vec![
                vec![vec![(F::one(), 4)], vec![(F::one(), 10)]],
                vec![vec![(F::one(), 6)], vec![(F::one(), 11)]],
                vec![vec![(F::one(), 10)], vec![(F::one(), 13)]],
            ],
        );
        map.insert(
            "poly-predicate-C".to_string(),
            vec![
                vec![vec![(F::one(), 7)], vec![(F::one(), 9), (F::one(), 10)]],
                vec![vec![(F::one(), 8)], vec![(F::one(), 13)]],
                vec![vec![(F::one(), 11)], vec![(F::one(), 5)]],
            ],
        );
        map
    }
}
impl<F: Field + core::convert::From<i8>> ConstraintSynthesizer<F> for Circuit1<F> {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> crate::utils::Result<()> {
        // Variable declarations -> Instance variables + Witness variables
        let input_variables_namespace = ns!(cs, "Input variables");
        let x1 = input_variables_namespace
            .cs()
            .new_input_variable(|| Ok(self.x1))
            .unwrap();
        let x2 = input_variables_namespace
            .cs()
            .new_input_variable(|| Ok(self.x2))
            .unwrap();
        let x3 = input_variables_namespace
            .cs()
            .new_input_variable(|| Ok(self.x3))
            .unwrap();
        let x4 = input_variables_namespace
            .cs()
            .new_input_variable(|| Ok(self.x4))
            .unwrap();
        let x5 = input_variables_namespace
            .cs()
            .new_input_variable(|| Ok(self.x5))
            .unwrap();
        ns!(cs, "Witness variables");
        let w1 = cs.new_witness_variable(|| Ok(self.w1)).unwrap();
        let w2 = cs.new_witness_variable(|| Ok(self.w2)).unwrap();
        let w3 = cs.new_witness_variable(|| Ok(self.w3)).unwrap();
        let w4 = cs.new_witness_variable(|| Ok(self.w4)).unwrap();
        let w5 = cs.new_witness_variable(|| Ok(self.w5)).unwrap();
        let w6 = cs.new_witness_variable(|| Ok(self.w6)).unwrap();
        let _w7 = cs.new_witness_variable(|| Ok(self.w7)).unwrap();
        let w8 = cs.new_witness_variable(|| Ok(self.w8)).unwrap();

        //  predicate declarations -> Polynomial predicates

        let one = F::ONE;
        let three = one + one + one;
        let seven = three + three + one;
        let predicate_a = PredicateConstraintSystem::new_polynomial_predicate_cs(
            4,
            vec![
                (one, vec![(0, 1), (1, 1)]),
                (three, vec![(2, 2)]),
                (-one, vec![(3, 1)]),
            ],
        );
        let predicate_b = PredicateConstraintSystem::new_polynomial_predicate_cs(
            3,
            vec![
                (seven, vec![(1, 1)]),
                (one, vec![(0, 3)]),
                (-one, vec![(2, 1)]),
            ],
        );

        let predicate_c = PredicateConstraintSystem::new_polynomial_predicate_cs(
            3,
            vec![(one, vec![(0, 1), (1, 1)]), (-one, vec![(2, 1)])],
        );
        cs.register_predicate("poly-predicate-A", predicate_a)?;
        cs.register_predicate("poly-predicate-B", predicate_b)?;
        cs.register_predicate("poly-predicate-C", predicate_c)?;

        // Enforing constraints to the predicates

        ns!(cs, "Predicate A constraints");
        cs.enforce_constraint_arity_4(
            "poly-predicate-A",
            || lc!() + x1,
            || lc!() + x2,
            || lc!() + x3,
            || lc!() + w4,
        )?;
        ns!(cs, "Predicate B constraints");
        cs.enforce_constraint_arity_3(
            "poly-predicate-B",
            || lc!() + x4,
            || lc!() + w1,
            || lc!() + w5,
        )?;
        cs.enforce_constraint_arity_3(
            "poly-predicate-B",
            || lc!() + w5,
            || lc!() + w6,
            || lc!() + w8,
        )?;
        ns!(cs, "Predicate C constraints");
        cs.enforce_constraint_arity_3(
            "poly-predicate-C",
            || lc!() + w2,
            || lc!() + w3,
            || lc!() + w6,
        )?;
        cs.enforce_constraint_arity_3(
            "poly-predicate-C",
            || lc!() + w5 + w4,
            || lc!() + w8,
            || lc!() + x5,
        )?;
        Ok(())
    }
}
