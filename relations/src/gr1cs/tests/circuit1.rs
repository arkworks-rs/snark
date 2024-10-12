use ark_ff::Field;
use ark_std::{
    collections::BTreeMap,
    string::ToString,
    vec::Vec,
};

use crate::{
    gr1cs::{
        local_predicate::PredicateConstraintSystem,
        ConstraintSynthesizer, ConstraintSystemRef,
    },
    lc,
};

use super::{Label, Matrix};

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
    fn generate_constraints(self, mut cs: ConstraintSystemRef<F>) -> crate::utils::Result<()> {
        // Variable declarations -> Instance variables + Witness variables

        let x1 = cs.new_input_variable(|| Ok(self.x1)).unwrap();
        let x2 = cs.new_input_variable(|| Ok(self.x2)).unwrap();
        let x3 = cs.new_input_variable(|| Ok(self.x3)).unwrap();
        let x4 = cs.new_input_variable(|| Ok(self.x4)).unwrap();
        let x5 = cs.new_input_variable(|| Ok(self.x5)).unwrap();
        let w1 = cs.new_witness_variable(|| Ok(self.w1)).unwrap();
        let w2 = cs.new_witness_variable(|| Ok(self.w2)).unwrap();
        let w3 = cs.new_witness_variable(|| Ok(self.w3)).unwrap();
        let w4 = cs.new_witness_variable(|| Ok(self.w4)).unwrap();
        let w5 = cs.new_witness_variable(|| Ok(self.w5)).unwrap();
        let w6 = cs.new_witness_variable(|| Ok(self.w6)).unwrap();
        let _w7 = cs.new_witness_variable(|| Ok(self.w7)).unwrap();
        let w8 = cs.new_witness_variable(|| Ok(self.w8)).unwrap();

        // Local predicate declarations -> Polynomial predicates

        let local_predicate_a = PredicateConstraintSystem::new_polynomial_predicate(
            cs.clone(),
            4,
            vec![
                (F::from(1u8), vec![(0, 1), (1, 1)]),
                (F::from(3u8), vec![(2, 2)]),
                (F::from(-1i8), vec![(3, 1)]),
            ],
        );
        let local_predicate_b = PredicateConstraintSystem::new_polynomial_predicate(
            cs.clone(),
            3,
            vec![
                (F::from(7u8), vec![(1, 1)]),
                (F::from(1u8), vec![(0, 3)]),
                (F::from(-1i8), vec![(2, 1)]),
            ],
        );

        let local_predicate_c = PredicateConstraintSystem::new_polynomial_predicate(
            cs.clone(),
            3,
            vec![
                (F::from(1u8), vec![(0, 1), (1, 1)]),
                (F::from(-1i8), vec![(2, 1)]),
            ],
        );
        cs.register_predicate("poly-predicate-A", local_predicate_a)?;
        cs.register_predicate("poly-predicate-B", local_predicate_b)?;
        cs.register_predicate("poly-predicate-C", local_predicate_c)?;

        // Enforing constraints to the predicates

        let predicate_a_constraint_1 = vec![lc!() + x1, lc!() + x2, lc!() + x3, lc!() + w4];
        let predicate_b_constraint_1 = vec![lc!() + x4, lc!() + w1, lc!() + w5];
        let predicate_c_constraint_1 = vec![lc!() + w2, lc!() + w3, lc!() + w6];
        let predicate_b_constraint_2 = vec![lc!() + w5, lc!() + w6, lc!() + w8];
        let predicate_c_constraint_2 = vec![lc!() + w5 + w4, lc!() + w8, lc!() + x5];

        cs.enforce_constraint("poly-predicate-A", predicate_a_constraint_1)?;
        cs.enforce_constraint("poly-predicate-B", predicate_b_constraint_1)?;
        cs.enforce_constraint("poly-predicate-C", predicate_c_constraint_1)?;
        cs.enforce_constraint("poly-predicate-B", predicate_b_constraint_2)?;
        cs.enforce_constraint("poly-predicate-C", predicate_c_constraint_2)?;
        Ok(())
    }
}
