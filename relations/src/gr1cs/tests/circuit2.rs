use ark_ff::Field;

use ark_std::{collections::BTreeMap, string::ToString, vec::Vec};

use crate::{
    gr1cs::{ConstraintSynthesizer, ConstraintSystemRef},
    lc,
};

use super::{
    predicate::polynomial_constraint::R1CS_PREDICATE_LABEL, Label, Matrix, Variable,
};

#[derive(Debug, Clone)]
pub struct Circuit2<F: Field> {
    pub a: F,
    pub b: F,
    pub c: F,
}

impl<F: Field> Circuit2<F> {
    pub fn get_matrices() -> BTreeMap<Label, Vec<Matrix<F>>> {
        let two = F::one() + F::one();
        let mut map: BTreeMap<Label, Vec<Matrix<F>>> = BTreeMap::new();
        map.insert(
            R1CS_PREDICATE_LABEL.to_string(),
            vec![
                vec![
                    vec![(F::one(), 1)],
                    vec![(F::one(), 1)],
                    vec![(F::one(), 0)],
                ],
                vec![
                    vec![(two, 2)],
                    vec![(F::one(), 1), (F::one(), 2)],
                    vec![(two, 1), (two, 2)],
                ],
                vec![
                    vec![(F::one(), 3)],
                    vec![(F::one(), 1), (F::one(), 2)],
                    vec![(two, 1), (two, 2)],
                ],
            ],
        );

        map
    }
}
impl<F: Field + core::convert::From<i8>> ConstraintSynthesizer<F> for Circuit2<F> {
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> crate::utils::Result<()> {
        // Variable declarations -> Instance variables + Witness variables
        let two = F::one() + F::one();
        let a = cs.new_input_variable(|| Ok(self.a)).unwrap();
        let b = cs.new_witness_variable(|| Ok(self.b)).unwrap();
        let c = cs.new_witness_variable(|| Ok(self.c)).unwrap();
        cs.enforce_r1cs_constraint(lc!() + a, lc!() + (two, b), lc!() + c)?;
        let d = cs.new_lc(lc!() + a + b)?;
        cs.enforce_r1cs_constraint(lc!() + a, lc!() + d, lc!() + d)?;
        let e = cs.new_lc(lc!() + d + d)?;
        cs.enforce_r1cs_constraint(lc!() + Variable::One, lc!() + e, lc!() + e)?;

        Ok(())
    }
}
