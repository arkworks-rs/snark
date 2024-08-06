pub mod lookup_constraint;
pub mod polynomial_constraint;

use core::{fmt::Debug, ops::Not, pin::Pin};

use ark_ff::Field;

use super::{
    constraint_system, Constraint, ConstraintSystemRef, LcIndex, LinearCombination, Matrix,
};
use crate::utils::{error::SynthesisError::ArityMismatch, variable::Variable::SymbolicLc};
use ark_std::{boxed::Box, rc::Rc, vec::Vec};
use lookup_constraint::LookupConstraint;
use polynomial_constraint::PolynomialConstraint;

// use crate::{gr1cs::ConstraintSystemRef, utils::Result};
pub trait Evaluatable<F> {
    fn evaluate(&self, point: Vec<F>) -> F;
}

#[derive(Debug, Clone)]
pub enum ConstraintType<F: Field>
where
    PolynomialConstraint<F>: Evaluatable<F>,
    LookupConstraint<F>: Evaluatable<F>,
{
    Polynomial(PolynomialConstraint<F>),
    Lookup(LookupConstraint<F>),
}

impl<F: Field> Evaluatable<F> for ConstraintType<F> {
    fn evaluate(&self, point: Vec<F>) -> F {
        match self {
            ConstraintType::Polynomial(p) => p.evaluate(point),
            ConstraintType::Lookup(l) => l.evaluate(point),
        }
    }
}

#[derive(Debug, Clone)]
pub struct LocalPredicate<F: Field> {
    arity: usize,
    cs: ConstraintSystemRef<F>,
    constraints: Vec<Constraint>,
    constraint_type: ConstraintType<F>,
}

impl<F: Field> LocalPredicate<F> {
    fn new(cs: ConstraintSystemRef<F>, arity: usize, constraint_type: ConstraintType<F>) -> Self {
        Self {
            arity,
            cs,
            constraints: Vec::new(),
            constraint_type,
        }
    }

    pub fn new_polynomial_predicate(
        cs: ConstraintSystemRef<F>,
        arity: usize,
        terms: Vec<(F, Vec<(usize, usize)>)>,
    ) -> Self {
        Self::new(
            cs,
            arity,
            ConstraintType::Polynomial(PolynomialConstraint::new(arity, terms)),
        )
    }
    pub fn get_arity(&self) -> usize {
        self.arity
    }

    pub fn num_constraints(&self) -> usize {
        self.constraints.len()
    }

    pub fn get_constraints(&self) -> &Vec<Constraint> {
        &self.constraints
    }

    pub fn get_cs(&self) -> ConstraintSystemRef<F> {
        self.cs.clone()
    }

    pub fn get_constraint_type(&self) -> &ConstraintType<F> {
        &self.constraint_type
    }

    pub fn enforce_constraint(&mut self, constraint: Constraint) -> crate::utils::Result<()> {
        if constraint.len() != self.arity {
            return Err(ArityMismatch);
        }
        self.constraints.push(constraint);
        Ok(())
    }

    pub fn which_constraint_is_unsatisfied(&self) -> Option<usize> {
        for (i, constraint) in self.constraints.iter().enumerate() {
            let point: Vec<F> = constraint
                .iter()
                .map(|lc_index| self.cs.assigned_value(SymbolicLc(*lc_index)).unwrap())
                .collect();
            let result = self.constraint_type.evaluate(point);
            if !result.is_zero() {
                return Some(i);
            }
        }
        None
    }
    pub fn to_matrices(&self) -> Vec<Matrix<F>> {
        let mut matrices: Vec<Matrix<F>> = vec![Vec::new(); self.arity];
        for constraint in self.constraints.iter() {
            for (matrix_ind, lc_index) in constraint.iter().enumerate() {
                let lc: LinearCombination<F> = self.cs.get_lc(*lc_index).unwrap();
                let row: Vec<(F, usize)> = self.make_row(&lc);
                matrices[matrix_ind].push(row);
            }
        }

        matrices
    }

    #[inline]
    fn make_row(&self, l: &LinearCombination<F>) -> Vec<(F, usize)> {
        let num_input = self.cs.num_instance_variables();
        l.0.iter()
            .filter_map(|(coeff, var)| {
                if coeff.is_zero() {
                    None
                } else {
                    Some((
                        *coeff,
                        var.get_index_unchecked(num_input).expect("no symbolic LCs"),
                    ))
                }
            })
            .collect()
    }

    pub fn build_r1cs_predicate(cs: ConstraintSystemRef<F>) -> Self {
        Self::new_polynomial_predicate(
            cs,
            3,
            vec![
                (F::from(1u8), vec![(0, 1),(1,1)]),
                (F::from(-1i8), vec![(2, 1)]),
            ],
        )
    }
}
