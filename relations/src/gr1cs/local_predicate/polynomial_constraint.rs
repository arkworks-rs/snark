use ark_ff::Field;
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial, Polynomial,
};
use ark_std::vec::Vec;

use crate::gr1cs::{Constraint, ConstraintSystemRef};

use super::Evaluatable;

#[derive(Debug, Clone)]
pub struct PolynomialConstraint<F: Field> {
    polynomial: SparsePolynomial<F, SparseTerm>,
}

impl<F: Field> PolynomialConstraint<F> {
    pub fn new(arity: usize, terms: Vec<(F, Vec<(usize, usize)>)>) -> Self {
        let sparse_terms = terms
            .iter()
            .map(|(coeff, term)| (coeff.clone(), SparseTerm::new(term.clone())))
            .collect();
        Self {
            polynomial: SparsePolynomial::from_coefficients_vec(arity, sparse_terms),
        }
    }
}

impl<F: Field> Evaluatable<F> for PolynomialConstraint<F> {
    fn evaluate(&self, point: Vec<F>) -> F {
        self.polynomial.evaluate(&point)
    }
}
pub const R1CS_PREDICATE_LABEL: &str = "R1CS";
