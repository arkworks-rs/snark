//! This module contains the implementation of the Polynomial Predicate struct.
//! A polynomial predicate is a kind of local predicate which is defined in https://eprint.iacr.org/2024/1245
//! Other kinds of local predicates can be added in the future such as lookup
//! table predicates.

use ark_ff::Field;
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial, Polynomial,
};
use ark_std::vec::Vec;

use super::Predicate;

/// A polynomial predicat is just a polynomial
#[derive(Debug, Clone)]
pub struct PolynomialPredicate<F: Field> {
    polynomial: SparsePolynomial<F, SparseTerm>,
}

impl<F: Field> PolynomialPredicate<F> {
    /// Create a new polynomial predicate with a given arity (number of
    /// variables) and terms
    pub fn new(arity: usize, terms: Vec<(F, Vec<(usize, usize)>)>) -> Self {
        let sparse_terms = terms
            .iter()
            .map(|(coeff, term)| (*coeff, SparseTerm::new(term.clone())))
            .collect();
        Self {
            polynomial: SparsePolynomial::from_coefficients_vec(arity, sparse_terms),
        }
    }
}

/// This is the implementation of the Predicate trait for PolynomialPredicate.
/// The evaluation of a polynomial predicate is the evaluation of the underlying
/// polynomial and the arity is the number of variables in the polynomial.
impl<F: Field> Predicate<F> for PolynomialPredicate<F> {
    fn evaluate(&self, variables: Vec<F>) -> bool {
        !self.polynomial.evaluate(&variables).is_zero()
    }

    fn arity(&self) -> usize {
        self.polynomial.num_vars()
    }
}

/// The label for the popular R1CS predicate
pub const R1CS_PREDICATE_LABEL: &str = "R1CS";
