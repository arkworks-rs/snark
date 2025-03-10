//! This module contains the implementation of the Polynomial Predicate struct.
//! A polynomial predicate is a kind of  predicate which is defined in https://eprint.iacr.org/2024/1245
//! Other kinds of  predicates can be added in the future such as lookup
//! table predicates.

use ark_ff::Field;
use ark_poly::{
    multivariate::{SparsePolynomial, SparseTerm, Term},
    DenseMVPolynomial, Polynomial,
};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::vec::Vec;

/// A polynomial predicat is just a polynomial
#[derive(Debug, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct PolynomialPredicate<F: Field> {
    /// The sparse polynomial for the predicate
    pub polynomial: SparsePolynomial<F, SparseTerm>,
}

impl<F: Field> PolynomialPredicate<F> {
    /// Create a new polynomial predicate with a given arity (number of
    /// variables) and terms
    pub fn new(arity: usize, terms: impl IntoIterator<Item = (F, Vec<(usize, usize)>)>) -> Self {
        let sparse_terms = terms
            .into_iter()
            .map(|(coeff, term)| (coeff, SparseTerm::new(term.clone())))
            .collect();
        Self {
            polynomial: SparsePolynomial::from_coefficients_vec(arity, sparse_terms),
        }
    }
}

/// This is the implementation of the Predicate trait for PolynomialPredicate.
/// The evaluation of a polynomial predicate is the evaluation of the underlying
/// polynomial and the arity is the number of variables in the polynomial.
impl<F: Field> PolynomialPredicate<F> {
    /// Check if the predicate is satisfied by the given variables
    pub fn is_satisfied(&self, variables: &[F]) -> bool {
        !self.polynomial.evaluate(variables).is_zero()
    }
    /// What is the evaluation of the polynomial predicate given the variables
    pub fn eval(&self, variables: &[F]) -> F {
        self.polynomial.evaluate(variables)
    }

    /// Get the arity of the polynomial predicate
    /// The arity of P(x1, x2, ..., xn) is n
    pub fn arity(&self) -> usize {
        self.polynomial.num_vars()
    }

    /// Get the degree of the polynomial predicate
    /// The degree of x1 + x2^4 + x3^2 is 4
    pub fn degree(&self) -> usize {
        self.polynomial.degree()
    }
}

/// The label for the popular R1CS predicate
pub const R1CS_PREDICATE_LABEL: &str = "R1CS";
