//! This module contains the implementation of the Polynomial Predicate struct.
//! A polynomial predicate is a kind of predicate which is defined in <https://eprint.iacr.org/2024/1245>.
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
    /// Create a new polynomial predicate with a given number of variables (= arity) and
    /// number of terms.
    ///
    /// Here `terms` is a list of tuples of the form `(variable`, `power`) where
    /// `variable` is the index of the variable and `power` is the exponent of the variable.
    ///
    /// So, for example, if `terms` is `[F::ONE, [(0, 1), (1, 2)]]`, then the polynomial is
    /// `x_0 + x_1^2`.
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
    /// Check if the predicate is satisfied by the given variables.
    pub fn is_satisfied(&self, variables: &[F]) -> bool {
        // TODO: Change the polynomial eval to get a slice as an evaluation point
        !self.polynomial.evaluate(&variables.to_vec()).is_zero()
    }

    /// Evaluate the predicate on the given variables.
    pub fn eval(&self, variables: &[F]) -> F {
        self.polynomial.evaluate(&variables.to_vec())
    }

    /// Get the arity of the polynomial predicate.
    /// For example, the arity of P(x1, x2, ..., xn) is n.
    pub fn arity(&self) -> usize {
        self.polynomial.num_vars()
    }

    /// Get the degree of the polynomial predicate.
    /// For example, the degree of x1 + x2^4 + x3^2 is 4.
    pub fn degree(&self) -> usize {
        self.polynomial.degree()
    }
}

/// A label for the R1CS predicate.
pub const R1CS_PREDICATE_LABEL: &str = "R1CS";

/// A label for the Square R1CS predicate introduced in
/// [\[Groth-Maller17\]](https://eprint.iacr.org/2017/540).
pub const SR1CS_PREDICATE_LABEL: &str = "SR1CS";
