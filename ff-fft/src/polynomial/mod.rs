//! Modules for working with univariate or multivariate polynomials.
use crate::Vec;
use algebra_core::Field;
use core::{
    fmt::Debug,
    hash::Hash,
    ops::{AddAssign, Index, SubAssign},
};
use rand::Rng;

pub mod multivariate;
pub mod univariate;

/// Describes the common interface for univariate and multivariate polynomials
pub trait Polynomial<F: Field>:
    Sized
    + Clone
    + Debug
    + Default
    + Hash
    + PartialEq
    + Eq
    + Send
    + Sync
    + for<'a> AddAssign<&'a Self>
    + for<'a> AddAssign<(F, &'a Self)>
    + for<'a> SubAssign<&'a Self>
{
    /// The domain of the polynomial.
    type Domain: Sized + Clone + Ord + Debug;

    /// Returns the zero polynomial.
    fn zero() -> Self;

    /// Checks if the given polynomial is zero.
    fn is_zero(&self) -> bool;

    /// Returns the total degree of the polynomial
    fn degree(&self) -> usize;

    /// Evaluates `self` at the given `point` in `Self::Domain`.
    fn evaluate(&self, point: &Self::Domain) -> F;

    /// If `num_vars` is `None`, outputs a polynomial a univariate polynomial
    /// of degree `d` where each coefficient is sampled uniformly at random.
    ///
    /// If `num_vars` is `Some(l)`, outputs an `l`-variate polynomial which
    /// is the sum of `l` `d`-degree univariate polynomials where each coefficient
    /// is sampled uniformly at random.
    fn rand<R: Rng>(d: usize, num_vars: Option<usize>, rng: &mut R) -> Self;

    /// Sample a random point from `Self::Domain`.  
    fn rand_domain_point<R: Rng>(domain_size: Option<usize>, rng: &mut R) -> Self::Domain;
}

/// Describes the interface for univariate polynomials
pub trait UVPolynomial<F: Field>: Polynomial<F> {
    /// Constructs a new polynomial from a list of coefficients.
    fn from_coefficients_slice(coeffs: &[F]) -> Self;

    /// Constructs a new polynomial from a list of coefficients.
    fn from_coefficients_vec(coeffs: Vec<F>) -> Self;

    /// Returns the coefficients of `self`
    fn coeffs(&self) -> &[F];
}

/// Describes the interface for univariate polynomials
pub trait MVPolynomial<F: Field>: Polynomial<F> {
    /// The type of the terms of `self`
    type Term: multivariate::Term;

    /// Constructs a new polynomial from a list of tuples of the form `(Self::Term, coeff)`
    fn from_coefficients_slice(num_vars: usize, terms: &[(Self::Term, F)]) -> Self {
        Self::from_coefficients_vec(num_vars, terms.to_vec())
    }

    /// Constructs a new polynomial from a list of tuples of the form `(Self::Term, coeff)`
    fn from_coefficients_vec(num_vars: usize, terms: Vec<(Self::Term, F)>) -> Self;

    /// Returns the terms of a `self` as a list of tuples of the form `(Self::Term, coeff)`
    fn terms(&self) -> &[(Self::Term, F)];

    /// Given some point `z`, compute the quotients `w_i(X)` s.t
    ///
    /// `p(X) - p(z) = (X_1-z_1)*w_1(X) + (X_2-z_2)*w_2(X) + ... + (X_l-z_l)*w_l(X)`
    ///
    /// These quotients can always be found with no remainder.
    fn divide_at_point(&self, point: &Self::Domain) -> Vec<Self>
    where
        Self::Domain: Index<usize, Output = F>;
}
