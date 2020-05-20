//! A sparse polynomial represented in coefficient form.

use core::fmt;

use crate::{
    BTreeMap, DenseOrSparsePolynomial, DensePolynomial, EvaluationDomain, Evaluations, Vec,
};
use algebra_core::{FftField, Field};

/// Stores a sparse polynomial in coefficient form.
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct SparsePolynomial<F: Field> {
    /// The coefficient a_i of `x^i` is stored as (i, a_i) in `self.coeffs`.
    /// the entries in `self.coeffs` *must*  be sorted in increasing order of
    /// `i`.
    coeffs: Vec<(usize, F)>,
}

impl<F: Field> fmt::Debug for SparsePolynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        for (i, coeff) in self.coeffs.iter().filter(|(_, c)| !c.is_zero()) {
            if *i == 0 {
                write!(f, "\n{:?}", coeff)?;
            } else if *i == 1 {
                write!(f, " + \n{:?} * x", coeff)?;
            } else {
                write!(f, " + \n{:?} * x^{}", coeff, i)?;
            }
        }
        Ok(())
    }
}

impl<F: Field> core::ops::Deref for SparsePolynomial<F> {
    type Target = [(usize, F)];

    fn deref(&self) -> &[(usize, F)] {
        &self.coeffs
    }
}

impl<F: Field> SparsePolynomial<F> {
    /// Returns the zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    /// Checks if the given polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || self.coeffs.iter().all(|(_, c)| c.is_zero())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_slice(coeffs: &[(usize, F)]) -> Self {
        Self::from_coefficients_vec(coeffs.to_vec())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_vec(mut coeffs: Vec<(usize, F)>) -> Self {
        // While there are zeros at the end of the coefficient vector, pop them off.
        while coeffs.last().map_or(false, |(_, c)| c.is_zero()) {
            coeffs.pop();
        }
        // Ensure that coeffs are in ascending order.
        coeffs.sort_by(|(c1, _), (c2, _)| c1.cmp(c2));
        // Check that either the coefficients vec is empty or that the last coeff is
        // non-zero.
        assert!(coeffs.last().map_or(true, |(_, c)| !c.is_zero()));

        Self { coeffs }
    }

    /// Returns the degree of the polynomial.
    pub fn degree(&self) -> usize {
        if self.is_zero() {
            0
        } else {
            assert!(self.coeffs.last().map_or(false, |(_, c)| !c.is_zero()));
            self.coeffs.last().unwrap().0
        }
    }

    /// Evaluates `self` at the given `point` in the field.
    pub fn evaluate(&self, point: F) -> F {
        if self.is_zero() {
            return F::zero();
        }
        let mut total = F::zero();
        for (i, c) in &self.coeffs {
            total += &(*c * &point.pow(&[*i as u64]));
        }
        total
    }

    /// Perform a naive n^2 multiplicatoin of `self` by `other`.
    pub fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            SparsePolynomial::zero()
        } else {
            let mut result = BTreeMap::new();
            for (i, self_coeff) in self.coeffs.iter() {
                for (j, other_coeff) in other.coeffs.iter() {
                    let cur_coeff = result.entry(i + j).or_insert(F::zero());
                    *cur_coeff += &(*self_coeff * other_coeff);
                }
            }
            let result = result.into_iter().collect::<Vec<_>>();
            SparsePolynomial::from_coefficients_vec(result)
        }
    }
}

impl<F: FftField> SparsePolynomial<F> {
    /// Evaluate `self` over `domain`.
    pub fn evaluate_over_domain_by_ref<D: EvaluationDomain<F>>(
        &self,
        domain: D,
    ) -> Evaluations<F, D> {
        let poly: DenseOrSparsePolynomial<'_, F> = self.into();
        DenseOrSparsePolynomial::<F>::evaluate_over_domain(poly, domain)
    }

    /// Evaluate `self` over `domain`.
    pub fn evaluate_over_domain<D: EvaluationDomain<F>>(self, domain: D) -> Evaluations<F, D> {
        let poly: DenseOrSparsePolynomial<'_, F> = self.into();
        DenseOrSparsePolynomial::<F>::evaluate_over_domain(poly, domain)
    }
}

impl<F: Field> Into<DensePolynomial<F>> for SparsePolynomial<F> {
    fn into(self) -> DensePolynomial<F> {
        let mut other = vec![F::zero(); self.degree() + 1];
        for (i, coeff) in self.coeffs {
            other[i] = coeff;
        }
        DensePolynomial::from_coefficients_vec(other)
    }
}

#[cfg(test)]
mod tests {
    use crate::{DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, SparsePolynomial};
    use algebra::bls12_381::fr::Fr;
    use algebra_core::One;

    #[test]
    fn evaluate_over_domain() {
        for size in 2..10 {
            let domain_size = 1 << size;
            let domain = GeneralEvaluationDomain::new(domain_size).unwrap();
            let two = Fr::one() + &Fr::one();
            let sparse_poly = SparsePolynomial::from_coefficients_vec(vec![(0, two), (1, two)]);
            let evals1 = sparse_poly.evaluate_over_domain_by_ref(domain);

            let dense_poly: DensePolynomial<Fr> = sparse_poly.into();
            let evals2 = dense_poly.clone().evaluate_over_domain(domain);
            assert_eq!(evals1.clone().interpolate(), evals2.clone().interpolate());
            assert_eq!(evals1.interpolate(), dense_poly);
            assert_eq!(evals2.interpolate(), dense_poly);
        }
    }
}
