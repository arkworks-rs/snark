//! A sparse multivariate polynomial represented in coefficient form.
use algebra_core::Field;
use core::{
    cmp::Ordering,
    fmt,
    ops::{Add, Deref, Neg, Sub},
};
use rand::Rng;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Represents a single term in a multivariate polynomial.
/// Variables are stored in ascending order.
#[derive(Clone, PartialOrd, PartialEq, Eq, Hash, Default)]
struct PolyVars(Vec<(usize, usize)>);

/// Stores a sparse multi-variate polynomial in coefficient form.
/// Terms are stored in ascending order of total degree
#[derive(Clone, Eq, Hash, Default)]
pub struct SparseMultiPolynomial<F: Field> {
    num_vars: usize,
    terms: Vec<(PolyVars, F)>,
}

impl<F: Field> SparseMultiPolynomial<F> {
    /// Return the number of variables
    pub fn num_vars(&self) -> usize {
        self.num_vars
    }

    /// Returns the zero polynomial.
    pub fn zero(num_vars: usize) -> Self {
        Self {
            num_vars,
            terms: Vec::new(),
        }
    }

    /// Checks if the given polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty() || self.terms.iter().all(|(_, c)| c.is_zero())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_slice(num_vars: usize, terms: &[(Vec<(usize, usize)>, F)]) -> Self {
        let terms: Vec<(PolyVars, F)> = terms
            .into_iter()
            .map(|c| (PolyVars::new((*c.0).to_vec()), c.1))
            .collect();
        Self::from_coefficients_vec(num_vars, terms)
    }

    fn from_coefficients_vec(num_vars: usize, mut terms: Vec<(PolyVars, F)>) -> Self {
        // Ensure that terms are in ascending order.
        terms.sort_by(|(c1, _), (c2, _)| c1.cmp(c2));
        // If any terms are duplicated, add them together
        let mut terms_dedup: Vec<(PolyVars, F)> = Vec::new();
        for term in terms {
            match terms_dedup.last_mut() {
                Some(prev) => {
                    if prev.0 == term.0 {
                        *prev = (prev.0.clone(), prev.1 + term.1);
                        continue;
                    }
                },
                _ => {},
            };
            // Assert correct number of indeterminates
            assert!(
                term.0.iter().all(|(var, _)| *var < num_vars),
                "Invalid number of indeterminates"
            );
            terms_dedup.push(term);
        }
        let mut result = Self {
            num_vars,
            terms: terms_dedup,
        };
        // Remove any terms with zero coefficients
        result.remove_zeros();
        result
    }

    fn remove_zeros(&mut self) {
        self.terms.retain(|(_, c)| c.is_zero());
    }

    /// Evaluates `self` at the given `point` in the field.
    pub fn evaluate(&self, point: &[F]) -> F {
        assert!(point.len() >= self.num_vars, "Invalid evaluation domain");
        if self.is_zero() {
            return F::zero();
        }
        cfg_into_iter!(&self.terms)
            .map(|(term, coeff)| *coeff * term.evaluate(point))
            .sum()
    }

    /// Outputs a polynomial which is the sum of `l` `d`-degree univariate
    /// polynomials where each coefficient is sampled uniformly at random
    /// from the field `F`.
    pub fn rand<R: Rng>(l: usize, d: usize, rng: &mut R) -> Self {
        let mut random_terms = Vec::new();
        for var in 0..l {
            for deg in 0..=d {
                random_terms.push((PolyVars(vec![(var, deg)]), F::rand(rng)));
            }
        }
        Self::from_coefficients_vec(l, random_terms)
    }
}

impl PolyVars {
    fn new(mut term: Vec<(usize, usize)>) -> Self {
        term.sort_by(|(v1, _), (v2, _)| v1.cmp(v2));
        Self(term)
    }

    /// Returns the sum of all variable powers in `self`
    fn total_degree(&self) -> usize {
        self.iter().fold(0, |sum, acc| sum + acc.1)
    }

    /// Evaluates `self` at the given `point` in the field.
    fn evaluate<F: Field>(&self, point: &[F]) -> F {
        cfg_into_iter!(self)
            .map(|var| {
                let term = point.get(var.0).unwrap();
                term.pow(&[var.1 as u64])
            })
            .product()
    }
}

impl<'a, 'b, F: Field> Add<&'a SparseMultiPolynomial<F>> for &'b SparseMultiPolynomial<F> {
    type Output = SparseMultiPolynomial<F>;

    fn add(self, other: &'a SparseMultiPolynomial<F>) -> SparseMultiPolynomial<F> {
        let mut result = Vec::new();
        let mut cur_iter = self.terms.iter().peekable();
        let mut other_iter = other.terms.iter().peekable();
        // Since both polynomials are sorted, iterate over them in ascending order,
        // combining any common terms
        loop {
            // Peek at iterators to decide which to take from
            let which = match (cur_iter.peek(), other_iter.peek()) {
                (Some(cur), Some(other)) => Some((cur.0).cmp(&other.0)),
                (Some(_), None) => Some(Ordering::Less),
                (None, Some(_)) => Some(Ordering::Greater),
                (None, None) => None,
            };
            // Push the smallest element to the `result` coefficient vec
            result.push(match which {
                Some(Ordering::Less) => cur_iter.next().unwrap().clone(),
                Some(Ordering::Equal) => {
                    let other = other_iter.next().unwrap();
                    let cur = cur_iter.next().unwrap();
                    (cur.0.clone(), cur.1 + other.1)
                },
                Some(Ordering::Greater) => other_iter.next().unwrap().clone(),
                None => break,
            });
        }
        return SparseMultiPolynomial::from_coefficients_vec(
            core::cmp::max(self.num_vars(), other.num_vars()),
            result,
        );
    }
}

impl<F: Field> Neg for SparseMultiPolynomial<F> {
    type Output = SparseMultiPolynomial<F>;

    #[inline]
    fn neg(mut self) -> SparseMultiPolynomial<F> {
        for coeff in &mut self.terms {
            (*coeff).1 = -coeff.1;
        }
        self
    }
}

impl<'a, 'b, F: Field> Sub<&'a SparseMultiPolynomial<F>> for &'b SparseMultiPolynomial<F> {
    type Output = SparseMultiPolynomial<F>;

    #[inline]
    fn sub(self, other: &'a SparseMultiPolynomial<F>) -> SparseMultiPolynomial<F> {
        let neg_other = other.clone().neg();
        self + &neg_other
    }
}

impl<F: Field> PartialEq for SparseMultiPolynomial<F> {
    fn eq(&self, other: &Self) -> bool {
        self.terms == other.terms
    }
}

impl<F: Field> fmt::Debug for SparseMultiPolynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        for (term, coeff) in self.terms.iter().filter(|(_, c)| !c.is_zero()) {
            if term.len() == 0 {
                write!(f, "\n{:?}", coeff)?;
            } else {
                write!(f, "\n{:?} {:?}", coeff, term)?;
            }
        }
        Ok(())
    }
}

impl fmt::Debug for PolyVars {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        for variable in self.iter() {
            if variable.1 == 1 {
                write!(f, " * x_{}", variable.0)?;
            } else {
                write!(f, " * x_{}^{}", variable.0, variable.1)?;
            }
        }
        Ok(())
    }
}

impl Deref for PolyVars {
    type Target = [(usize, usize)];

    fn deref(&self) -> &[(usize, usize)] {
        &self.0
    }
}

impl Ord for PolyVars {
    /// Sort by total degree. If total degree is equal then ordering
    /// is given by exponent weight in lower-numbered variables
    /// ie. `x_1 > x_2`, `x_1^2 > x_1 * x_2`, etc.
    fn cmp(&self, other: &Self) -> Ordering {
        if self.total_degree() != other.total_degree() {
            self.total_degree().cmp(&other.total_degree())
        } else {
            // Iterate through all variables and return the corresponding ordering
            // if they differ in variable numbering or power
            for (cur, other) in self.iter().zip(other.iter()) {
                if other.0 == cur.0 {
                    if cur.1 != other.1 {
                        return (cur.1).cmp(&other.1);
                    }
                } else {
                    return (other.0).cmp(&cur.0);
                }
            }
            Ordering::Equal
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::bls12_381::fr::Fr;
    use algebra_core::{test_rng, Field, UniformRand, Zero};

    /// Generate random `l`-variate polynomial of maximum individual degree `d`
    fn rand_poly<R: Rng>(l: usize, d: usize, rng: &mut R) -> SparseMultiPolynomial<Fr> {
        let mut random_terms = Vec::new();
        let num_terms = rng.gen_range(1, 1000);
        // For each term, randomly select up to `l` variables with degree
        // in [1,d] and random coefficient
        for _ in 0..num_terms {
            let term = (0..l)
                .map(|i| {
                    if rng.gen_bool(0.5) {
                        Some((i, rng.gen_range(1, d + 1)))
                    } else {
                        None
                    }
                })
                .filter(|t| t.is_some())
                .map(|t| t.unwrap())
                .collect();
            let coeff = Fr::rand(rng);
            random_terms.push((term, coeff));
        }
        SparseMultiPolynomial::from_coefficients_slice(l, random_terms.as_slice())
    }

    #[test]
    fn add_polynomials() {
        let rng = &mut test_rng();
        let max_degree = 10;
        for a_var_count in 1..20 {
            for b_var_count in 1..20 {
                let p1 = rand_poly(a_var_count, max_degree, rng);
                let p2 = rand_poly(b_var_count, max_degree, rng);
                let res1 = &p1 + &p2;
                let res2 = &p2 + &p1;
                assert_eq!(res1, res2);
            }
        }
    }

    #[test]
    fn sub_polynomials() {
        let rng = &mut test_rng();
        let max_degree = 10;
        for a_var_count in 1..20 {
            for b_var_count in 1..20 {
                let p1 = rand_poly(a_var_count, max_degree, rng);
                let p2 = rand_poly(b_var_count, max_degree, rng);
                let res1 = &p1 - &p2;
                let res2 = &p2 - &p1;
                assert_eq!(&res1 + &p2, p1);
                assert_eq!(res1, -res2.clone());
            }
        }
    }

    #[test]
    fn evaluate_polynomials() {
        let rng = &mut test_rng();
        let max_degree = 10;
        for var_count in 1..20 {
            let p = rand_poly(var_count, max_degree, rng);
            let mut point = Vec::with_capacity(var_count);
            for _ in 0..var_count {
                point.push(Fr::rand(rng));
            }
            let mut total = Fr::zero();
            for (term, coeff) in p.terms.iter() {
                let mut summand = *coeff;
                for var in term.iter() {
                    let eval = point.get(var.0).unwrap();
                    summand *= eval.pow(&[var.1 as u64]);
                }
                total += summand;
            }
            assert_eq!(p.evaluate(&point), total);
        }
    }

    #[test]
    fn add_and_evaluate_polynomials() {
        let rng = &mut test_rng();
        let max_degree = 10;
        for a_var_count in 1..20 {
            for b_var_count in 1..20 {
                let p1 = rand_poly(a_var_count, max_degree, rng);
                let p2 = rand_poly(b_var_count, max_degree, rng);
                let mut point = Vec::new();
                for _ in 0..core::cmp::max(a_var_count, b_var_count) {
                    point.push(Fr::rand(rng));
                }
                // Evaluate both polynomials at a given point
                let eval1 = p1.evaluate(&point);
                let eval2 = p2.evaluate(&point);
                // Add polynomials
                let sum = &p1 + &p2;
                // Evaluate result at same point
                let eval3 = sum.evaluate(&point);
                assert_eq!(eval1 + eval2, eval3);
            }
        }
    }
}
