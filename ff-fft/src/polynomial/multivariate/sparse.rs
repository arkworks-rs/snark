//! A sparse multivariate polynomial represented in coefficient form.
use crate::{
    multivariate::{SparseTerm, Term},
    MVPolynomial, Polynomial, Vec,
};
use algebra_core::Field;
use core::{
    cmp::Ordering,
    fmt,
    ops::{Add, AddAssign, Neg, Sub, SubAssign},
};
use rand::Rng;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Stores a sparse multivariate polynomial in coefficient form.
#[derive(Derivative)]
#[derivative(Clone, PartialEq, Eq, Hash, Default)]
pub struct SparsePolynomial<F: Field, T: Term> {
    /// The number of variables the polynomial supports
    #[derivative(PartialEq = "ignore")]
    pub num_vars: usize,
    /// List of every term along with its coefficient
    pub terms: Vec<(T, F)>,
}

impl<F: Field, T: Term> SparsePolynomial<F, T> {
    fn remove_zeros(&mut self) {
        self.terms.retain(|(_, c)| !c.is_zero());
    }
}

impl<F: Field> Polynomial<F> for SparsePolynomial<F, SparseTerm> {
    type Domain = Vec<F>;

    /// Returns the zero polynomial.
    fn zero() -> Self {
        Self {
            num_vars: 0,
            terms: Vec::new(),
        }
    }

    /// Checks if the given polynomial is zero.
    fn is_zero(&self) -> bool {
        self.terms.is_empty() || self.terms.iter().all(|(_, c)| c.is_zero())
    }

    /// Returns the total degree of the polynomial
    fn degree(&self) -> usize {
        self.terms
            .iter()
            .map(|(term, _)| (*term).degree())
            .max()
            .unwrap_or(0)
    }

    /// Evaluates `self` at the given `point` in `Self::Domain`.
    fn evaluate(&self, point: &Vec<F>) -> F {
        assert!(point.len() >= self.num_vars, "Invalid evaluation domain");
        if self.is_zero() {
            return F::zero();
        }
        cfg_into_iter!(&self.terms)
            .map(|(term, coeff)| *coeff * term.evaluate(point))
            .sum()
    }

    /// Outputs an `l`-variate polynomial which is the sum of `l` `d`-degree
    /// univariate polynomials where each coefficient is sampled uniformly at random.
    fn rand<R: Rng>(d: usize, l: Option<usize>, rng: &mut R) -> Self {
        let l = l.unwrap();
        let mut random_terms = Vec::new();
        random_terms.push((SparseTerm::new(vec![]), F::rand(rng)));
        for var in 0..l {
            for deg in 1..=d {
                random_terms.push((SparseTerm::new(vec![(var, deg)]), F::rand(rng)));
            }
        }
        Self::from_coefficients_vec(l, random_terms)
    }

    /// Sample a random point from `Self::Domain`.  
    fn rand_domain_point<R: Rng>(domain_size: Option<usize>, rng: &mut R) -> Vec<F> {
        let mut point = Vec::with_capacity(domain_size.unwrap());
        for _ in 0..domain_size.unwrap() {
            point.push(F::rand(rng))
        }
        point
    }
}

impl<F: Field> MVPolynomial<F> for SparsePolynomial<F, SparseTerm> {
    type Term = SparseTerm;

    /// Constructs a new polynomial from a list of tuples of the form `(Self::Term, coeff)`
    fn from_coefficients_vec(num_vars: usize, mut terms: Vec<(SparseTerm, F)>) -> Self {
        // Ensure that terms are in ascending order.
        terms.sort_by(|(c1, _), (c2, _)| c1.cmp(c2));
        // If any terms are duplicated, add them together
        let mut terms_dedup: Vec<(SparseTerm, F)> = Vec::new();
        for term in terms {
            match terms_dedup.last_mut() {
                Some(prev) => {
                    if prev.0 == term.0 {
                        *prev = (prev.0.clone(), prev.1 + term.1);
                        continue;
                    }
                }
                _ => {}
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

    /// Returns the terms of a `self` as a list of tuples of the form `(Self::Term, coeff)`
    fn terms(&self) -> &[(Self::Term, F)] {
        self.terms.as_slice()
    }

    /// Given some point `z`, compute the quotients `w_i(X)` s.t
    ///
    /// `p(X) - p(z) = (X_1-z_1)*w_1(X) + (X_2-z_2)*w_2(X) + ... + (X_l-z_l)*w_l(X)`
    ///
    /// These quotients can always be found with no remainder.
    fn divide_at_point(&self, point: &Self::Domain) -> Vec<Self> {
        if self.is_zero() {
            return vec![SparsePolynomial::zero(); self.num_vars];
        }
        assert_eq!(point.len(), self.num_vars, "Invalid evaluation point");
        let mut quotients = Vec::with_capacity(self.num_vars);
        // `cur` represents the current dividend
        let mut cur = self.clone();
        // Divide `cur` by `X_i - z_i`
        for i in 0..self.num_vars {
            let mut quotient_terms = Vec::new();
            let mut remainder_terms = Vec::new();
            for (mut term, mut coeff) in cur.terms {
                // Since the final remainder is guaranteed to be 0, all the constant terms
                // cancel out so we don't need to keep track of them
                if term.is_constant() {
                    continue;
                }
                // If the current term contains `X_i` then divide appropiately,
                // otherwise add it to the remainder
                match term.binary_search_by(|(var, _)| var.cmp(&i)) {
                    Ok(idx) => {
                        // Repeatedly divide the term by `X_i - z_i` until the remainder
                        // doesn't contain any `X_i`s
                        while term[idx].1 > 1 {
                            // First divide by `X_i` and add the term to the quotient
                            term.0[idx] = (i, term[idx].1 - 1);
                            quotient_terms.push((term.clone(), coeff));
                            // Then compute the remainder term in-place
                            coeff *= &point[i];
                        }
                        // Since `X_i` is power 1, we can remove it entirely
                        term.0.remove(idx);
                        quotient_terms.push((term.clone(), coeff));
                        remainder_terms.push((term, point[i] * coeff));
                    }
                    Err(_) => remainder_terms.push((term, coeff)),
                }
            }
            quotients.push(SparsePolynomial::from_coefficients_vec(
                self.num_vars,
                quotient_terms,
            ));
            // Set the current dividend to be the remainder of this division
            cur = SparsePolynomial::from_coefficients_vec(self.num_vars, remainder_terms);
        }
        quotients
    }
}

impl<'a, 'b, F: Field, T: Term> Add<&'a SparsePolynomial<F, T>> for &'b SparsePolynomial<F, T> {
    type Output = SparsePolynomial<F, T>;

    fn add(self, other: &'a SparsePolynomial<F, T>) -> SparsePolynomial<F, T> {
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
                }
                Some(Ordering::Greater) => other_iter.next().unwrap().clone(),
                None => break,
            });
        }
        // Remove any zero terms
        result.retain(|(_, c)| !c.is_zero());
        SparsePolynomial {
            num_vars: core::cmp::max(self.num_vars, other.num_vars),
            terms: result,
        }
    }
}

impl<'a, 'b, F: Field, T: Term> AddAssign<&'a SparsePolynomial<F, T>> for SparsePolynomial<F, T> {
    fn add_assign(&mut self, other: &'a SparsePolynomial<F, T>) {
        *self = &*self + other;
    }
}

impl<'a, 'b, F: Field, T: Term> AddAssign<(F, &'a SparsePolynomial<F, T>)>
    for SparsePolynomial<F, T>
{
    fn add_assign(&mut self, (f, other): (F, &'a SparsePolynomial<F, T>)) {
        let other = Self {
            num_vars: other.num_vars,
            terms: other
                .terms
                .iter()
                .map(|(term, coeff)| (term.clone(), *coeff * &f))
                .collect(),
        };
        *self = &*self + &other;
        self.remove_zeros()
    }
}

impl<F: Field, T: Term> Neg for SparsePolynomial<F, T> {
    type Output = SparsePolynomial<F, T>;

    #[inline]
    fn neg(mut self) -> SparsePolynomial<F, T> {
        for coeff in &mut self.terms {
            (*coeff).1 = -coeff.1;
        }
        self
    }
}

impl<'a, 'b, F: Field, T: Term> Sub<&'a SparsePolynomial<F, T>> for &'b SparsePolynomial<F, T> {
    type Output = SparsePolynomial<F, T>;

    #[inline]
    fn sub(self, other: &'a SparsePolynomial<F, T>) -> SparsePolynomial<F, T> {
        let neg_other = other.clone().neg();
        self + &neg_other
    }
}

impl<'a, 'b, F: Field, T: Term> SubAssign<&'a SparsePolynomial<F, T>> for SparsePolynomial<F, T> {
    #[inline]
    fn sub_assign(&mut self, other: &'a SparsePolynomial<F, T>) {
        *self = &*self - other;
    }
}

impl<F: Field, T: Term> fmt::Debug for SparsePolynomial<F, T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        for (term, coeff) in self.terms.iter().filter(|(_, c)| !c.is_zero()) {
            if term.is_constant() {
                write!(f, "\n{:?}", coeff)?;
            } else {
                write!(f, "\n{:?} {:?}", coeff, term)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::bls12_381::fr::Fr;
    use algebra_core::{test_rng, Field, One, UniformRand, Zero};

    // TODO: Make tests generic over term type

    /// Generate random `l`-variate polynomial of maximum individual degree `d`
    fn rand_poly<R: Rng>(l: usize, d: usize, rng: &mut R) -> SparsePolynomial<Fr, SparseTerm> {
        let mut random_terms = Vec::new();
        let num_terms = rng.gen_range(1, 1000);
        // For each term, randomly select up to `l` variables with degree
        // in [1,d] and random coefficient
        random_terms.push((SparseTerm::new(vec![]), Fr::rand(rng)));
        for _ in 1..num_terms {
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
            random_terms.push((SparseTerm::new(term), coeff));
        }
        SparsePolynomial::from_coefficients_slice(l, &random_terms)
    }

    /// Perform a naive n^2 multiplication of `self` by `other`.
    fn naive_mul(
        cur: &SparsePolynomial<Fr, SparseTerm>,
        other: &SparsePolynomial<Fr, SparseTerm>,
    ) -> SparsePolynomial<Fr, SparseTerm> {
        if cur.is_zero() || other.is_zero() {
            SparsePolynomial::zero()
        } else {
            let mut result_terms = Vec::new();
            for (cur_term, cur_coeff) in cur.terms.iter() {
                for (other_term, other_coeff) in other.terms.iter() {
                    let mut term = cur_term.0.clone();
                    term.extend(other_term.0.clone());
                    result_terms.push((SparseTerm::new(term), *cur_coeff * *other_coeff));
                }
            }
            SparsePolynomial::from_coefficients_slice(cur.num_vars, result_terms.as_slice())
        }
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
                assert_eq!(res1, -res2);
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

    #[test]
    /// This is just to make sure naive_mul works as expected
    fn mul_polynomials_fixed() {
        let a = SparsePolynomial::from_coefficients_slice(
            4,
            &[
                (SparseTerm(vec![]), "2".parse().unwrap()),
                (SparseTerm(vec![(0, 1), (1, 2)]), "4".parse().unwrap()),
                (SparseTerm(vec![(0, 1), (0, 1)]), "8".parse().unwrap()),
                (SparseTerm(vec![(3, 0)]), "1".parse().unwrap()),
            ],
        );
        let b = SparsePolynomial::from_coefficients_slice(
            4,
            &[
                (SparseTerm(vec![(0, 1), (1, 2)]), "1".parse().unwrap()),
                (SparseTerm(vec![(2, 1)]), "2".parse().unwrap()),
                (SparseTerm(vec![(3, 1)]), "1".parse().unwrap()),
            ],
        );
        let result = naive_mul(&a, &b);
        let expected = SparsePolynomial::from_coefficients_slice(
            4,
            &[
                (SparseTerm(vec![(0, 1), (1, 2)]), "3".parse().unwrap()),
                (SparseTerm(vec![(2, 1)]), "6".parse().unwrap()),
                (SparseTerm(vec![(3, 1)]), "3".parse().unwrap()),
                (SparseTerm(vec![(0, 2), (1, 4)]), "4".parse().unwrap()),
                (
                    SparseTerm(vec![(0, 1), (1, 2), (2, 1)]),
                    "8".parse().unwrap(),
                ),
                (
                    SparseTerm(vec![(0, 1), (1, 2), (3, 1)]),
                    "4".parse().unwrap(),
                ),
                (SparseTerm(vec![(0, 3), (1, 2)]), "8".parse().unwrap()),
                (SparseTerm(vec![(0, 2), (2, 1)]), "16".parse().unwrap()),
                (SparseTerm(vec![(0, 2), (3, 1)]), "8".parse().unwrap()),
            ],
        );
        assert_eq!(expected, result);
    }

    #[test]
    fn divide_at_point_fixed() {
        let dividend = SparsePolynomial::<Fr, _>::from_coefficients_slice(
            3,
            &[
                (SparseTerm(vec![]), "4".parse().unwrap()),
                (SparseTerm(vec![(0, 1)]), "1".parse().unwrap()),
                (SparseTerm(vec![(0, 1), (1, 1)]), "3".parse().unwrap()),
                (
                    SparseTerm(vec![(0, 2), (1, 2), (2, 3)]),
                    "5".parse().unwrap(),
                ),
            ],
        );
        let point = vec![
            "4".parse().unwrap(),
            "2".parse().unwrap(),
            "1".parse().unwrap(),
        ];
        let result = dividend.divide_at_point(&point);
        let expected_result = vec![
            SparsePolynomial::from_coefficients_slice(
                3,
                &[
                    (SparseTerm(vec![]), "1".parse().unwrap()),
                    (SparseTerm(vec![(1, 1)]), "3".parse().unwrap()),
                    (SparseTerm(vec![(1, 2), (2, 3)]), "20".parse().unwrap()),
                    (
                        SparseTerm(vec![(0, 1), (1, 2), (2, 3)]),
                        "5".parse().unwrap(),
                    ),
                ],
            ),
            SparsePolynomial::from_coefficients_slice(
                3,
                &[
                    (SparseTerm(vec![]), "12".parse().unwrap()),
                    (SparseTerm(vec![(2, 3)]), "160".parse().unwrap()),
                    (SparseTerm(vec![(1, 1), (2, 3)]), "80".parse().unwrap()),
                ],
            ),
            SparsePolynomial::from_coefficients_slice(
                3,
                &[
                    (SparseTerm(vec![]), "320".parse().unwrap()),
                    (SparseTerm(vec![(2, 1)]), "320".parse().unwrap()),
                    (SparseTerm(vec![(2, 2)]), "320".parse().unwrap()),
                ],
            ),
        ];
        assert_eq!(expected_result, result);
    }

    #[test]
    fn divide_at_point_random() {
        let rng = &mut test_rng();
        let max_degree = 10;
        for var_count in 0..20 {
            let dividend = SparsePolynomial::rand(max_degree, Some(var_count), rng);
            let mut point = Vec::new();
            for _ in 0..var_count {
                point.push(Fr::rand(rng));
            }
            let dividend_eval = SparsePolynomial::from_coefficients_slice(
                0,
                &vec![(SparseTerm(vec![]), dividend.evaluate(&point))],
            );
            let quotients = dividend.divide_at_point(&point);
            let mut result = SparsePolynomial::zero();
            for (i, q) in quotients.iter().enumerate() {
                let test_terms = vec![
                    (SparseTerm(vec![(i, 1)]), Fr::one()),
                    (SparseTerm(vec![]), -point[i]),
                ];
                let test = SparsePolynomial::from_coefficients_slice(var_count, &test_terms);
                result += &naive_mul(q, &test);
            }
            assert_eq!(&dividend - &dividend_eval, result);
        }
    }
}
