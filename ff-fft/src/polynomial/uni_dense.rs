//! A polynomial represented in coefficient form.

use crate::{GeneralEvaluationDomain, Vec};
use core::{
    fmt,
    ops::{Add, AddAssign, Deref, DerefMut, Div, Mul, Neg, Sub, SubAssign},
};

use crate::{DenseOrSparseUniPolynomial, EvaluationDomain, Evaluations};
use algebra_core::{FftField, Field};
use rand::Rng;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Stores a polynomial in coefficient form.
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct DenseUniPolynomial<F: Field> {
    /// The coefficient of `x^i` is stored at location `i` in `self.coeffs`.
    pub coeffs: Vec<F>,
}

impl<F: Field> fmt::Debug for DenseUniPolynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        for (i, coeff) in self.coeffs.iter().enumerate().filter(|(_, c)| !c.is_zero()) {
            if i == 0 {
                write!(f, "\n{:?}", coeff)?;
            } else if i == 1 {
                write!(f, " + \n{:?} * x", coeff)?;
            } else {
                write!(f, " + \n{:?} * x^{}", coeff, i)?;
            }
        }
        Ok(())
    }
}

impl<F: Field> Deref for DenseUniPolynomial<F> {
    type Target = [F];

    fn deref(&self) -> &[F] {
        &self.coeffs
    }
}

impl<F: Field> DerefMut for DenseUniPolynomial<F> {
    fn deref_mut(&mut self) -> &mut [F] {
        &mut self.coeffs
    }
}

impl<F: Field> DenseUniPolynomial<F> {
    /// Returns the zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    /// Checks if the given polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || self.coeffs.iter().all(|coeff| coeff.is_zero())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_slice(coeffs: &[F]) -> Self {
        Self::from_coefficients_vec(coeffs.to_vec())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_vec(coeffs: Vec<F>) -> Self {
        let mut result = Self { coeffs };
        // While there are zeros at the end of the coefficient vector, pop them off.
        result.truncate_leading_zeros();
        // Check that either the coefficients vec is empty or that the last coeff is
        // non-zero.
        assert!(result.coeffs.last().map_or(true, |coeff| !coeff.is_zero()));

        result
    }

    /// Returns the degree of the polynomial.
    pub fn degree(&self) -> usize {
        if self.is_zero() {
            0
        } else {
            assert!(self.coeffs.last().map_or(false, |coeff| !coeff.is_zero()));
            self.coeffs.len() - 1
        }
    }

    fn truncate_leading_zeros(&mut self) {
        while self.coeffs.last().map_or(false, |c| c.is_zero()) {
            self.coeffs.pop();
        }
    }

    /// Evaluates `self` at the given `point` in the field.
    pub fn evaluate(&self, point: F) -> F {
        if self.is_zero() {
            return F::zero();
        }
        let mut powers_of_point = vec![F::one()];
        let mut cur = point;
        for _ in 0..self.degree() {
            powers_of_point.push(cur);
            cur *= &point;
        }
        assert_eq!(powers_of_point.len(), self.coeffs.len());
        cfg_into_iter!(powers_of_point)
            .zip(&self.coeffs)
            .map(|(power, coeff)| power * coeff)
            .sum()
    }

    /// Perform a naive n^2 multiplication of `self` by `other`.
    pub fn naive_mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            DenseUniPolynomial::zero()
        } else {
            let mut result = vec![F::zero(); self.degree() + other.degree() + 1];
            for (i, self_coeff) in self.coeffs.iter().enumerate() {
                for (j, other_coeff) in other.coeffs.iter().enumerate() {
                    result[i + j] += &(*self_coeff * other_coeff);
                }
            }
            DenseUniPolynomial::from_coefficients_vec(result)
        }
    }

    /// Outputs a polynomial of degree `d` where each coefficient is sampled
    /// uniformly at random from the field `F`.
    pub fn rand<R: Rng>(d: usize, rng: &mut R) -> Self {
        let mut random_coeffs = Vec::new();
        for _ in 0..=d {
            random_coeffs.push(F::rand(rng));
        }
        Self::from_coefficients_vec(random_coeffs)
    }
}

impl<F: FftField> DenseUniPolynomial<F> {
    /// Multiply `self` by the vanishing polynomial for the domain `domain`.
    /// Returns the result of the multiplication.
    pub fn mul_by_vanishing_poly<D: EvaluationDomain<F>>(&self, domain: D) -> DenseUniPolynomial<F> {
        let mut shifted = vec![F::zero(); domain.size()];
        shifted.extend_from_slice(&self.coeffs);
        cfg_iter_mut!(shifted)
            .zip(&self.coeffs)
            .for_each(|(s, c)| *s -= c);
        DenseUniPolynomial::from_coefficients_vec(shifted)
    }

    /// Divide `self` by the vanishing polynomial for the domain `domain`.
    /// Returns the quotient and remainder of the division.
    pub fn divide_by_vanishing_poly<D: EvaluationDomain<F>>(
        &self,
        domain: D,
    ) -> Option<(DenseUniPolynomial<F>, DenseUniPolynomial<F>)> {
        let self_poly: DenseOrSparseUniPolynomial<F> = self.into();
        let vanishing_poly: DenseOrSparseUniPolynomial<F> = domain.vanishing_polynomial().into();
        self_poly.divide_with_q_and_r(&vanishing_poly)
    }
}

impl<'a, 'b, F: Field> Add<&'a DenseUniPolynomial<F>> for &'b DenseUniPolynomial<F> {
    type Output = DenseUniPolynomial<F>;

    fn add(self, other: &'a DenseUniPolynomial<F>) -> DenseUniPolynomial<F> {
        let mut result = if self.is_zero() {
            other.clone()
        } else if other.is_zero() {
            self.clone()
        } else if self.degree() >= other.degree() {
            let mut result = self.clone();
            for (a, b) in result.coeffs.iter_mut().zip(&other.coeffs) {
                *a += b
            }
            result
        } else {
            let mut result = other.clone();
            for (a, b) in result.coeffs.iter_mut().zip(&self.coeffs) {
                *a += b
            }
            result
        };
        result.truncate_leading_zeros();
        result
    }
}

impl<'a, 'b, F: Field> AddAssign<&'a DenseUniPolynomial<F>> for DenseUniPolynomial<F> {
    fn add_assign(&mut self, other: &'a DenseUniPolynomial<F>) {
        if self.is_zero() {
            self.coeffs.truncate(0);
            self.coeffs.extend_from_slice(&other.coeffs);
        } else if other.is_zero() {
        } else if self.degree() >= other.degree() {
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a += b
            }
        } else {
            // Add the necessary number of zero coefficients.
            self.coeffs.resize(other.coeffs.len(), F::zero());
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a += b
            }
            self.truncate_leading_zeros();
        }
    }
}

impl<'a, 'b, F: Field> AddAssign<(F, &'a DenseUniPolynomial<F>)> for DenseUniPolynomial<F> {
    fn add_assign(&mut self, (f, other): (F, &'a DenseUniPolynomial<F>)) {
        if self.is_zero() {
            self.coeffs.truncate(0);
            self.coeffs.extend_from_slice(&other.coeffs);
            self.coeffs.iter_mut().for_each(|c| *c *= &f);
        } else if other.is_zero() {
        } else if self.degree() >= other.degree() {
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a += &(f * b);
            }
        } else {
            // Add the necessary number of zero coefficients.
            self.coeffs.resize(other.coeffs.len(), F::zero());
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a += &(f * b);
            }
            self.truncate_leading_zeros();
        }
    }
}

impl<F: FftField> DenseUniPolynomial<F> {
    /// Evaluate `self` over `domain`.
    pub fn evaluate_over_domain_by_ref<D: EvaluationDomain<F>>(
        &self,
        domain: D,
    ) -> Evaluations<F, D> {
        let poly: DenseOrSparseUniPolynomial<'_, F> = self.into();
        DenseOrSparseUniPolynomial::<F>::evaluate_over_domain(poly, domain)
    }

    /// Evaluate `self` over `domain`.
    pub fn evaluate_over_domain<D: EvaluationDomain<F>>(self, domain: D) -> Evaluations<F, D> {
        let poly: DenseOrSparseUniPolynomial<'_, F> = self.into();
        DenseOrSparseUniPolynomial::<F>::evaluate_over_domain(poly, domain)
    }
}

impl<F: Field> Neg for DenseUniPolynomial<F> {
    type Output = DenseUniPolynomial<F>;

    #[inline]
    fn neg(mut self) -> DenseUniPolynomial<F> {
        for coeff in &mut self.coeffs {
            *coeff = -*coeff;
        }
        self
    }
}

impl<'a, 'b, F: Field> Sub<&'a DenseUniPolynomial<F>> for &'b DenseUniPolynomial<F> {
    type Output = DenseUniPolynomial<F>;

    #[inline]
    fn sub(self, other: &'a DenseUniPolynomial<F>) -> DenseUniPolynomial<F> {
        let mut result = if self.is_zero() {
            let mut result = other.clone();
            for coeff in &mut result.coeffs {
                *coeff = -(*coeff);
            }
            result
        } else if other.is_zero() {
            self.clone()
        } else if self.degree() >= other.degree() {
            let mut result = self.clone();
            for (a, b) in result.coeffs.iter_mut().zip(&other.coeffs) {
                *a -= b
            }
            result
        } else {
            let mut result = self.clone();
            result.coeffs.resize(other.coeffs.len(), F::zero());
            for (a, b) in result.coeffs.iter_mut().zip(&other.coeffs) {
                *a -= b;
            }
            result
        };
        result.truncate_leading_zeros();
        result
    }
}

impl<'a, 'b, F: Field> SubAssign<&'a DenseUniPolynomial<F>> for DenseUniPolynomial<F> {
    #[inline]
    fn sub_assign(&mut self, other: &'a DenseUniPolynomial<F>) {
        if self.is_zero() {
            self.coeffs.resize(other.coeffs.len(), F::zero());
            for (i, coeff) in other.coeffs.iter().enumerate() {
                self.coeffs[i] -= coeff;
            }
        } else if other.is_zero() {
        } else if self.degree() >= other.degree() {
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a -= b
            }
        } else {
            // Add the necessary number of zero coefficients.
            self.coeffs.resize(other.coeffs.len(), F::zero());
            for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                *a -= b
            }
            // If the leading coefficient ends up being zero, pop it off.
            self.truncate_leading_zeros();
        }
    }
}

impl<'a, 'b, F: Field> Div<&'a DenseUniPolynomial<F>> for &'b DenseUniPolynomial<F> {
    type Output = DenseUniPolynomial<F>;

    #[inline]
    fn div(self, divisor: &'a DenseUniPolynomial<F>) -> DenseUniPolynomial<F> {
        let a: DenseOrSparseUniPolynomial<_> = self.into();
        let b: DenseOrSparseUniPolynomial<_> = divisor.into();
        a.divide_with_q_and_r(&b).expect("division failed").0
    }
}

/// Performs O(nlogn) multiplication of polynomials if F is smooth.
impl<'a, 'b, F: FftField> Mul<&'a DenseUniPolynomial<F>> for &'b DenseUniPolynomial<F> {
    type Output = DenseUniPolynomial<F>;

    #[inline]
    fn mul(self, other: &'a DenseUniPolynomial<F>) -> DenseUniPolynomial<F> {
        if self.is_zero() || other.is_zero() {
            DenseUniPolynomial::zero()
        } else {
            let domain = GeneralEvaluationDomain::new(self.coeffs.len() + other.coeffs.len())
                .expect("field is not smooth enough to construct domain");
            let mut self_evals = self.evaluate_over_domain_by_ref(domain);
            let other_evals = other.evaluate_over_domain_by_ref(domain);
            self_evals *= &other_evals;
            self_evals.interpolate()
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::polynomial::*;
    use crate::{EvaluationDomain, GeneralEvaluationDomain};
    use algebra::bls12_381::fr::Fr;
    use algebra_core::{test_rng, Field, One, UniformRand, Zero};

    #[test]
    fn double_polynomials_random() {
        let rng = &mut test_rng();
        for degree in 0..70 {
            let p = DenseUniPolynomial::<Fr>::rand(degree, rng);
            let p_double = &p + &p;
            let p_quad = &p_double + &p_double;
            assert_eq!(&(&(&p + &p) + &p) + &p, p_quad);
        }
    }

    #[test]
    fn add_polynomials() {
        let rng = &mut test_rng();
        for a_degree in 0..70 {
            for b_degree in 0..70 {
                let p1 = DenseUniPolynomial::<Fr>::rand(a_degree, rng);
                let p2 = DenseUniPolynomial::<Fr>::rand(b_degree, rng);
                let res1 = &p1 + &p2;
                let res2 = &p2 + &p1;
                assert_eq!(res1, res2);
            }
        }
    }

    #[test]
    fn add_polynomials_with_mul() {
        let rng = &mut test_rng();
        for a_degree in 0..70 {
            for b_degree in 0..70 {
                let mut p1 = DenseUniPolynomial::rand(a_degree, rng);
                let p2 = DenseUniPolynomial::rand(b_degree, rng);
                let f = Fr::rand(rng);
                let f_p2 = DenseUniPolynomial::from_coefficients_vec(
                    p2.coeffs.iter().map(|c| f * c).collect(),
                );
                let res2 = &f_p2 + &p1;
                p1 += (f, &p2);
                let res1 = p1;
                assert_eq!(res1, res2);
            }
        }
    }

    #[test]
    fn sub_polynomials() {
        let rng = &mut test_rng();
        let p1 = DenseUniPolynomial::<Fr>::rand(5, rng);
        let p2 = DenseUniPolynomial::<Fr>::rand(3, rng);
        let res1 = &p1 - &p2;
        let res2 = &p2 - &p1;
        assert_eq!(
            &res1 + &p2,
            p1,
            "Subtraction should be inverse of addition!"
        );
        assert_eq!(res1, -res2, "p2 - p1 = -(p1 - p2)");
    }

    #[test]
    fn divide_polynomials_fixed() {
        let dividend = DenseUniPolynomial::from_coefficients_slice(&[
            "4".parse().unwrap(),
            "8".parse().unwrap(),
            "5".parse().unwrap(),
            "1".parse().unwrap(),
        ]);
        let divisor = DenseUniPolynomial::from_coefficients_slice(&[Fr::one(), Fr::one()]); // Construct a monic linear polynomial.
        let result = &dividend / &divisor;
        let expected_result = DenseUniPolynomial::from_coefficients_slice(&[
            "4".parse().unwrap(),
            "4".parse().unwrap(),
            "1".parse().unwrap(),
        ]);
        assert_eq!(expected_result, result);
    }

    #[test]
    fn divide_polynomials_random() {
        let rng = &mut test_rng();

        for a_degree in 0..70 {
            for b_degree in 0..70 {
                let dividend = DenseUniPolynomial::<Fr>::rand(a_degree, rng);
                let divisor = DenseUniPolynomial::<Fr>::rand(b_degree, rng);
                if let Some((quotient, remainder)) = DenseOrSparseUniPolynomial::divide_with_q_and_r(
                    &(&dividend).into(),
                    &(&divisor).into(),
                ) {
                    assert_eq!(dividend, &(&divisor * &quotient) + &remainder)
                }
            }
        }
    }

    #[test]
    fn evaluate_polynomials() {
        let rng = &mut test_rng();
        for a_degree in 0..70 {
            let p = DenseUniPolynomial::rand(a_degree, rng);
            let point: Fr = Fr::from(10u64);
            let mut total = Fr::zero();
            for (i, coeff) in p.coeffs.iter().enumerate() {
                total += &(point.pow(&[i as u64]) * coeff);
            }
            assert_eq!(p.evaluate(point), total);
        }
    }

    #[test]
    fn mul_polynomials_random() {
        let rng = &mut test_rng();
        for a_degree in 0..70 {
            for b_degree in 0..70 {
                let a = DenseUniPolynomial::<Fr>::rand(a_degree, rng);
                let b = DenseUniPolynomial::<Fr>::rand(b_degree, rng);
                assert_eq!(&a * &b, a.naive_mul(&b))
            }
        }
    }

    #[test]
    fn mul_by_vanishing_poly() {
        let rng = &mut test_rng();
        for size in 1..10 {
            let domain = GeneralEvaluationDomain::new(1 << size).unwrap();
            for degree in 0..70 {
                let p = DenseUniPolynomial::<Fr>::rand(degree, rng);
                let ans1 = p.mul_by_vanishing_poly(domain);
                let ans2 = &p * &domain.vanishing_polynomial().into();
                assert_eq!(ans1, ans2);
            }
        }
    }

    #[test]
    fn test_leading_zero() {
        let n = 10;
        let rand_poly = DenseUniPolynomial::rand(n, &mut test_rng());
        let coefficients = rand_poly.coeffs.clone();
        let leading_coefficient: Fr = coefficients[n];

        let negative_leading_coefficient = -leading_coefficient;
        let inverse_leading_coefficient = leading_coefficient.inverse().unwrap();

        let mut inverse_coefficients = coefficients.clone();
        inverse_coefficients[n] = inverse_leading_coefficient;

        let mut negative_coefficients = coefficients;
        negative_coefficients[n] = negative_leading_coefficient;

        let negative_poly = DenseUniPolynomial::from_coefficients_vec(negative_coefficients);
        let inverse_poly = DenseUniPolynomial::from_coefficients_vec(inverse_coefficients);

        let x = &inverse_poly * &rand_poly;
        assert_eq!(x.degree(), 2 * n);
        assert!(!x.coeffs.last().unwrap().is_zero());

        let y = &negative_poly + &rand_poly;
        assert_eq!(y.degree(), n - 1);
        assert!(!y.coeffs.last().unwrap().is_zero());
    }
}
