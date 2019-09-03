use std::fmt;
use std::ops::{Add, AddAssign, Deref, DerefMut, Div, Mul, Neg, Sub, SubAssign};

use crate::{Field, PrimeField};
use crate::fft::domain::EvaluationDomain;
use crate::fft::evaluations::Evaluations;
use rand::Rng;
use rayon::prelude::*;

/// Stores a polynomial in coefficient form.
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct Polynomial<F: Field> {
    /// The coefficient of `x^i` is stored at location `i` in the `coeffs` vector.
    pub coeffs: Vec<F>,
}

impl<F: Field> fmt::Debug for Polynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        for (i, coeff) in self.coeffs.iter().enumerate() {
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

impl<F: Field> Deref for Polynomial<F> {
    type Target = [F];

    fn deref(&self) -> &[F] {
        &self.coeffs
    }
}

impl<F: Field> DerefMut for Polynomial<F> {
    fn deref_mut(&mut self) -> &mut [F] {
        &mut self.coeffs
    }
}

impl<F: Field> Polynomial<F> {
    /// Returns the zero polynomial.
    pub fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    /// Checks if the given polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.len() == 0 || self.coeffs.iter().all(|coeff| coeff.is_zero())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_slice(coeffs: &[F]) -> Self {
        Self::from_coefficients_vec(coeffs.to_vec())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_vec(mut coeffs: Vec<F>) -> Self {
        // While there are zeros at the end of the coefficient vector, pop them off.
        while coeffs.last().map_or(false, |c| c.is_zero()) {
            coeffs.pop();
        }
        // Check that either the coefficients vec is empty or that the last coeff is non-zero.
        //
        assert!(coeffs.last().map_or(true, |coeff| !coeff.is_zero()));

        Self { coeffs }
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
        let zero = F::zero();
        powers_of_point
            .into_par_iter()
            .zip(&self.coeffs)
            .map(|(power, coeff)| power * coeff)
            .reduce(|| zero, |a, b| a + &b)
    }

    pub fn naive_mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            Polynomial::zero()
        } else {
            let mut result = vec![F::zero(); self.degree() + other.degree() + 1];
            for (i, self_coeff) in self.coeffs.iter().enumerate() {
                for (j, other_coeff) in other.coeffs.iter().enumerate() {
                    result[i + j] += &(*self_coeff * other_coeff);
                }
            }
            Polynomial::from_coefficients_vec(result)
        }
    }

    /// Outputs a polynomial of degree `d` where each coefficient is sampled uniformly at random
    /// from the field `F`.
    pub fn rand<R: Rng>(d: usize, rng: &mut R) -> Self {
        let mut random_coeffs = Vec::new();
        for _ in 0..(d + 1) {
            random_coeffs.push(F::rand(rng));
        }
        Self::from_coefficients_vec(random_coeffs)
    }
}

impl<'a, 'b, F: Field> Add<&'a Polynomial<F>> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    fn add(self, other: &'a Polynomial<F>) -> Polynomial<F> {
        if self.is_zero() {
            other.clone()
        } else if other.is_zero() {
            self.clone()
        } else {
            if self.degree() >= other.degree() {
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
                // If the leading coefficient ends up being zero, pop it off.
                while result.coeffs.last().unwrap().is_zero() {
                    result.coeffs.pop();
                }
                result
            }
        }
    }
}

impl<'a, 'b, F: Field> AddAssign<&'a Polynomial<F>> for Polynomial<F> {
    fn add_assign(&mut self, other: &'a Polynomial<F>) {
        if self.is_zero() {
            self.coeffs.truncate(0);
            self.coeffs.extend_from_slice(&other.coeffs);
        } else if other.is_zero() {
            return;
        } else {
            if self.degree() >= other.degree() {
                for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                    *a += b
                }
            } else {
                // Add the necessary number of zero coefficients.
                self.coeffs.resize(other.coeffs.len(), F::zero());
                for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                    *a += b
                }
                // If the leading coefficient ends up being zero, pop it off.
                while self.coeffs.last().unwrap().is_zero() {
                    self.coeffs.pop();
                }
            }
        }
    }
}

impl<'a, 'b, F: Field> AddAssign<(F, &'a Polynomial<F>)> for Polynomial<F> {
    fn add_assign(&mut self, (f, other): (F, &'a Polynomial<F>)) {
        if self.is_zero() {
            self.coeffs.truncate(0);
            self.coeffs.extend_from_slice(&other.coeffs);
            self.coeffs.iter_mut().for_each(|c| *c *= &f);
        } else if other.is_zero() {
            return;
        } else {
            if self.degree() >= other.degree() {
                for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                    *a += &(f * b);
                }
            } else {
                // Add the necessary number of zero coefficients.
                self.coeffs.resize(other.coeffs.len(), F::zero());
                for (a, b) in self.coeffs.iter_mut().zip(&other.coeffs) {
                    *a += &(f * b);
                }
                // If the leading coefficient ends up being zero, pop it off.
                while self.coeffs.last().unwrap().is_zero() {
                    self.coeffs.pop();
                }
            }
        }
    }
}

impl<F: Field> Neg for Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn neg(mut self) -> Polynomial<F> {
        for coeff in &mut self.coeffs {
            *coeff = -*coeff;
        }
        self
    }
}

impl<'a, 'b, F: Field> Sub<&'a Polynomial<F>> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn sub(self, other: &'a Polynomial<F>) -> Polynomial<F> {
        if self.is_zero() {
            let mut result = other.clone();
            for coeff in &mut result.coeffs {
                *coeff = -(*coeff);
            }
            result
        } else if other.is_zero() {
            self.clone()
        } else {
            if self.degree() >= other.degree() {
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
                if !result.is_zero() {
                    // If the leading coefficient ends up being zero, pop it off.
                    while result.coeffs.last().unwrap().is_zero() {
                        result.coeffs.pop();
                    }
                }

                result
            }
        }
    }
}

impl<'a, 'b, F: Field> SubAssign<&'a Polynomial<F>> for Polynomial<F> {
    #[inline]
    fn sub_assign(&mut self, other: &'a Polynomial<F>) {
        if self.is_zero() {
            self.coeffs.resize(other.coeffs.len(), F::zero());
            for (i, coeff) in other.coeffs.iter().enumerate() {
                self.coeffs[i] -= coeff;
            }
        } else if other.is_zero() {
            return;
        } else {
            if self.degree() > other.degree() {
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
                while self.coeffs.last().unwrap().is_zero() {
                    self.coeffs.pop();
                }
            }
        }
    }
}

impl<'a, 'b, F: Field> Div<&'a Polynomial<F>> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn div(self, divisor: &'a Polynomial<F>) -> Polynomial<F> {
        if self.is_zero() {
            Polynomial::zero()
        } else if divisor.is_zero() {
            panic!("Dividing by zero polynomial")
        } else {
            if self.degree() < divisor.degree() {
                Polynomial::zero()
            } else {
                let mut quotient = vec![F::zero(); self.degree() - divisor.degree() + 1];
                let mut remainder = self.clone();
                // Can unwrap here because we know it is not zero.
                let divisor_leading = divisor.last().unwrap();
                let divisor_is_monic = divisor_leading.is_one();
                while !remainder.is_zero() && remainder.degree() >= divisor.degree() {
                    let cur_q_coeff = if divisor_is_monic {
                        *remainder.last().unwrap()
                    } else {
                        *remainder.last().unwrap() / divisor_leading
                    };
                    let cur_q_degree = remainder.degree() - divisor.degree();
                    quotient[cur_q_degree] = cur_q_coeff;

                    for (i, coeff) in divisor.coeffs.iter().enumerate() {
                        remainder[cur_q_degree + i] -= &(cur_q_coeff * coeff);
                    }
                    while remainder
                        .coeffs
                        .last()
                        .map_or(false, |coeff| coeff.is_zero())
                    {
                        remainder.coeffs.pop();
                    }
                }
                Polynomial::from_coefficients_vec(quotient)
            }
        }
    }
}

/// Performs O(nlogn) multiplication of polynomials if F is smooth.
impl<'a, 'b, F: PrimeField> Mul<&'a Polynomial<F>> for &'b Polynomial<F> {
    type Output = Polynomial<F>;

    #[inline]
    fn mul(self, other: &'a Polynomial<F>) -> Polynomial<F> {
        if self.is_zero() || other.is_zero() {
            Polynomial::zero()
        } else {
            let domain = EvaluationDomain::new(self.coeffs.len() + other.coeffs.len()).expect("field is not smooth enough to construct domain");
            let mut self_evals = Evaluations::evaluate_polynomial_over_domain(self, domain);
            let other_evals = Evaluations::evaluate_polynomial_over_domain(other, domain);
            self_evals *= &other_evals;
            self_evals.interpolate_over_domain(domain)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::fft::polynomial::*;
    use crate::fields::{bls12_381::fr::Fr, Field};
    use rand::{thread_rng, Rand};

    #[test]
    fn double_polynomials_random() {
        let rng = &mut thread_rng();
        let p = Polynomial::<Fr>::rand(5, rng);
        let p_double = &p + &p;
        let p_quad = &p_double + &p_double;
        assert_eq!(&(&(&p + &p) + &p) + &p, p_quad);
    }

    #[test]
    fn add_polynomials() {
        let rng = &mut thread_rng();
        let p1 = Polynomial::<Fr>::rand(5, rng);
        let p2 = Polynomial::<Fr>::rand(4, rng);;
        let res1 = &p1 + &p2;
        let res2 = &p2 + &p1;
        assert_eq!(res1, res2);
    }

    #[test]
    fn add_polynomials_with_mul() {
        let rng = &mut thread_rng();
        let mut p1 = Polynomial::rand(5, rng);
        let p2 = Polynomial::rand(4, rng);;
        let f = Fr::rand(rng);

        let f_p2 = Polynomial::from_coefficients_vec(p2.coeffs.iter().map(|c| f * c).collect());

        let res2 = &f_p2 + &p1;
        p1 += (f, &p2);
        let res1 = p1;
        assert_eq!(res1, res2);
    }

    #[test]
    fn sub_polynomials() {
        let rng = &mut thread_rng();
        let p1 = Polynomial::<Fr>::rand(5, rng);
        let p2 = Polynomial::<Fr>::rand(3, rng);
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
        let dividend = Polynomial::from_coefficients_slice(&[
            "4".parse().unwrap(),
            "8".parse().unwrap(),
            "5".parse().unwrap(),
            "1".parse().unwrap(),
        ]);
        let divisor = Polynomial::from_coefficients_slice(&[Fr::one(), Fr::one()]); // Construct a monic linear polynomial.
        let result = &dividend / &divisor;
        let expected_result = Polynomial::from_coefficients_slice(&[
            "4".parse().unwrap(),
            "4".parse().unwrap(),
            "1".parse().unwrap(),
        ]);
        assert_eq!(expected_result, result);
    }

    #[test]
    fn divide_polynomials_random() {
        let rng = &mut thread_rng();

        let factor1 = Polynomial::<Fr>::rand(5, rng);
        let factor2 = Polynomial::<Fr>::rand(5, rng);
        let prod = &factor1 * &factor2;
        let result1 = &prod / &factor1;
        let result2 = &prod / &factor2;
        assert_eq!(factor2, result1);
        assert_eq!(factor1, result2);
    }

    #[test]
    fn evaluate_polynomials() {
        let rng = &mut thread_rng();
        let p = Polynomial::rand(5, rng);
        let point: Fr = "10".parse().unwrap();
        let mut total = Fr::zero();
        for (i, coeff) in p.coeffs.iter().enumerate() {
            total += &(point.pow(&[i as u64]) * coeff);
        }
        assert_eq!(p.evaluate(point), total);
    }

    #[test]
    fn mul_polynomials_random() {
        let rng = &mut thread_rng();
        for a_degree in 0..70 {
            for b_degree in 0..70 {
                let a = Polynomial::<Fr>::rand(a_degree, rng);
                let b = Polynomial::<Fr>::rand(b_degree, rng);
                assert_eq!(&a * &b, a.naive_mul(&b))
            }
        }
    }
}
