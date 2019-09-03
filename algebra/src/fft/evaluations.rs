use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};
use crate::PrimeField;
use crate::fft::{polynomial::Polynomial, domain::EvaluationDomain};

/// Stores a polynomial in coefficient form.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Evaluations<F: PrimeField> {
    /// The evaluations of a polynomial over the domain `D`
    pub evals: Vec<F>,
    #[doc(hidden)]
    domain: EvaluationDomain<F>,
}


impl<F: PrimeField> Evaluations<F> {
    pub fn from_vec_and_domain(evals: Vec<F>, domain: EvaluationDomain<F>) -> Self {
        Self {
            evals,
            domain,
        }
    }

    pub fn evaluate_polynomial_over_domain(poly: &Polynomial<F>, domain: EvaluationDomain<F>) -> Self {
        Self::from_vec_and_domain(domain.fft(&poly.coeffs), domain)
    }

    pub fn evaluate_polynomial_over_domain_in_place(mut poly: Polynomial<F>, domain: EvaluationDomain<F>) -> Self {
        domain.fft_in_place(&mut poly.coeffs);
        Self::from_vec_and_domain(poly.coeffs, domain)
    }

    /// Interpolate a polynomial from a list of evaluations
    pub fn interpolate_over_domain(&self, domain: EvaluationDomain<F>) -> Polynomial<F> {
        Polynomial::from_coefficients_vec(domain.ifft(&self.evals))
    }
}

impl<'a, 'b, F: PrimeField> Mul<&'a Evaluations<F>> for &'b Evaluations<F> {

    type Output = Evaluations<F>;

    #[inline]
    fn mul(self, other: &'a Evaluations<F>) -> Evaluations<F> {
        let mut result = self.clone();
        result *= other;
        result
    }
}

impl<'a, F: PrimeField> MulAssign<&'a Evaluations<F>> for Evaluations<F> {
    #[inline]
    fn mul_assign(&mut self, other: &'a Evaluations<F>) {
        assert_eq!(self.domain, other.domain, "domains are unequal");
        self.evals.iter_mut().zip(&other.evals).for_each(|(a, b)| *a *= b);
    }
}

impl<'a, 'b, F: PrimeField> Add<&'a Evaluations<F>> for &'b Evaluations<F> {

    type Output = Evaluations<F>;

    #[inline]
    fn add(self, other: &'a Evaluations<F>) -> Evaluations<F> {
        let mut result = self.clone();
        result += other;
        result
    }
}

impl<'a, F: PrimeField> AddAssign<&'a Evaluations<F>> for Evaluations<F> {
    #[inline]
    fn add_assign(&mut self, other: &'a Evaluations<F>) {
        assert_eq!(self.domain, other.domain, "domains are unequal");
        self.evals.iter_mut().zip(&other.evals).for_each(|(a, b)| *a += b);
    }
}

impl<'a, 'b, F: PrimeField> Sub<&'a Evaluations<F>> for &'b Evaluations<F> {

    type Output = Evaluations<F>;

    #[inline]
    fn sub(self, other: &'a Evaluations<F>) -> Evaluations<F> {
        let mut result = self.clone();
        result -= other;
        result
    }
}

impl<'a, F: PrimeField> SubAssign<&'a Evaluations<F>> for Evaluations<F> {
    #[inline]
    fn sub_assign(&mut self, other: &'a Evaluations<F>) {
        assert_eq!(self.domain, other.domain, "domains are unequal");
        self.evals.iter_mut().zip(&other.evals).for_each(|(a, b)| *a -= b);
    }
}

impl<'a, 'b, F: PrimeField> Div<&'a Evaluations<F>> for &'b Evaluations<F> {

    type Output = Evaluations<F>;

    #[inline]
    fn div(self, other: &'a Evaluations<F>) -> Evaluations<F> {
        let mut result = self.clone();
        result /= other;
        result
    }
}

impl<'a, F: PrimeField> DivAssign<&'a Evaluations<F>> for Evaluations<F> {
    #[inline]
    fn div_assign(&mut self, other: &'a Evaluations<F>) {
        assert_eq!(self.domain, other.domain, "domains are unequal");
        self.evals.iter_mut().zip(&other.evals).for_each(|(a, b)| *a /= b);
    }
}
