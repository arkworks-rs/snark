//! A polynomial represented in evaluations form.

use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};
use crate::{PrimeField, ToBytes, FromBytes};
use crate::{DensePolynomial, EvaluationDomain, get_best_evaluation_domain};

/// Stores a polynomial in evaluation form.
pub struct Evaluations<F: PrimeField> {
    /// The evaluations of a polynomial over the domain `D`
    pub evals: Vec<F>,
    /// Evaluation domain
    pub domain: Box<dyn EvaluationDomain<F>>,
}

impl<F: PrimeField> Evaluations<F> {
    /// Construct `Self` from evaluations and a domain.
    pub fn from_vec_and_domain(evals: Vec<F>, domain: Box<dyn EvaluationDomain<F>>) -> Self {
        Self {
            evals,
            domain,
        }
    }

    /// Interpolate a polynomial from a list of evaluations
    pub fn interpolate_by_ref(&self) -> DensePolynomial<F> {
        DensePolynomial::from_coefficients_vec(self.domain.ifft(&self.evals))
    }

    /// Interpolate a polynomial from a list of evaluations
    pub fn interpolate(self) -> DensePolynomial<F> {
        let Self { mut evals, domain } = self;
        domain.ifft_in_place(&mut evals);
        DensePolynomial::from_coefficients_vec(evals)
    }
}

impl<F: PrimeField> ToBytes for Evaluations<F>
{
    fn write<W: std::io::Write>(&self, mut w: W) -> std::io::Result<()> {
        (self.evals.len() as u64).write(&mut w)?;
        for e in self.evals.iter() {
            e.write(&mut w)?;
        }
        Ok(())
    }
}

impl<F: PrimeField> FromBytes for Evaluations<F>
{
    #[inline]
    fn read<Read: std::io::Read>(mut reader: Read) -> std::io::Result<Evaluations<F>> {
        let evals_count = u64::read(&mut reader)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
        let mut evals = vec![];
        for _ in 0..evals_count {
            let eval = F::read(&mut reader)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            evals.push(eval);
        }
        let domain = get_best_evaluation_domain::<F>(evals.len()).unwrap();
        Ok(Evaluations::<F> {
            evals,
            domain
        })
    }
}

impl<F: PrimeField> std::ops::Index<usize> for Evaluations<F> {
    type Output = F;

    fn index(&self, index: usize) -> &F {
        &self.evals[index]
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
        assert_eq!(self.domain.as_ref(), other.domain.as_ref(), "domains are unequal");
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
        assert_eq!(self.domain.as_ref(), other.domain.as_ref(), "domains are unequal");
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
        assert_eq!(self.domain.as_ref(), other.domain.as_ref(), "domains are unequal");
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
        assert_eq!(self.domain.as_ref(), other.domain.as_ref(), "domains are unequal");
        self.evals.iter_mut().zip(&other.evals).for_each(|(a, b)| *a /= b);
    }
}

impl<F: PrimeField> Clone for Evaluations<F> {
    fn clone(&self) -> Self {
        Self{
            evals: self.evals.clone(),
            domain: self.domain.clone(),
        }
    }
}

impl<F: PrimeField> PartialEq for Evaluations<F>{
    fn eq(&self, other: &Self) -> bool {
        self.evals.eq(&other.evals) &&
            self.domain.eq(&other.domain)
    }
}
