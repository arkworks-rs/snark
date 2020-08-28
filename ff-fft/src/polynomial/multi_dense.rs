//! A dense multivariate polynomial represented in coefficient form.
use algebra_core::Field;
use core::{
    fmt,
    ops::{Add, AddAssign, Neg, Sub, SubAssign},
};
use ndarray::{Array, Dimension, IxDyn};
use rand::Rng;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Stores a dense multi-variate polynomial in coefficient form.
#[derive(Clone, PartialEq, Eq, Hash, Default)]
pub struct DenseMultiPolynomial<F: Field> {
    coeffs: Array<F, IxDyn>,
}

impl<F: Field> DenseMultiPolynomial<F> {
    /// Convert variable degree bounds to coefficient matrix dimensions
    ///
    /// TODO: Add one
    fn to_dim(deg_bounds: &[usize]) -> IxDyn {
        IxDyn(&deg_bounds.iter().map(|e| e + 1).collect::<Vec<_>>())
    }

    /// Return the number of variables
    pub fn num_vars(&self) -> usize {
        self.coeffs.dim().ndim()
    }

    /// Returns the zero polynomial.
    pub fn zero(deg_bounds: &[usize]) -> Self {
        Self {
            coeffs: Array::zeros(DenseMultiPolynomial::<F>::to_dim(deg_bounds)),
        }
    }

    /// Checks if the given polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty() || self.coeffs.iter().all(|c| c.is_zero())
    }

    /// Constructs a new polynomial from a list of coefficients.
    pub fn from_coefficients_slice(
        deg_bounds: &[usize],
        coeffs: &[(Vec<(usize, usize)>, F)],
    ) -> Self {
        let mut poly = DenseMultiPolynomial::<F>::zero(deg_bounds);
        for (term, coeff) in coeffs {
            let mut idx = vec![0; poly.num_vars()];
            term.iter().for_each(|(var, exp)| idx[*var] = *exp);
            poly.coeffs[IxDyn(&idx)] += *coeff;
        }
        poly
    }

    /// Evaluates `self` at the given `point` in the field.
    pub fn evaluate(&self, point: &[F]) -> F {
        assert!(point.len() >= self.num_vars(), "Invalid evaluation domain");
        if self.is_zero() {
            return F::zero();
        }
        cfg_into_bridged_iter!(self.coeffs.indexed_iter())
            .map(|(idx, coeff)| {
                let term = idx
                    .slice()
                    .iter()
                    .enumerate()
                    .map(|(i, exp)| {
                        let term = point.get(i).expect("Invalid evaluation point");
                        term.pow(&[*exp as u64])
                    })
                    .product::<F>();
                *coeff * term
            })
            .sum()
    }

    /// Outputs a polynomial which is the sum of `l` `d`-degree univariate
    /// polynomials where each coefficient is sampled uniformly at random
    /// from the field `F`.
    pub fn rand<R: Rng>(l: usize, d: usize, rng: &mut R) -> Self {
        let mut random_coeffs = Vec::new();
        for var in 0..l {
            for deg in 0..=d {
                random_coeffs.push((vec![(var, deg)], F::rand(rng)));
            }
        }
        Self::from_coefficients_slice(&vec![d; l], &random_coeffs)
    }
}

impl<'a, 'b, F: Field> Add<&'a DenseMultiPolynomial<F>> for &'b DenseMultiPolynomial<F> {
    type Output = DenseMultiPolynomial<F>;

    #[inline]
    fn add(self, other: &'a DenseMultiPolynomial<F>) -> DenseMultiPolynomial<F> {
        assert_eq!(self.coeffs.dim(), other.coeffs.dim(), "Invalid dimensions");
        DenseMultiPolynomial {
            coeffs: &self.coeffs + &other.coeffs,
        }
    }
}

impl<'a, F: Field> AddAssign<&'a DenseMultiPolynomial<F>> for DenseMultiPolynomial<F> {
    
    #[inline]
    fn add_assign(&mut self, other: &'a DenseMultiPolynomial<F>) {
        assert_eq!(self.coeffs.dim(), other.coeffs.dim(), "Invalid dimensions");

        cfg_into_bridged_iter!(self.coeffs.iter_mut().zip(other.coeffs.iter()))
            .for_each(|(e, o)| *e += o);
    }
}

impl<F: Field> Neg for DenseMultiPolynomial<F> {
    type Output = DenseMultiPolynomial<F>;

    #[inline]
    fn neg(mut self) -> DenseMultiPolynomial<F> {
        self.coeffs = -self.coeffs;
        self
    }
}

impl<'a, 'b, F: Field> Sub<&'a DenseMultiPolynomial<F>> for &'b DenseMultiPolynomial<F> {
    type Output = DenseMultiPolynomial<F>;

    #[inline]
    fn sub(self, other: &'a DenseMultiPolynomial<F>) -> DenseMultiPolynomial<F> {
        let neg_other = other.clone().neg();
        self + &neg_other
    }
}

impl<'a, F: Field> SubAssign<&'a DenseMultiPolynomial<F>> for DenseMultiPolynomial<F> {
    
    #[inline]
    fn sub_assign(&mut self, other: &'a DenseMultiPolynomial<F>) {
        let neg_other = other.clone().neg();
        *self += &neg_other
    }
}

impl<F: Field> fmt::Debug for DenseMultiPolynomial<F> {
    fn fmt(&self, f: &mut fmt::Formatter) -> Result<(), fmt::Error> {
        for (idx, coeff) in self.coeffs.indexed_iter().filter(|(_, c)| !c.is_zero()) {
            write!(f, "\n{:?}", coeff)?;
            for (var, exp) in idx.slice().iter().enumerate() {
                if *exp == 1 {
                    write!(f, " * x_{}", var)?;
                } else if *exp != 0 {
                    write!(f, " * x_{}^{}", var, exp)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebra::bls12_381::fr::Fr;
    use algebra_core::{test_rng, Field, UniformRand, Zero};

    /// Generate random `l`-variate polynomial of maximum individual degree `d`
    fn rand_poly<R: Rng>(deg_bounds: &[usize], rng: &mut R) -> DenseMultiPolynomial<Fr> {
        let mut random_coeffs = Vec::new();
        let num_terms = rng.gen_range(1, 1000);
        // For each term, randomly select up to `l` variables with degree
        // in [1,d] and random coefficient
        let l = deg_bounds.len();
        for _ in 0..num_terms {
            let term = (0..l)
                .map(|i| {
                    if rng.gen_bool(0.5) {
                        Some((i, rng.gen_range(1, deg_bounds[i] + 1)))
                    } else {
                        None
                    }
                })
                .filter(|t| t.is_some())
                .map(|t| t.unwrap())
                .collect();
            let coeff = Fr::rand(rng);
            random_coeffs.push((term, coeff));
        }
        DenseMultiPolynomial::from_coefficients_slice(deg_bounds, &random_coeffs)
    }

    #[test]
    fn add_polynomials() {
        let rng = &mut test_rng();
        let max_degree = 5;
        for var_count in 0..10 {
            let mut deg_bounds = Vec::with_capacity(var_count);
            for _ in 0..var_count {
                deg_bounds.push(rng.gen_range(1, max_degree + 1));
            }
            let p1 = rand_poly(&deg_bounds, rng);
            let p2 = rand_poly(&deg_bounds, rng);
            // Test Add
            let res1 = &p1 + &p2;
            let res2 = &p2 + &p1;
            assert_eq!(res1, res2);
            // Test AddAssign
            let mut cloned_p1 = p1.clone();
            let mut cloned_p2 = p2.clone();
            cloned_p1 += &p2;
            cloned_p2 += &p1;
            assert_eq!(cloned_p1, cloned_p2);
            assert_eq!(cloned_p1, res1);
        }
    }

    #[test]
    fn sub_polynomials() {
        let rng = &mut test_rng();
        let max_degree = 5;
        for var_count in 0..10 {
            let mut deg_bounds = Vec::with_capacity(var_count);
            for _ in 0..var_count {
                deg_bounds.push(rng.gen_range(1, max_degree + 1));
            }
            let p1 = rand_poly(&deg_bounds, rng);
            let p2 = rand_poly(&deg_bounds, rng);
            // Test Sub
            let res1 = &p1 - &p2;
            let res2 = &p2 - &p1;
            assert_eq!(&res1 + &p2, p1);
            assert_eq!(res1, -res2.clone());
            // Test SubAssign
            let mut cloned_p1 = p1.clone();
            let mut cloned_p2 = p2.clone();
            cloned_p1 -= &p2;
            cloned_p2 -= &p1;
            assert_eq!(cloned_p1, res1);
            assert_eq!(cloned_p2, res2);
        }
    }

    #[test]
    fn evaluate_polynomials() {
        let rng = &mut test_rng();
        let max_degree = 5;
        for var_count in 0..10 {
            let mut deg_bounds = Vec::with_capacity(var_count);
            for _ in 0..var_count {
                deg_bounds.push(rng.gen_range(1, max_degree + 1));
            }
            let p = rand_poly(&deg_bounds, rng);
            let mut point = Vec::with_capacity(var_count);
            for _ in 0..var_count {
                point.push(Fr::rand(rng));
            }
            let mut total = Fr::zero();
            for (idx, coeff) in p.coeffs.indexed_iter() {
                let mut summand = *coeff;
                for (var, exp) in idx.slice().iter().enumerate() {
                    let eval = point.get(var).unwrap();
                    summand *= eval.pow(&[*exp as u64]);
                }
                total += summand;
            }
            assert_eq!(p.evaluate(&point), total);
        }
    }

    #[test]
    fn add_and_evaluate_polynomials() {
        let rng = &mut test_rng();
        let max_degree = 5;
        for var_count in 0..10 {
            let mut deg_bounds = Vec::with_capacity(var_count);
            for _ in 0..var_count {
                deg_bounds.push(rng.gen_range(1, max_degree + 1));
            }

            let p1 = rand_poly(&deg_bounds, rng);
            let p2 = rand_poly(&deg_bounds, rng);
            let mut point = Vec::new();
            for _ in 0..var_count {
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
