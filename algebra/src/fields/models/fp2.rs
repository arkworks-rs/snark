use super::quadratic_extension::*;
use crate::fields::{PrimeField, SquareRootField};
use std::marker::PhantomData;

pub trait Fp2Parameters: 'static + Send + Sync {
    type Fp: PrimeField + SquareRootField;

    //alpha
    const NONRESIDUE: Self::Fp;
    //quadratic nonresidue for square root algorithm
    const QUADRATIC_NONRESIDUE: (Self::Fp, Self::Fp);
    //coefficients of the powers of the Frobenius automorphism as linear map over F
    // (pi^0(X), pi^1(X)) = (C1_0*X, C1_1*X),
    const FROBENIUS_COEFF_FP2_C1: &'static [Self::Fp];

    #[inline(always)]
    fn mul_fp_by_nonresidue(fe: &Self::Fp) -> Self::Fp {
        Self::NONRESIDUE * fe
    }
}

pub struct Fp2ParamsWrapper<P: Fp2Parameters>(PhantomData<P>);

impl<P: Fp2Parameters> QuadExtParameters for Fp2ParamsWrapper<P> {
    type BasePrimeField = P::Fp;
    type BaseField = P::Fp;
    type FrobCoeff = P::Fp;

    const DEGREE_OVER_BASE_PRIME_FIELD: usize = 2;

    const NONRESIDUE: Self::BaseField = P::NONRESIDUE;

    const FROBENIUS_COEFF_C1: &'static [Self::FrobCoeff] = P::FROBENIUS_COEFF_FP2_C1;

    #[inline(always)]
    fn mul_base_field_by_nonresidue(fe: &Self::BaseField) -> Self::BaseField {
        P::mul_fp_by_nonresidue(fe)
    }

    fn mul_base_field_by_frob_coeff(fe: &mut Self::BaseField, power: usize) {
        *fe *= &Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}

pub type Fp2<P> = QuadExtField<Fp2ParamsWrapper<P>>;

impl<P: Fp2Parameters> Fp2<P> {
    pub fn mul_assign_by_fp(&mut self, other: &P::Fp) {
        self.c0 *= other;
        self.c1 *= other;
    }
}
