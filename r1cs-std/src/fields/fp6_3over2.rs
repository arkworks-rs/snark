use crate::fields::{cubic_extension::*, fp2::*};
use algebra::fields::{fp6_3over2::*, CubicExtParameters, Fp2};
use core::ops::MulAssign;
use r1cs_core::SynthesisError;

/// A sextic extension field constructed as the tower of a
/// cubic extension over a quadratic extension field.
/// This is the R1CS equivalent of `algebra_core::fp6_3over3::Fp6<P>`.
pub type Fp6Var<P> = CubicExtVar<Fp2Var<<P as Fp6Parameters>::Fp2Params>, Fp6ParamsWrapper<P>>;

impl<P: Fp6Parameters> CubicExtVarParams<Fp2Var<P::Fp2Params>> for Fp6ParamsWrapper<P> {
    fn mul_base_field_vars_by_frob_coeff(
        c1: &mut Fp2Var<P::Fp2Params>,
        c2: &mut Fp2Var<P::Fp2Params>,
        power: usize,
    ) {
        *c1 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        *c2 *= Self::FROBENIUS_COEFF_C2[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}

impl<P: Fp6Parameters> Fp6Var<P> {
    /// Multiplies `self` by a sparse element which has `c0 == c2 == zero`.
    pub fn mul_by_0_c1_0(&self, c1: &Fp2Var<P::Fp2Params>) -> Result<Self, SynthesisError> {
        // Karatsuba multiplication
        // v0 = a0 * b0 = 0

        // v1 = a1 * b1
        let v1 = &self.c1 * c1;

        // v2 = a2 * b2 = 0

        let a1_plus_a2 = &self.c1 + &self.c2;
        let b1_plus_b2 = c1.clone();

        let a0_plus_a1 = &self.c0 + &self.c1;

        // c0 = (NONRESIDUE * ((a1 + a2)*(b1 + b2) - v1 - v2)) + v0
        //    = NONRESIDUE * ((a1 + a2) * b1 - v1)
        let c0 = &(a1_plus_a2 * &b1_plus_b2 - &v1) * P::NONRESIDUE;

        // c1 = (a0 + a1) * (b0 + b1) - v0 - v1 + NONRESIDUE * v2
        //    = (a0 + a1) * b1 - v1
        let c1 = a0_plus_a1 * c1 - &v1;
        // c2 = (a0 + a2) * (b0 + b2) - v0 - v2 + v1
        //    = v1
        let c2 = v1;
        Ok(Self::new(c0, c1, c2))
    }

    /// Multiplies `self` by a sparse element which has `c2 == zero`.
    pub fn mul_by_c0_c1_0(
        &self,
        c0: &Fp2Var<P::Fp2Params>,
        c1: &Fp2Var<P::Fp2Params>,
    ) -> Result<Self, SynthesisError> {
        let v0 = &self.c0 * c0;
        let v1 = &self.c1 * c1;
        // v2 = 0.

        let a1_plus_a2 = &self.c1 + &self.c2;
        let a0_plus_a1 = &self.c0 + &self.c1;
        let a0_plus_a2 = &self.c0 + &self.c2;

        let b1_plus_b2 = c1.clone();
        let b0_plus_b1 = c0 + c1;
        let b0_plus_b2 = c0.clone();

        let c0 = (&a1_plus_a2 * &b1_plus_b2 - &v1) * P::NONRESIDUE + &v0;

        let c1 = a0_plus_a1 * &b0_plus_b1 - &v0 - &v1;

        let c2 = a0_plus_a2 * &b0_plus_b2 - &v0 + &v1;

        Ok(Self::new(c0, c1, c2))
    }
}

impl<P: Fp6Parameters> MulAssign<Fp2<P::Fp2Params>> for Fp6Var<P> {
    fn mul_assign(&mut self, other: Fp2<P::Fp2Params>) {
        self.c0 *= other;
        self.c1 *= other;
        self.c2 *= other;
    }
}
