use crate::fields::{fp3::Fp3Var, quadratic_extension::*};
use algebra::fields::{fp6_2over3::*, QuadExtParameters};

/// A sextic extension field constructed as the tower of a
/// quadratic extension over a cubic extension field.
/// This is the R1CS equivalent of `algebra_core::fp6_2over3::Fp6<P>`.
pub type Fp6Var<P> = QuadExtVar<Fp3Var<<P as Fp6Parameters>::Fp3Params>, Fp6ParamsWrapper<P>>;

impl<P: Fp6Parameters> QuadExtVarParams<Fp3Var<P::Fp3Params>> for Fp6ParamsWrapper<P> {
    fn mul_base_field_var_by_frob_coeff(fe: &mut Fp3Var<P::Fp3Params>, power: usize) {
        fe.c0 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        fe.c1 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        fe.c2 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}
