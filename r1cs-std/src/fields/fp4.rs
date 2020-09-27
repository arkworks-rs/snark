use crate::fields::{fp2::Fp2Var, quadratic_extension::*};
use algebra::fields::{Fp4Parameters, Fp4ParamsWrapper, QuadExtParameters};

/// A quartic extension field constructed as the tower of a
/// quadratic extension over a quadratic extension field.
/// This is the R1CS equivalent of `algebra_core::Fp4<P>`.
pub type Fp4Var<P> = QuadExtVar<Fp2Var<<P as Fp4Parameters>::Fp2Params>, Fp4ParamsWrapper<P>>;

impl<P: Fp4Parameters> QuadExtVarParams<Fp2Var<P::Fp2Params>> for Fp4ParamsWrapper<P> {
    fn mul_base_field_var_by_frob_coeff(fe: &mut Fp2Var<P::Fp2Params>, power: usize) {
        fe.c0 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        fe.c1 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}
