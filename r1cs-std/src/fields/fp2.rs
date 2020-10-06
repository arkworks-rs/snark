use crate::fields::{fp::FpVar, quadratic_extension::*};
use algebra::fields::{Fp2Parameters, Fp2ParamsWrapper, QuadExtParameters};

/// A quadratic extension field constructed over a prime field.
/// This is the R1CS equivalent of `algebra_core::Fp2<P>`.
pub type Fp2Var<P> = QuadExtVar<FpVar<<P as Fp2Parameters>::Fp>, Fp2ParamsWrapper<P>>;

impl<P: Fp2Parameters> QuadExtVarParams<FpVar<P::Fp>> for Fp2ParamsWrapper<P> {
    fn mul_base_field_var_by_frob_coeff(fe: &mut FpVar<P::Fp>, power: usize) {
        *fe *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}
