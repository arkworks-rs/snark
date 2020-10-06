use crate::fields::{cubic_extension::*, fp::FpVar};
use algebra::fields::{CubicExtParameters, Fp3Parameters, Fp3ParamsWrapper};

/// A cubic extension field constructed over a prime field.
/// This is the R1CS equivalent of `algebra_core::Fp3<P>`.
pub type Fp3Var<P> = CubicExtVar<FpVar<<P as Fp3Parameters>::Fp>, Fp3ParamsWrapper<P>>;

impl<P: Fp3Parameters> CubicExtVarParams<FpVar<P::Fp>> for Fp3ParamsWrapper<P> {
    fn mul_base_field_vars_by_frob_coeff(
        c1: &mut FpVar<P::Fp>,
        c2: &mut FpVar<P::Fp>,
        power: usize,
    ) {
        *c1 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        *c2 *= Self::FROBENIUS_COEFF_C2[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}
