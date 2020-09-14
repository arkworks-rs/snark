use algebra::bls12_377::{Fq, Fq12Parameters, Fq2Parameters, Fq6Parameters};

use crate::fields::{fp::FpVar, fp12::Fp12Var, fp2::Fp2Var, fp6_3over2::Fp6Var};

/// A variable that is the R1CS equivalent of `algebra::bls12_377::Fq`.
pub type FqVar = FpVar<Fq>;

/// A variable that is the R1CS equivalent of `algebra::bls12_377::Fq2`.
pub type Fq2Var = Fp2Var<Fq2Parameters>;
/// A variable that is the R1CS equivalent of `algebra::bls12_377::Fq6`.
pub type Fq6Var = Fp6Var<Fq6Parameters>;
/// A variable that is the R1CS equivalent of `algebra::bls12_377::Fq12`.
pub type Fq12Var = Fp12Var<Fq12Parameters>;

#[test]
fn bls12_377_field_test() {
    use super::*;
    use crate::fields::tests::*;
    use algebra::bls12_377::{Fq, Fq12, Fq2, Fq6};

    field_test::<_, _, FqVar>().unwrap();
    frobenius_tests::<Fq, _, FqVar>(13).unwrap();

    field_test::<_, _, Fq2Var>().unwrap();
    frobenius_tests::<Fq2, _, Fq2Var>(13).unwrap();

    field_test::<_, _, Fq6Var>().unwrap();
    frobenius_tests::<Fq6, _, Fq6Var>(13).unwrap();

    field_test::<_, _, Fq12Var>().unwrap();
    frobenius_tests::<Fq12, _, Fq12Var>(13).unwrap();
}
