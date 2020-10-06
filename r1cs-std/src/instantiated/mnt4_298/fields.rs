use algebra::mnt4_298::{Fq, Fq2Parameters, Fq4Parameters};

use crate::fields::{fp::FpVar, fp2::Fp2Var, fp4::Fp4Var};

/// A variable that is the R1CS equivalent of `algebra::mnt4_298::Fq`.
pub type FqVar = FpVar<Fq>;
/// A variable that is the R1CS equivalent of `algebra::mnt4_298::Fq2`.
pub type Fq2Var = Fp2Var<Fq2Parameters>;
/// A variable that is the R1CS equivalent of `algebra::mnt4_298::Fq4`.
pub type Fq4Var = Fp4Var<Fq4Parameters>;

#[test]
fn mnt4_298_field_gadgets_test() {
    use super::*;
    use crate::fields::tests::*;
    use algebra::mnt4_298::{Fq, Fq2, Fq4};

    field_test::<_, _, FqVar>().unwrap();
    frobenius_tests::<Fq, _, FqVar>(13).unwrap();

    field_test::<_, _, Fq2Var>().unwrap();
    frobenius_tests::<Fq2, _, Fq2Var>(13).unwrap();

    field_test::<_, _, Fq4Var>().unwrap();
    frobenius_tests::<Fq4, _, Fq4Var>(13).unwrap();
}
