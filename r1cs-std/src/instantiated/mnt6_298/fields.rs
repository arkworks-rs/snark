use algebra::mnt6_298::{Fq, Fq3Parameters, Fq6Parameters};

use crate::fields::{fp::FpVar, fp3::Fp3Var, fp6_2over3::Fp6Var};

/// A variable that is the R1CS equivalent of `algebra::mnt4_298::Fq`.
pub type FqVar = FpVar<Fq>;
/// A variable that is the R1CS equivalent of `algebra::mnt4_298::Fq3`.
pub type Fq3Var = Fp3Var<Fq3Parameters>;
/// A variable that is the R1CS equivalent of `algebra::mnt4_298::Fq6`.
pub type Fq6Var = Fp6Var<Fq6Parameters>;

#[test]
fn mnt6_298_field_gadgets_test() {
    use super::*;
    use crate::fields::tests::*;
    use algebra::mnt6_298::{Fq, Fq3, Fq6};

    field_test::<_, _, FqVar>().unwrap();
    frobenius_tests::<Fq, _, FqVar>(13).unwrap();

    field_test::<_, _, Fq3Var>().unwrap();
    frobenius_tests::<Fq3, _, Fq3Var>(13).unwrap();

    field_test::<_, _, Fq6Var>().unwrap();
    frobenius_tests::<Fq6, _, Fq6Var>(13).unwrap();
}
