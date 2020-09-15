use algebra::mnt6_753::{Fq, Fq3Parameters, Fq6Parameters};

use crate::fields::{fp::FpVar, fp3::Fp3Var, fp6_2over3::Fp6Var};

pub type FqVar = FpVar<Fq>;
pub type Fq3Var = Fp3Var<Fq3Parameters>;
pub type Fq6Var = Fp6Var<Fq6Parameters>;

#[test]
fn mnt6_753_field_gadgets_test() {
    use super::*;
    use crate::fields::tests::*;
    use algebra::mnt6_753::{Fq, Fq3, Fq6};

    field_test::<_, _, FqVar>().unwrap();
    frobenius_tests::<Fq, _, FqVar>(13).unwrap();

    field_test::<_, _, Fq3Var>().unwrap();
    frobenius_tests::<Fq3, _, Fq3Var>(13).unwrap();

    field_test::<_, _, Fq6Var>().unwrap();
    frobenius_tests::<Fq6, _, Fq6Var>(13).unwrap();
}
