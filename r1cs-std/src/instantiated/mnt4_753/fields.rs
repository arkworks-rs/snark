use algebra::mnt4_753::{Fq, Fq2Parameters, Fq4Parameters};

use crate::fields::{fp::FpGadget, fp2::Fp2Gadget, fp4::Fp4Gadget};

pub type FqGadget = FpGadget<Fq>;
pub type Fq2Gadget = Fp2Gadget<Fq2Parameters, Fq>;
pub type Fq4Gadget = Fp4Gadget<Fq4Parameters, Fq>;

#[test]
fn mnt4_753_field_gadgets_test() {
    use super::*;
    use crate::fields::tests::*;
    use algebra::mnt4_753::{Fq, Fq2, Fq4};

    field_test::<_, Fq, FqGadget>();
    frobenius_tests::<Fq, Fq, FqGadget>(13);

    field_test::<_, Fq, Fq2Gadget>();
    frobenius_tests::<Fq2, Fq, Fq2Gadget>(13);

    field_test::<_, Fq, Fq4Gadget>();
    frobenius_tests::<Fq4, Fq, Fq4Gadget>(13);
}
