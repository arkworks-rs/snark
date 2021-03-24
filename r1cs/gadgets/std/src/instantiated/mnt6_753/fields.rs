use algebra::fields::mnt6753::{Fq, Fq3Parameters, Fq6Parameters};

use crate::fields::{fp::FpGadget, fp3::Fp3Gadget, fp6_2over3::Fp6Gadget};

pub type FqGadget = FpGadget<Fq>;
pub type Fq3Gadget = Fp3Gadget<Fq3Parameters, Fq>;
pub type Fq6Gadget = Fp6Gadget<Fq6Parameters, Fq>;

#[test]
fn mnt6_753_field_gadgets_test() {
    use super::*;
    use crate::fields::tests::*;
    use algebra::fields::mnt6753::{Fq, Fq3, Fq6};

    field_test::<_, Fq, FqGadget>();
    frobenius_tests::<Fq, Fq, FqGadget>(13);
    equ_verdict_fp_gadget_test::<Fq>();
    from_bits_fp_gadget_test::<Fq>();
    bit_fp_gadgets_test::<Fq>();

    field_test::<_, Fq, Fq3Gadget>();
    frobenius_tests::<Fq3, Fq, Fq3Gadget>(13);

    field_test::<_, Fq, Fq6Gadget>();
    frobenius_tests::<Fq6, Fq, Fq6Gadget>(13);
}
