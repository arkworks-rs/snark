use crate::fields::fp::FpGadget;
use algebra::fields::tweedle::{Fq, Fr};

pub type FqGadget = FpGadget<Fq>;
pub type FrGadget = FpGadget<Fr>;

#[test]
fn test_fq() {
    use crate::fields::tests::*;

    field_test::<_, Fq, FqGadget>();
    frobenius_tests::<Fq, Fq, FqGadget>(13);
    equ_verdict_fp_gadget_test::<Fq>();
    from_bits_fp_gadget_test::<Fq>();
    bit_fp_gadgets_test::<Fq>();
}

#[test]
fn test_fr() {
    use crate::fields::tests::*;

    field_test::<_, Fr, FrGadget>();
    frobenius_tests::<Fr, Fr, FrGadget>(13);
    equ_verdict_fp_gadget_test::<Fr>();
    from_bits_fp_gadget_test::<Fr>();
    bit_fp_gadgets_test::<Fr>();
}
