use algebra::fields::bn_382::Fr;

use crate::fields::fp::FpGadget;

pub type FqGadget = FpGadget<Fr>;

#[test]
fn bn382_g_field_gadgets_test() {
    use super::*;
    use crate::fields::tests::*;

    field_test::<_, Fr, FqGadget>();
    frobenius_tests::<Fr, Fr, FqGadget>(13);
    equ_verdict_fp_gadget_test::<Fr>();
    from_bits_fp_gadget_test::<Fr>();
    bit_fp_gadgets_test::<Fr>();
}
