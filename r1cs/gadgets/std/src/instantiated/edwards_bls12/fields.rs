use crate::fields::fp::FpGadget;
use algebra::fields::edwards_bls12::fq::Fq;

pub type FqGadget = FpGadget<Fq>;

#[test]
fn test() {
    use crate::fields::tests::*;

    field_test::<_, Fq, FqGadget>();
    frobenius_tests::<Fq, Fq, FqGadget>(13);
    equ_verdict_fp_gadget_test::<Fq>();
    from_bits_fp_gadget_test::<Fq>();
    bit_fp_gadgets_test::<Fq>();
}
