use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::edwards_on_bls12_381::*;

use crate::edwards_on_bls12_381::FqGadget;

pub type EdwardsGadget = AffineGadget<EdwardsParameters, Fq, FqGadget>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<Fq, _, EdwardsGadget>();
}
