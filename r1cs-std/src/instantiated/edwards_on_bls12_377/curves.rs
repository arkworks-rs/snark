use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::edwards_on_bls12_377::*;

use crate::edwards_on_bls12_377::FqGadget;

pub type EdwardsGadget = AffineGadget<EdwardsParameters, Fq, FqGadget>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<_, EdwardsParameters, EdwardsGadget>();
}
