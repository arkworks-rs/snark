use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::edwards_bls12::*;

use crate::edwards_bls12::FqGadget;

pub type EdwardsBlsGadget = AffineGadget<EdwardsParameters, Fq, FqGadget>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<_, EdwardsParameters, EdwardsBlsGadget>();
}
