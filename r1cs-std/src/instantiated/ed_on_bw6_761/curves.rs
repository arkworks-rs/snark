use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::ed_on_cp6_782::*;

use crate::ed_on_cp6_782::FqGadget;

pub type EdwardsGadget = AffineGadget<EdwardsParameters, Fq, FqGadget>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<_, EdwardsParameters, EdwardsGadget>();
}
