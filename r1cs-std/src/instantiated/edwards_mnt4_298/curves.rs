use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::edwards_mnt4_298::*;

use crate::instantiated::edwards_mnt4_298::fields::FqGadget;

pub type EdwardsMnt4298gadget = AffineGadget<EdwardsParameters, Fq, FqGadget>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<_, EdwardsParameters, EdwardsMnt4298gadget>();
}
