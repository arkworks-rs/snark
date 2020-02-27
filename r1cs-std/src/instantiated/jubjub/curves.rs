use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::jubjub::*;

use crate::jubjub::FqGadget;

pub type JubJubGadget = AffineGadget<JubJubParameters, Fq, FqGadget>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<Fq, _, JubJubGadget>();
}
