use crate::groups::curves::twisted_edwards::AffineGadget;
use algebra::edwards_sw6::*;

use crate::edwards_sw6::FqGadget;

pub type EdwardsSWGadget = AffineGadget<EdwardsParameters, Fq, FqGadget>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<_, EdwardsParameters, EdwardsSWGadget>();
}
