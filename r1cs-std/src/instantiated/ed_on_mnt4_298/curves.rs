use crate::groups::curves::twisted_edwards::AffineVar;
use algebra::ed_on_mnt4_298::*;

use crate::instantiated::ed_on_mnt4_298::fields::FqVar;

pub type EdwardsVar = AffineVar<EdwardsParameters, FqVar>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<EdwardsParameters, EdwardsVar>().unwrap();
}
