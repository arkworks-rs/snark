use crate::groups::curves::twisted_edwards::AffineVar;
use algebra::ed_on_bls12_377::*;

use crate::ed_on_bls12_377::FqVar;

pub type EdwardsVar = AffineVar<EdwardsParameters, FqVar>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<EdwardsParameters, EdwardsVar>().unwrap();
}
