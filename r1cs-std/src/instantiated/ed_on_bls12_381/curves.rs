use crate::groups::curves::twisted_edwards::AffineVar;
use algebra::ed_on_bls12_381::*;

use crate::ed_on_bls12_381::FqVar;

pub type EdwardsVar = AffineVar<EdwardsParameters, FqVar>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<_, EdwardsVar>().unwrap();
}
