use crate::groups::curves::twisted_edwards::AffineVar;
use algebra::ed_on_mnt4_753::*;

use crate::instantiated::ed_on_mnt4_753::fields::FqVar;

pub type EdwardsVar = AffineVar<EdwardsParameters, FqVar>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<EdwardsParameters, EdwardsVar>().unwrap();
}
