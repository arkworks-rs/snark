use crate::groups::curves::twisted_edwards::AffineVar;
use algebra::ed_on_bn254::*;

use crate::ed_on_bn254::FqVar;

/// A variable that is the R1CS equivalent of `algebra::ed_on_bn254::EdwardsAffine`.
pub type EdwardsVar = AffineVar<EdwardsParameters, FqVar>;

#[test]
fn test() {
    crate::groups::curves::twisted_edwards::test::<_, EdwardsVar>().unwrap();
}
