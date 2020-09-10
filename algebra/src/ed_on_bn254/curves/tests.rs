#![allow(unused_imports)]
use crate::ed_on_bn254::*;
use algebra_core::{FromBytes, ToBytes, Zero};
use core::str::FromStr;
edwards_curve_tests!();

#[test]
#[cfg(feature = "all_tests")]
fn test_scalar_multiplication() {
    println!("Started getting field elements");
    let f1 = Fr::from_str(
        "4691331900926794624732159288782398864809513177368446695323460897088210774597",
    )
    .unwrap();
    let f2 = Fr::from_str(
        "1305028103380024953477151132159456965337646722479526711736847301646466538045",
    )
    .unwrap();

    println!("Finished getting field elements");
    let g = EdwardsAffine::from_str(
        "(15863623088992515880085393097393553694825975317405843389771115419751650972659, \
         16950150798460657717958625567821834550301663161624707787222815936182638968203)",
    )
    .unwrap();
    let f1f2g = EdwardsAffine::from_str(
        "(20773645713088336957786354488799297695596635653208610804806657050882264237947, \
         19987327827845206670850937090314462639017692512983955920885166014935289314257)",
    )
    .unwrap();

    println!("Finished getting group elements");

    assert!(!g.is_zero());
    assert!(!f1f2g.is_zero());

    let f1g = g.mul(f1).into_affine();
    assert_eq!(g.mul(f1 * &f2).into_affine(), f1f2g);
    assert_eq!(f1g.mul(f2).into_affine(), f1f2g);
}

#[test]
#[cfg(feature = "all_tests")]
fn test_bytes() {
    let g_from_repr = EdwardsAffine::from_str(
        "(15863623088992515880085393097393553694825975317405843389771115419751650972659, \
         16950150798460657717958625567821834550301663161624707787222815936182638968203)",
    )
    .unwrap();

    let g_bytes = algebra_core::to_bytes![g_from_repr].unwrap();
    let g = EdwardsAffine::read(g_bytes.as_slice()).unwrap();
    assert_eq!(g_from_repr, g);
}
