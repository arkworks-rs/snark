use crate::groups::mnt6;
use algebra::mnt6_753::Parameters;

/// An element of G1 in the MNT6-753 bilinear group.
pub type G1Var = mnt6::G1Var<Parameters>;
/// An element of G2 in the MNT6-753 bilinear group.
pub type G2Var = mnt6::G2Var<Parameters>;

/// Represents the cached precomputation that can be performed on a G1 element
/// which enables speeding up pairing computation.
pub type G1PreparedVar = mnt6::G1PreparedVar<Parameters>;
/// Represents the cached precomputation that can be performed on a G2 element
/// which enables speeding up pairing computation.
pub type G2PreparedVar = mnt6::G2PreparedVar<Parameters>;

#[test]
fn test() {
    use algebra::curves::models::mnt6::MNT6Parameters;
    crate::groups::curves::short_weierstrass::test::<
        <Parameters as MNT6Parameters>::G1Parameters,
        G1Var,
    >()
    .unwrap();
    crate::groups::curves::short_weierstrass::test::<
        <Parameters as MNT6Parameters>::G2Parameters,
        G2Var,
    >()
    .unwrap();
}
