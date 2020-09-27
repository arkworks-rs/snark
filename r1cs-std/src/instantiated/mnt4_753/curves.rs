use crate::groups::mnt4;
use algebra::mnt4_753::Parameters;

/// An element of G1 in the MNT4-753 bilinear group.
pub type G1Var = mnt4::G1Var<Parameters>;
/// An element of G2 in the MNT4-753 bilinear group.
pub type G2Var = mnt4::G2Var<Parameters>;

/// Represents the cached precomputation that can be performed on a G1 element
/// which enables speeding up pairing computation.
pub type G1PreparedVar = mnt4::G1PreparedVar<Parameters>;
/// Represents the cached precomputation that can be performed on a G2 element
/// which enables speeding up pairing computation.
pub type G2PreparedVar = mnt4::G2PreparedVar<Parameters>;

#[test]
fn test() {
    use algebra::curves::models::mnt4::MNT4Parameters;
    crate::groups::curves::short_weierstrass::test::<
        <Parameters as MNT4Parameters>::G1Parameters,
        G1Var,
    >()
    .unwrap();
    crate::groups::curves::short_weierstrass::test::<
        <Parameters as MNT4Parameters>::G2Parameters,
        G2Var,
    >()
    .unwrap();
}
