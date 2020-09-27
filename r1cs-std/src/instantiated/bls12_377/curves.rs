use crate::groups::bls12;
use algebra::bls12_377::Parameters;

/// An element of G1 in the BLS12-377 bilinear group.
pub type G1Var = bls12::G1Var<Parameters>;
/// An element of G2 in the BLS12-377 bilinear group.
pub type G2Var = bls12::G2Var<Parameters>;

/// Represents the cached precomputation that can be performed on a G1 element
/// which enables speeding up pairing computation.
pub type G1PreparedVar = bls12::G1PreparedVar<Parameters>;
/// Represents the cached precomputation that can be performed on a G2 element
/// which enables speeding up pairing computation.
pub type G2PreparedVar = bls12::G2PreparedVar<Parameters>;

#[test]
fn test() {
    use algebra::curves::models::bls12::Bls12Parameters;
    crate::groups::curves::short_weierstrass::test::<
        <Parameters as Bls12Parameters>::G1Parameters,
        G1Var,
    >()
    .unwrap();
    crate::groups::curves::short_weierstrass::test::<
        <Parameters as Bls12Parameters>::G2Parameters,
        G2Var,
    >()
    .unwrap();
}
