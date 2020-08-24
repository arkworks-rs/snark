use crate::groups::bls12;
use algebra::bls12_377::Parameters;

pub type G1Var = bls12::G1Var<Parameters>;
pub type G2Var = bls12::G2Var<Parameters>;

pub type G1PreparedVar = bls12::G1PreparedVar<Parameters>;
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
