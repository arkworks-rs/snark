use crate::pairing::bls12::PairingGadget as Bls12PG;
use algebra::curves::bls12_377::Bls12_377Parameters;

pub type PairingGadget = Bls12PG<Bls12_377Parameters>;
