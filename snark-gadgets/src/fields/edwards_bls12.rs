use algebra::curves::bls12_377::Bls12_377;

use crate::fields::fp::FpGadget;

// JubJub Fq uses BLS12-377 Fr.
pub type FqGadget = FpGadget<Bls12_377>;
