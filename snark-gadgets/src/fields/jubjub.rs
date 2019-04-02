use algebra::curves::bls12_381::Bls12_381;

use crate::fields::fp::FpGadget;

// JubJub Fq uses BLS12-381 Fr.
pub type FqGadget = FpGadget<Bls12_381>;
