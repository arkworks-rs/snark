use algebra::fields::jubjub::fq::Fq;

use crate::fields::fp::FpGadget;

// JubJub Fq uses BLS12-381 Fr.
pub type FqGadget = FpGadget<Fq>;
