use algebra::curves::sw6::SW6;

use crate::fields::fp::FpGadget;

// Edwards_SW6 Fq uses SW6 Fr.
pub type FqGadget = FpGadget<SW6>;
