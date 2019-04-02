use algebra::{
    curves::sw6::SW6,
    fields::bls12_377::{Fq12Parameters, Fq2Parameters, Fq6Parameters},
};

use super::{fp::FpGadget, fp12::Fp12Gadget, fp2::Fp2Gadget, fp6_3over2::Fp6Gadget};

pub type FqGadget = FpGadget<SW6>;
pub type Fq2Gadget = Fp2Gadget<Fq2Parameters, SW6>;
pub type Fq6Gadget = Fp6Gadget<Fq6Parameters, SW6>;
pub type Fq12Gadget = Fp12Gadget<Fq12Parameters, SW6>;
