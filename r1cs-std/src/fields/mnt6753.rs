use algebra::fields::mnt6753::{Fq, Fq3b7Parameters, Fq6b7Parameters};

use super::{fp::FpGadget, fp3::Fp3Gadget, fp6_2over3::Fp6Gadget};

pub type FqGadget = FpGadget<Fq>;
pub type Fq3Gadget = Fp3Gadget<Fq3b7Parameters, Fq>;
pub type Fq6Gadget = Fp6Gadget<Fq6b7Parameters, Fq>;
