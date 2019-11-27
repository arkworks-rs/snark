use algebra::fields::mnt4753::{Fq, Fq2Parameters, Fq4Parameters};

use super::{fp::FpGadget, fp4::Fp4Gadget, fp2::Fp2Gadget};

pub type FqGadget = FpGadget<Fq>;
pub type Fq2Gadget = Fp2Gadget<Fq2Parameters, Fq>;
pub type Fq4Gadget = Fp4Gadget<Fq4Parameters, Fq>;
