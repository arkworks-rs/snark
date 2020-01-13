use algebra::Field;
use std::fmt::Debug;

use crate::prf::PRF;
use r1cs_core::{R1CS, SynthesisError};

use r1cs_std::prelude::*;

pub trait PRFGadget<P: PRF, ConstraintF: Field> {
    type OutputGadget: EqGadget<ConstraintF>
        + ToBytesGadget<ConstraintF>
        + AllocGadget<P::Output, ConstraintF>
        + Clone
        + Debug;

    fn new_seed<CS: R1CS<ConstraintF>>(cs: CS, output: &P::Seed) -> Vec<UInt8>;

    fn check_evaluation_gadget<CS: R1CS<ConstraintF>>(
        cs: CS,
        seed: &[UInt8],
        input: &[UInt8],
    ) -> Result<Self::OutputGadget, SynthesisError>;
}
