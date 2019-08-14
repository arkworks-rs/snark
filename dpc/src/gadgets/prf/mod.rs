use algebra::Field;
use std::fmt::Debug;

use crate::crypto_primitives::prf::PRF;
use r1cs_core::{ConstraintSystem, SynthesisError};

use r1cs_std::prelude::*;

pub mod blake2s;

pub trait PRFGadget<P: PRF, ConstraintF: Field> {
    type OutputGadget: EqGadget<ConstraintF> + ToBytesGadget<ConstraintF> + AllocGadget<P::Output, ConstraintF> + Clone + Debug;

    fn new_seed<CS: ConstraintSystem<ConstraintF>>(cs: CS, output: &P::Seed) -> Vec<UInt8>;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        seed: &[UInt8],
        input: &[UInt8],
    ) -> Result<Self::OutputGadget, SynthesisError>;
}
