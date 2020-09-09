use algebra_core::Field;
use core::fmt::Debug;

use crate::{prf::PRF, Vec};
use r1cs_core::{Namespace, SynthesisError};

use r1cs_std::prelude::*;

pub trait PRFGadget<P: PRF, F: Field> {
    type OutputVar: EqGadget<F>
        + ToBytesGadget<F>
        + AllocVar<P::Output, F>
        + R1CSVar<F, Value = P::Output>
        + Clone
        + Debug;

    fn new_seed(cs: impl Into<Namespace<F>>, seed: &P::Seed) -> Vec<UInt8<F>>;

    fn evaluate(seed: &[UInt8<F>], input: &[UInt8<F>]) -> Result<Self::OutputVar, SynthesisError>;
}
