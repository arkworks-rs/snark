use algebra::PairingEngine;
use std::fmt::Debug;

use crate::crypto_primitives::prf::PRF;
use snark::{ConstraintSystem, SynthesisError};

use snark_gadgets::{
    uint8::UInt8,
    utils::{AllocGadget, EqGadget, ToBytesGadget},
};

pub mod blake2s;

pub trait PRFGadget<P: PRF, E: PairingEngine> {
    type OutputGadget: EqGadget<E> + ToBytesGadget<E> + AllocGadget<P::Output, E> + Clone + Debug;

    fn new_seed<CS: ConstraintSystem<E>>(cs: CS, output: &P::Seed) -> Vec<UInt8>;

    fn check_evaluation_gadget<CS: ConstraintSystem<E>>(
        cs: CS,
        seed: &[UInt8],
        input: &[UInt8],
    ) -> Result<Self::OutputGadget, SynthesisError>;
}
