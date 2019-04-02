use algebra::PairingEngine;
use snark::{ConstraintSystem, SynthesisError};
use snark_gadgets::{
    uint8::UInt8,
    utils::{AllocGadget, EqGadget, ToBytesGadget},
};

use crate::crypto_primitives::signature::SignatureScheme;

pub mod schnorr;

pub trait SigRandomizePkGadget<S: SignatureScheme, E: PairingEngine> {
    type ParametersGadget: AllocGadget<S::Parameters, E> + Clone;

    type PublicKeyGadget: ToBytesGadget<E> + EqGadget<E> + AllocGadget<S::PublicKey, E> + Clone;

    fn check_randomization_gadget<CS: ConstraintSystem<E>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        public_key: &Self::PublicKeyGadget,
        randomness: &[UInt8],
    ) -> Result<Self::PublicKeyGadget, SynthesisError>;
}
