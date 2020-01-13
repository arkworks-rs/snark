use algebra::Field;
use r1cs_core::{R1CS, SynthesisError};
use r1cs_std::prelude::*;

use crate::signature::SignatureScheme;

pub trait SigRandomizePkGadget<S: SignatureScheme, ConstraintF: Field> {
    type ParametersGadget: AllocGadget<S::Parameters, ConstraintF> + Clone;

    type PublicKeyGadget: ToBytesGadget<ConstraintF>
        + EqGadget<ConstraintF>
        + AllocGadget<S::PublicKey, ConstraintF>
        + Clone;

    fn check_randomization_gadget<CS: R1CS<ConstraintF>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        public_key: &Self::PublicKeyGadget,
        randomness: &[UInt8],
    ) -> Result<Self::PublicKeyGadget, SynthesisError>;
}
