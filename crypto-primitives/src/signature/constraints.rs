use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use crate::signature::SignatureScheme;

pub trait SigRandomizePkGadget<S: SignatureScheme, ConstraintF: Field> {
    type ParametersGadget: AllocGadget<S::Parameters, ConstraintF> + Clone;

    type PublicKeyGadget: ToBytesGadget<ConstraintF>
        + EqGadget<ConstraintF>
        + AllocGadget<S::PublicKey, ConstraintF>
        + Clone;

    fn check_randomization_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        parameters: &Self::ParametersGadget,
        public_key: &Self::PublicKeyGadget,
        randomness: &[UInt8],
    ) -> Result<Self::PublicKeyGadget, SynthesisError>;
}

use crate::signature::FieldBasedSignatureScheme;

pub trait FieldBasedSigGadget<S: FieldBasedSignatureScheme, ConstraintF: Field> {

    type DataGadget:      FieldGadget<ConstraintF, ConstraintF>;
    type SignatureGadget: AllocGadget<S::Signature, ConstraintF>;
    type PublicKeyGadget: AllocGadget<S::PublicKey, ConstraintF>;

    fn check_verify_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        public_key: &Self::PublicKeyGadget,
        signature:  &Self::SignatureGadget,
        message:    &[Self::DataGadget],
    ) -> Result<(), SynthesisError>;
}