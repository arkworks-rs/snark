use crate::vrf::FieldBasedVrf;
use r1cs_std::fields::FieldGadget;
use r1cs_std::alloc::AllocGadget;
use r1cs_core::{ConstraintSystem, SynthesisError};
use algebra::Field;

pub trait FieldBasedVrfGadget<S: FieldBasedVrf, ConstraintF: Field> {

    type DataGadget:      FieldGadget<ConstraintF, ConstraintF>;
    type ProofGadget:     AllocGadget<S::Proof, ConstraintF>;
    type PublicKeyGadget: AllocGadget<S::PublicKey, ConstraintF>;

    fn check_verify_gadget<CS: ConstraintSystem<ConstraintF>>(
        cs:         CS,
        public_key: &Self::PublicKeyGadget,
        proof:      &Self::ProofGadget,
        message:    &[Self::DataGadget],
    ) -> Result<Self::DataGadget, SynthesisError>;
}