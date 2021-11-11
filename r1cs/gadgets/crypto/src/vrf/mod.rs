use algebra::Field;
use primitives::vrf::FieldBasedVrf;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::alloc::AllocGadget;
use r1cs_std::fields::FieldGadget;

pub mod ecvrf;

pub trait FieldBasedVrfGadget<S: FieldBasedVrf, ConstraintF: Field> {
    type DataGadget: FieldGadget<ConstraintF, ConstraintF>;
    type ProofGadget: AllocGadget<S::Proof, ConstraintF>;
    type PublicKeyGadget: AllocGadget<S::PublicKey, ConstraintF>;
    type GHParametersGadget: AllocGadget<S::GHParams, ConstraintF>;

    fn enforce_proof_to_hash_verification<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        pp: &Self::GHParametersGadget,
        public_key: &Self::PublicKeyGadget,
        proof: &Self::ProofGadget,
        message: Self::DataGadget,
    ) -> Result<Self::DataGadget, SynthesisError>;
}
