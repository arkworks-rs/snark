use algebra::PairingEngine;

use crate::gadgets::commitment::{
    pedersen::{
        PedersenCommitmentGadget, PedersenCommitmentGadgetParameters, PedersenRandomnessGadget,
    },
    CommitmentGadget,
};
pub use crate::gadgets::crh::injective_map::InjectiveMapGadget;
use algebra::groups::Group;
use snark::{ConstraintSystem, SynthesisError};
use snark_gadgets::{groups::GroupGadget, uint8::UInt8};

use crate::crypto_primitives::commitment::{
    injective_map::{InjectiveMap, PedersenCommCompressor},
    pedersen::PedersenWindow,
};
use std::marker::PhantomData;

pub struct PedersenCommitmentCompressorGadget<
    G: Group,
    I: InjectiveMap<G>,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
    IG: InjectiveMapGadget<G, I, E, GG>,
> {
    _compressor:        PhantomData<I>,
    _compressor_gadget: PhantomData<IG>,
    _crh:               PedersenCommitmentGadget<G, E, GG>,
}

impl<G, I, E, GG, IG, W> CommitmentGadget<PedersenCommCompressor<G, I, W>, E>
    for PedersenCommitmentCompressorGadget<G, I, E, GG, IG>
where
    G: Group,
    I: InjectiveMap<G>,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
    IG: InjectiveMapGadget<G, I, E, GG>,
    W: PedersenWindow,
{
    type OutputGadget = IG::OutputGadget;
    type ParametersGadget = PedersenCommitmentGadgetParameters<G, W, E>;
    type RandomnessGadget = PedersenRandomnessGadget;

    fn check_commitment_gadget<CS: ConstraintSystem<E>>(
        mut cs: CS,
        parameters: &Self::ParametersGadget,
        input: &[UInt8],
        r: &Self::RandomnessGadget,
    ) -> Result<Self::OutputGadget, SynthesisError> {
        let result = PedersenCommitmentGadget::<G, E, GG>::check_commitment_gadget(
            cs.ns(|| "PedersenComm"),
            parameters,
            input,
            r,
        )?;
        IG::evaluate_map(cs.ns(|| "InjectiveMap"), &result)
    }
}
