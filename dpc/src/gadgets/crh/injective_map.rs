use algebra::PairingEngine;
use std::{fmt::Debug, marker::PhantomData};

use crate::gadgets::crh::{
    pedersen::{PedersenCRHGadget, PedersenCRHGadgetParameters},
    FixedLengthCRHGadget,
};
use algebra::{
    curves::{
        models::{ModelParameters, TEModelParameters},
        twisted_edwards_extended::{GroupAffine as TEAffine, GroupProjective as TEProjective},
    },
    groups::Group,
};
use snark::{ConstraintSystem, SynthesisError};
use snark_gadgets::{
    fields::fp::FpGadget,
    groups::{curves::twisted_edwards::AffineGadget as TwistedEdwardsGadget, GroupGadget},
    uint8::UInt8,
    utils::{AllocGadget, CondSelectGadget, EqGadget, ToBytesGadget},
};

use crate::crypto_primitives::crh::{
    injective_map::{InjectiveMap, PedersenCRHCompressor, TECompressor},
    pedersen::{PedersenCRH, PedersenWindow},
};

pub trait InjectiveMapGadget<G: Group, I: InjectiveMap<G>, E: PairingEngine, GG: GroupGadget<G, E>>
{
    type OutputGadget: EqGadget<E>
        + ToBytesGadget<E>
        + CondSelectGadget<E>
        + AllocGadget<I::Output, E>
        + Debug
        + Clone
        + Sized;

    fn evaluate_map<CS: ConstraintSystem<E>>(
        cs: CS,
        ge: &GG,
    ) -> Result<Self::OutputGadget, SynthesisError>;
    fn cost() -> usize;
}

pub struct TECompressorGadget;

impl<E, P> InjectiveMapGadget<TEAffine<P>, TECompressor, E, TwistedEdwardsGadget<P, E, FpGadget<E>>>
    for TECompressorGadget
where
    E: PairingEngine,
    P: TEModelParameters + ModelParameters<BaseField = E::Fr>,
{
    type OutputGadget = FpGadget<E>;

    fn evaluate_map<CS: ConstraintSystem<E>>(
        _cs: CS,
        ge: &TwistedEdwardsGadget<P, E, FpGadget<E>>,
    ) -> Result<Self::OutputGadget, SynthesisError> {
        Ok(ge.x.clone())
    }

    fn cost() -> usize {
        0
    }
}

impl<E, P>
    InjectiveMapGadget<TEProjective<P>, TECompressor, E, TwistedEdwardsGadget<P, E, FpGadget<E>>>
    for TECompressorGadget
where
    E: PairingEngine,
    P: TEModelParameters + ModelParameters<BaseField = E::Fr>,
{
    type OutputGadget = FpGadget<E>;

    fn evaluate_map<CS: ConstraintSystem<E>>(
        _cs: CS,
        ge: &TwistedEdwardsGadget<P, E, FpGadget<E>>,
    ) -> Result<Self::OutputGadget, SynthesisError> {
        Ok(ge.x.clone())
    }

    fn cost() -> usize {
        0
    }
}

pub struct PedersenCRHCompressorGadget<
    G: Group,
    I: InjectiveMap<G>,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
    IG: InjectiveMapGadget<G, I, E, GG>,
> {
    _compressor:        PhantomData<I>,
    _compressor_gadget: PhantomData<IG>,
    _crh:               PedersenCRHGadget<G, E, GG>,
}

impl<G, I, E, GG, IG, W> FixedLengthCRHGadget<PedersenCRHCompressor<G, I, W>, E>
    for PedersenCRHCompressorGadget<G, I, E, GG, IG>
where
    G: Group,
    I: InjectiveMap<G>,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
    IG: InjectiveMapGadget<G, I, E, GG>,
    W: PedersenWindow,
{
    type OutputGadget = IG::OutputGadget;
    type ParametersGadget = PedersenCRHGadgetParameters<G, W, E, GG>;

    fn check_evaluation_gadget<CS: ConstraintSystem<E>>(
        mut cs: CS,
        parameters: &Self::ParametersGadget,
        input: &[UInt8],
    ) -> Result<Self::OutputGadget, SynthesisError> {
        let result = PedersenCRHGadget::<G, E, GG>::check_evaluation_gadget(
            cs.ns(|| "PedCRH"),
            parameters,
            input,
        )?;
        IG::evaluate_map(cs.ns(|| "InjectiveMap"), &result)
    }

    fn cost() -> usize {
        <PedersenCRHGadget<G, E, GG> as FixedLengthCRHGadget<PedersenCRH<G, W>, E>>::cost()
            + IG::cost()
    }
}
