use std::{fmt::Debug, marker::PhantomData};

use crate::crh::{
    FixedLengthCRHGadget,
    injective_map::{InjectiveMap, PedersenCRHCompressor, TECompressor},
    pedersen::{
        PedersenWindow,
        constraints::{PedersenCRHGadget, PedersenCRHGadgetParameters},
    }
};

use algebra::{
    curves::{
        models::{ModelParameters, TEModelParameters},
        twisted_edwards_extended::{GroupAffine as TEAffine, GroupProjective as TEProjective},
    },
    fields::{Field, PrimeField, SquareRootField},
    groups::Group,
};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{
    fields::fp::FpGadget,
    groups::{curves::twisted_edwards::AffineGadget as TwistedEdwardsGadget, GroupGadget},
    prelude::*,
};

pub trait InjectiveMapGadget<G: Group, I: InjectiveMap<G>, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>> {
    type OutputGadget: EqGadget<ConstraintF>
        + ToBytesGadget<ConstraintF>
        + CondSelectGadget<ConstraintF>
        + AllocGadget<I::Output, ConstraintF>
        + Debug
        + Clone
        + Sized;

    fn evaluate_map<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        ge: &GG,
    ) -> Result<Self::OutputGadget, SynthesisError>;
}

pub struct TECompressorGadget;

impl<ConstraintF, P> InjectiveMapGadget<TEAffine<P>, TECompressor, ConstraintF, TwistedEdwardsGadget<P, ConstraintF, FpGadget<ConstraintF>>>
    for TECompressorGadget
where
    ConstraintF: PrimeField + SquareRootField,
    P: TEModelParameters + ModelParameters<BaseField = ConstraintF>,
{
    type OutputGadget = FpGadget<ConstraintF>;

    fn evaluate_map<CS: ConstraintSystem<ConstraintF>>(
        _cs: CS,
        ge: &TwistedEdwardsGadget<P, ConstraintF, FpGadget<ConstraintF>>,
    ) -> Result<Self::OutputGadget, SynthesisError> {
        Ok(ge.x.clone())
    }
}

impl<ConstraintF, P>
    InjectiveMapGadget<TEProjective<P>, TECompressor, ConstraintF, TwistedEdwardsGadget<P, ConstraintF, FpGadget<ConstraintF>>>
    for TECompressorGadget
where
    ConstraintF: PrimeField + SquareRootField,
    P: TEModelParameters + ModelParameters<BaseField = ConstraintF>,
{
    type OutputGadget = FpGadget<ConstraintF>;

    fn evaluate_map<CS: ConstraintSystem<ConstraintF>>(
        _cs: CS,
        ge: &TwistedEdwardsGadget<P, ConstraintF, FpGadget<ConstraintF>>,
    ) -> Result<Self::OutputGadget, SynthesisError> {
        Ok(ge.x.clone())
    }
}

pub struct PedersenCRHCompressorGadget<G, I, ConstraintF, GG, IG>
where
    G: Group,
    I: InjectiveMap<G>,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
    IG: InjectiveMapGadget<G, I, ConstraintF, GG>,
{
    _compressor:        PhantomData<I>,
    _compressor_gadget: PhantomData<IG>,
    _crh:               PedersenCRHGadget<G, ConstraintF, GG>,
}

impl<G, I, ConstraintF, GG, IG, W> FixedLengthCRHGadget<PedersenCRHCompressor<G, I, W>, ConstraintF>
    for PedersenCRHCompressorGadget<G, I, ConstraintF, GG, IG>
where
    G: Group,
    I: InjectiveMap<G>,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
    IG: InjectiveMapGadget<G, I, ConstraintF, GG>,
    W: PedersenWindow,
{
    type OutputGadget = IG::OutputGadget;
    type ParametersGadget = PedersenCRHGadgetParameters<G, W, ConstraintF, GG>;

    fn check_evaluation_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        parameters: &Self::ParametersGadget,
        input: &[UInt8],
    ) -> Result<Self::OutputGadget, SynthesisError> {
        let result = PedersenCRHGadget::<G, ConstraintF, GG>::check_evaluation_gadget(
            cs.ns(|| "PedCRH"),
            parameters,
            input,
        )?;
        IG::evaluate_map(cs.ns(|| "InjectiveMap"), &result)
    }
}
