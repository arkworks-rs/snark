use crate::commitment::{
    injective_map::{InjectiveMap, PedersenCommCompressor},
    pedersen::{
        constraints::{CommGadget, ParametersVar, RandomnessVar},
        Window,
    },
};

pub use crate::crh::injective_map::constraints::InjectiveMapGadget;
use algebra_core::{Field, PrimeField, ProjectiveCurve};
use r1cs_core::SynthesisError;
use r1cs_std::{
    groups::{CurveVar, GroupOpsBounds},
    uint8::UInt8,
};

use core::marker::PhantomData;

type ConstraintF<C> = <<C as ProjectiveCurve>::BaseField as Field>::BasePrimeField;

pub struct CommitmentCompressorGadget<C, I, W, GG, IG>
where
    C: ProjectiveCurve,
    I: InjectiveMap<C>,
    W: Window,
    GG: CurveVar<C, ConstraintF<C>>,
    IG: InjectiveMapGadget<C, I, GG>,
    for<'a> &'a GG: GroupOpsBounds<'a, C, GG>,
{
    _compressor: PhantomData<I>,
    _compressor_gadget: PhantomData<IG>,
    _comm: PhantomData<CommGadget<C, GG, W>>,
}

impl<C, I, GG, IG, W>
    crate::commitment::CommitmentGadget<PedersenCommCompressor<C, I, W>, ConstraintF<C>>
    for CommitmentCompressorGadget<C, I, W, GG, IG>
where
    C: ProjectiveCurve,
    I: InjectiveMap<C>,
    GG: CurveVar<C, ConstraintF<C>>,
    ConstraintF<C>: PrimeField,
    IG: InjectiveMapGadget<C, I, GG>,
    W: Window,
    for<'a> &'a GG: GroupOpsBounds<'a, C, GG>,
{
    type OutputVar = IG::OutputVar;
    type ParametersVar = ParametersVar<C, GG>;
    type RandomnessVar = RandomnessVar<ConstraintF<C>>;

    fn commit(
        parameters: &Self::ParametersVar,
        input: &[UInt8<ConstraintF<C>>],
        r: &Self::RandomnessVar,
    ) -> Result<Self::OutputVar, SynthesisError> {
        let result = CommGadget::<C, GG, W>::commit(parameters, input, r)?;
        IG::evaluate(&result)
    }
}
