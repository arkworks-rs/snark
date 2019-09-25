use crate::CryptoError;
use algebra::bytes::ToBytes;
use crate::Error;
use rand::Rng;
use std::{fmt::Debug, hash::Hash, marker::PhantomData};

use super::{
    pedersen::{PedersenCRH, PedersenParameters, PedersenWindow},
    FixedLengthCRH,
};
use algebra::{
    curves::{
        models::{ModelParameters, TEModelParameters},
        twisted_edwards_extended::{GroupAffine as TEAffine, GroupProjective as TEProjective},
        ProjectiveCurve,
    },
    groups::Group,
};


#[cfg(feature = "r1cs")]
pub mod constraints;

pub trait InjectiveMap<G: Group> {
    type Output: ToBytes + Clone + Eq + Hash + Default + Debug;
    fn injective_map(ge: &G) -> Result<Self::Output, CryptoError>;
}

pub struct TECompressor;

impl<P: TEModelParameters> InjectiveMap<TEAffine<P>> for TECompressor {
    type Output = <P as ModelParameters>::BaseField;

    fn injective_map(ge: &TEAffine<P>) -> Result<Self::Output, CryptoError> {
        debug_assert!(ge.is_in_correct_subgroup_assuming_on_curve());
        Ok(ge.x)
    }
}

impl<P: TEModelParameters> InjectiveMap<TEProjective<P>> for TECompressor {
    type Output = <P as ModelParameters>::BaseField;

    fn injective_map(ge: &TEProjective<P>) -> Result<Self::Output, CryptoError> {
        let ge = ge.into_affine();
        debug_assert!(ge.is_in_correct_subgroup_assuming_on_curve());
        Ok(ge.x)
    }
}

pub struct PedersenCRHCompressor<G: Group, I: InjectiveMap<G>, W: PedersenWindow> {
    _group:      PhantomData<G>,
    _compressor: PhantomData<I>,
    _crh:        PedersenCRH<G, W>,
}

impl<G: Group, I: InjectiveMap<G>, W: PedersenWindow> FixedLengthCRH
    for PedersenCRHCompressor<G, I, W>
{
    const INPUT_SIZE_BITS: usize = PedersenCRH::<G, W>::INPUT_SIZE_BITS;
    type Output = I::Output;
    type Parameters = PedersenParameters<G>;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Parameters, Error> {
        let time = start_timer!(|| format!("PedersenCRHCompressor::Setup"));
        let params = PedersenCRH::<G, W>::setup(rng);
        end_timer!(time);
        params
    }

    fn evaluate(parameters: &Self::Parameters, input: &[u8]) -> Result<Self::Output, Error> {
        let eval_time = start_timer!(|| "PedersenCRHCompressor::Eval");
        let result = I::injective_map(&PedersenCRH::<G, W>::evaluate(parameters, input)?)?;
        end_timer!(eval_time);
        Ok(result)
    }
}
