use algebra::{groups::Group, Field};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use crate::signature::SigRandomizePkGadget;

use std::{borrow::Borrow, marker::PhantomData};

use digest::Digest;
use primitives::signature::schnorr::{SchnorrPublicKey, SchnorrSigParameters, SchnorrSignature};

pub mod field_based_schnorr;

pub struct SchnorrSigGadgetParameters<G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>>
{
    generator: GG,
    _group: PhantomData<*const G>,
    _engine: PhantomData<*const ConstraintF>,
}

impl<G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>> Clone
    for SchnorrSigGadgetParameters<G, ConstraintF, GG>
{
    fn clone(&self) -> Self {
        Self {
            generator: self.generator.clone(),
            _group: PhantomData,
            _engine: PhantomData,
        }
    }
}

#[derive(Derivative)]
#[derivative(
    Debug(bound = "G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>"),
    Clone(bound = "G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>"),
    PartialEq(bound = "G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>"),
    Eq(bound = "G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>")
)]
pub struct SchnorrSigGadgetPk<G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>> {
    pub_key: GG,
    #[doc(hidden)]
    _group: PhantomData<*const G>,
    #[doc(hidden)]
    _engine: PhantomData<*const ConstraintF>,
}

pub struct SchnorrRandomizePkGadget<G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>> {
    #[doc(hidden)]
    _group: PhantomData<*const G>,
    #[doc(hidden)]
    _group_gadget: PhantomData<*const GG>,
    #[doc(hidden)]
    _engine: PhantomData<*const ConstraintF>,
}

impl<G, GG, D, ConstraintF> SigRandomizePkGadget<SchnorrSignature<G, D>, ConstraintF>
    for SchnorrRandomizePkGadget<G, ConstraintF, GG>
where
    G: Group,
    GG: GroupGadget<G, ConstraintF>,
    D: Digest + Send + Sync,
    ConstraintF: Field,
{
    type ParametersGadget = SchnorrSigGadgetParameters<G, ConstraintF, GG>;
    type PublicKeyGadget = SchnorrSigGadgetPk<G, ConstraintF, GG>;

    fn check_randomization_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        parameters: &Self::ParametersGadget,
        public_key: &Self::PublicKeyGadget,
        randomness: &[UInt8],
    ) -> Result<Self::PublicKeyGadget, SynthesisError> {
        let base = parameters.generator.clone();
        let randomness = randomness
            .iter()
            .flat_map(|b| b.into_bits_le())
            .collect::<Vec<_>>();
        let rand_pk = base.mul_bits(
            &mut cs.ns(|| "Compute Randomizer"),
            &public_key.pub_key,
            randomness.iter(),
        )?;
        Ok(SchnorrSigGadgetPk {
            pub_key: rand_pk,
            _group: PhantomData,
            _engine: PhantomData,
        })
    }
}

impl<G, ConstraintF, GG, D> AllocGadget<SchnorrSigParameters<G, D>, ConstraintF>
    for SchnorrSigGadgetParameters<G, ConstraintF, GG>
where
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
    D: Digest,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SchnorrSigParameters<G, D>>,
    {
        let generator = GG::alloc_checked(cs, || f().map(|pp| pp.borrow().generator))?;
        Ok(Self {
            generator,
            _engine: PhantomData,
            _group: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SchnorrSigParameters<G, D>>,
    {
        let generator = GG::alloc_input(cs, || f().map(|pp| pp.borrow().generator))?;
        Ok(Self {
            generator,
            _engine: PhantomData,
            _group: PhantomData,
        })
    }
}

impl<G, ConstraintF, GG> AllocGadget<SchnorrPublicKey<G>, ConstraintF>
    for SchnorrSigGadgetPk<G, ConstraintF, GG>
where
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SchnorrPublicKey<G>>,
    {
        let pub_key = GG::alloc_input(cs, || f().map(|pk| *pk.borrow()))?;
        Ok(Self {
            pub_key,
            _engine: PhantomData,
            _group: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SchnorrPublicKey<G>>,
    {
        let pub_key = GG::alloc_input(cs, || f().map(|pk| *pk.borrow()))?;
        Ok(Self {
            pub_key,
            _engine: PhantomData,
            _group: PhantomData,
        })
    }
}

impl<G, ConstraintF, GG> EqGadget<ConstraintF> for SchnorrSigGadgetPk<G, ConstraintF, GG>
where
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
{
    #[inline]
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        self.pub_key.is_eq(cs, &other.pub_key)
    }

    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.pub_key
            .conditional_enforce_equal(cs, &other.pub_key, condition)
    }

    #[inline]
    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.pub_key
            .conditional_enforce_not_equal(cs, &other.pub_key, condition)
    }
}

impl<G, ConstraintF, GG> ToBytesGadget<ConstraintF> for SchnorrSigGadgetPk<G, ConstraintF, GG>
where
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
{
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.pub_key.to_bytes(&mut cs.ns(|| "PubKey To Bytes"))
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.pub_key
            .to_bytes_strict(&mut cs.ns(|| "PubKey To Bytes"))
    }
}
