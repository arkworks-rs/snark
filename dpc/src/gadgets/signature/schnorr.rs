use algebra::{groups::Group, PairingEngine};
use snark::{ConstraintSystem, SynthesisError};

use crate::gadgets::signature::SigRandomizePkGadget;
use snark_gadgets::{
    boolean::Boolean,
    groups::GroupGadget,
    uint8::UInt8,
    utils::{AllocGadget, ConditionalEqGadget, EqGadget, ToBytesGadget},
};

use std::{borrow::Borrow, marker::PhantomData};

use crate::crypto_primitives::signature::schnorr::{
    SchnorrPublicKey, SchnorrSigParameters, SchnorrSignature,
};
use digest::Digest;

pub struct SchnorrSigGadgetParameters<G: Group, E: PairingEngine, GG: GroupGadget<G, E>> {
    generator: GG,
    _group:    PhantomData<*const G>,
    _engine:   PhantomData<*const E>,
}

impl<G: Group, E: PairingEngine, GG: GroupGadget<G, E>> Clone
    for SchnorrSigGadgetParameters<G, E, GG>
{
    fn clone(&self) -> Self {
        Self {
            generator: self.generator.clone(),
            _group:    PhantomData,
            _engine:   PhantomData,
        }
    }
}

#[derive(Derivative)]
#[derivative(
    Debug(bound = "G: Group, E: PairingEngine, GG: GroupGadget<G, E>"),
    Clone(bound = "G: Group, E: PairingEngine, GG: GroupGadget<G, E>"),
    PartialEq(bound = "G: Group, E: PairingEngine, GG: GroupGadget<G, E>"),
    Eq(bound = "G: Group, E: PairingEngine, GG: GroupGadget<G, E>")
)]
pub struct SchnorrSigGadgetPk<G: Group, E: PairingEngine, GG: GroupGadget<G, E>> {
    pub_key: GG,
    _group:  PhantomData<*const G>,
    _engine: PhantomData<*const E>,
}

pub struct SchnorrRandomizePkGadget<G: Group, E: PairingEngine, GG: GroupGadget<G, E>> {
    _group:        PhantomData<*const G>,
    _group_gadget: PhantomData<*const GG>,
    _engine:       PhantomData<*const E>,
}

impl<G, GG, D, E> SigRandomizePkGadget<SchnorrSignature<G, D>, E>
    for SchnorrRandomizePkGadget<G, E, GG>
where
    G: Group,
    GG: GroupGadget<G, E>,
    D: Digest + Send + Sync,
    E: PairingEngine,
{
    type ParametersGadget = SchnorrSigGadgetParameters<G, E, GG>;
    type PublicKeyGadget = SchnorrSigGadgetPk<G, E, GG>;

    fn check_randomization_gadget<CS: ConstraintSystem<E>>(
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
            _group:  PhantomData,
            _engine: PhantomData,
        })
    }
}

impl<G, E, GG, D> AllocGadget<SchnorrSigParameters<G, D>, E>
    for SchnorrSigGadgetParameters<G, E, GG>
where
    G: Group,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
    D: Digest,
{
    fn alloc<F, T, CS: ConstraintSystem<E>>(cs: CS, f: F) -> Result<Self, SynthesisError>
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

    fn alloc_input<F, T, CS: ConstraintSystem<E>>(cs: CS, f: F) -> Result<Self, SynthesisError>
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

impl<G, E, GG> AllocGadget<SchnorrPublicKey<G>, E> for SchnorrSigGadgetPk<G, E, GG>
where
    G: Group,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
{
    fn alloc<F, T, CS: ConstraintSystem<E>>(cs: CS, f: F) -> Result<Self, SynthesisError>
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

    fn alloc_input<F, T, CS: ConstraintSystem<E>>(cs: CS, f: F) -> Result<Self, SynthesisError>
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

impl<G, E, GG> ConditionalEqGadget<E> for SchnorrSigGadgetPk<G, E, GG>
where
    G: Group,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
{
    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.pub_key.conditional_enforce_equal(
            &mut cs.ns(|| "PubKey equality"),
            &other.pub_key,
            condition,
        )?;
        Ok(())
    }

    fn cost() -> usize {
        <GG as ConditionalEqGadget<E>>::cost()
    }
}

impl<G, E, GG> EqGadget<E> for SchnorrSigGadgetPk<G, E, GG>
where
    G: Group,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
{
}

impl<G, E, GG> ToBytesGadget<E> for SchnorrSigGadgetPk<G, E, GG>
where
    G: Group,
    E: PairingEngine,
    GG: GroupGadget<G, E>,
{
    fn to_bytes<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        self.pub_key.to_bytes(&mut cs.ns(|| "PubKey To Bytes"))
    }

    fn to_bytes_strict<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.pub_key
            .to_bytes_strict(&mut cs.ns(|| "PubKey To Bytes"))
    }
}
