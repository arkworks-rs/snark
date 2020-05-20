use crate::dpc::{delegable_dpc::DelegableDPCComponents, Predicate};
use crypto_primitives::NIZK;
use std::marker::PhantomData;

pub struct PrivatePredInput<C: DelegableDPCComponents> {
    pub vk: <C::PredicateNIZK as NIZK>::VerificationParameters,
    pub proof: <C::PredicateNIZK as NIZK>::Proof,
}

impl<C: DelegableDPCComponents> Default for PrivatePredInput<C> {
    fn default() -> Self {
        Self {
            vk: <C::PredicateNIZK as NIZK>::VerificationParameters::default(),
            proof: <C::PredicateNIZK as NIZK>::Proof::default(),
        }
    }
}

impl<C: DelegableDPCComponents> Clone for PrivatePredInput<C> {
    fn clone(&self) -> Self {
        Self {
            vk: self.vk.clone(),
            proof: self.proof.clone(),
        }
    }
}

#[derive(Derivative)]
#[derivative(
    Clone(bound = "C: DelegableDPCComponents"),
    Default(bound = "C: DelegableDPCComponents")
)]
pub struct DPCPredicate<C: DelegableDPCComponents> {
    #[derivative(Default(value = "vec![0u8; 32]"))]
    identity: Vec<u8>,
    _components: PhantomData<C>,
}

impl<C: DelegableDPCComponents> DPCPredicate<C> {
    pub fn new(identity: Vec<u8>) -> Self {
        Self {
            identity,
            _components: PhantomData,
        }
    }
}

impl<C: DelegableDPCComponents> Predicate for DPCPredicate<C> {
    type PublicInput = ();
    type PrivateWitness = PrivatePredInput<C>;

    fn evaluate(&self, _p: &Self::PublicInput, _w: &Self::PrivateWitness) -> bool {
        unimplemented!()
    }

    fn into_compact_repr(&self) -> Vec<u8> {
        self.identity.clone()
    }
}
