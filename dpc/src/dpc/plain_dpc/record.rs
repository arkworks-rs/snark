use crypto_primitives::{CommitmentScheme, FixedLengthCRH, PRF};
use crate::{
    dpc::{
        plain_dpc::{AddressPublicKey, DPCPredicate, PlainDPCComponents},
        Record,
    },
};
use algebra::to_bytes;
use std::marker::PhantomData;

#[derive(Derivative)]
#[derivative(
    Default(bound = "C: PlainDPCComponents"),
    Clone(bound = "C: PlainDPCComponents")
)]
pub struct DPCRecord<C: PlainDPCComponents> {
    pub(super) address_public_key: AddressPublicKey<C>,

    pub(super) is_dummy: bool,
    pub(super) payload:  [u8; 32],

    #[derivative(Default(value = "default_predicate_hash::<C::PredVkH>()"))]
    pub(super) birth_predicate_repr: Vec<u8>,
    #[derivative(Default(value = "default_predicate_hash::<C::PredVkH>()"))]
    pub(super) death_predicate_repr: Vec<u8>,

    pub(super) serial_number_nonce: <C::SnNonceH as FixedLengthCRH>::Output,

    pub(super) commitment:            <C::RecC as CommitmentScheme>::Output,
    pub(super) commitment_randomness: <C::RecC as CommitmentScheme>::Randomness,

    pub(super) _components: PhantomData<C>,
}

fn default_predicate_hash<C: FixedLengthCRH>() -> Vec<u8> {
    use algebra::bytes::ToBytes;
    to_bytes![C::Output::default()].unwrap()
}

impl<C: PlainDPCComponents> Record for DPCRecord<C> {
    type AddressPublicKey = AddressPublicKey<C>;
    type Commitment = <C::RecC as CommitmentScheme>::Output;
    type CommitmentRandomness = <C::RecC as CommitmentScheme>::Randomness;

    type Payload = [u8; 32];
    type Predicate = DPCPredicate<C>;
    type SerialNumberNonce = <C::SnNonceH as FixedLengthCRH>::Output;
    type SerialNumber = <C::P as PRF>::Output;

    fn address_public_key(&self) -> &Self::AddressPublicKey {
        &self.address_public_key
    }

    fn is_dummy(&self) -> bool {
        self.is_dummy
    }

    fn payload(&self) -> &Self::Payload {
        &self.payload
    }

    fn birth_predicate_repr(&self) -> &[u8] {
        &self.birth_predicate_repr
    }

    fn death_predicate_repr(&self) -> &[u8] {
        &self.death_predicate_repr
    }

    fn serial_number_nonce(&self) -> &Self::SerialNumberNonce {
        &self.serial_number_nonce
    }

    fn commitment(&self) -> Self::Commitment {
        self.commitment.clone()
    }

    fn commitment_randomness(&self) -> Self::CommitmentRandomness {
        self.commitment_randomness.clone()
    }
}
