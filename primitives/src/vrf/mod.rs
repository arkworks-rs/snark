use crate::Error;
use algebra::{Field, FromBytes, FromBytesChecked, SemanticallyValid, ToBytes, UniformRand};
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::{fmt::Debug, hash::Hash};

pub mod ecvrf;

pub trait FieldBasedVrf {
    type Data: Field;
    type PublicKey: FromBytes
        + FromBytesChecked
        + ToBytes
        + Hash
        + Eq
        + Copy
        + Clone
        + Default
        + Debug
        + Send
        + Sync
        + UniformRand
        + Serialize
        + for<'a> Deserialize<'a>;
    type SecretKey: ToBytes + Clone + Default + Serialize + for<'a> Deserialize<'a>;
    type Proof: Copy
        + Clone
        + Default
        + Send
        + Sync
        + Debug
        + Eq
        + PartialEq
        + ToBytes
        + FromBytes
        + FromBytesChecked
        + SemanticallyValid
        + Serialize
        + for<'a> Deserialize<'a>;
    type GHParams: Clone + Default;

    fn keygen<R: Rng>(rng: &mut R) -> (Self::PublicKey, Self::SecretKey);

    fn get_public_key(sk: &Self::SecretKey) -> Self::PublicKey;

    fn prove<R: Rng>(
        rng: &mut R,
        pp: &Self::GHParams,
        pk: &Self::PublicKey,
        sk: &Self::SecretKey,
        message: Self::Data,
    ) -> Result<Self::Proof, Error>;

    // Verifies the VRF proof and returns the VRF output
    fn proof_to_hash(
        pp: &Self::GHParams,
        pk: &Self::PublicKey,
        message: Self::Data,
        proof: &Self::Proof,
    ) -> Result<Self::Data, Error>;

    fn keyverify(pk: &Self::PublicKey) -> bool;
}
