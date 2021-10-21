use crate::Error;
use algebra::{
    bytes::{FromBytes, ToBytes},
    Field, FromBytesChecked, UniformRand,
};
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::fmt::Debug;
use std::hash::Hash;

pub mod schnorr;

pub trait SignatureScheme {
    type Parameters: Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>;
    type PublicKey: ToBytes
        + Serialize
        + for<'a> Deserialize<'a>
        + Hash
        + Eq
        + Clone
        + Default
        + Send
        + Sync;
    type SecretKey: ToBytes + Serialize + for<'a> Deserialize<'a> + Clone + Default;
    type Signature: Serialize + for<'a> Deserialize<'a> + Clone + Default + Send + Sync;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Parameters, Error>;

    fn keygen<R: Rng>(
        pp: &Self::Parameters,
        rng: &mut R,
    ) -> Result<(Self::PublicKey, Self::SecretKey), Error>;

    fn sign<R: Rng>(
        pp: &Self::Parameters,
        sk: &Self::SecretKey,
        message: &[u8],
        rng: &mut R,
    ) -> Result<Self::Signature, Error>;

    fn verify(
        pp: &Self::Parameters,
        pk: &Self::PublicKey,
        message: &[u8],
        signature: &Self::Signature,
    ) -> Result<bool, Error>;

    fn randomize_public_key(
        pp: &Self::Parameters,
        public_key: &Self::PublicKey,
        randomness: &[u8],
    ) -> Result<Self::PublicKey, Error>;

    fn randomize_signature(
        pp: &Self::Parameters,
        signature: &Self::Signature,
        randomness: &[u8],
    ) -> Result<Self::Signature, Error>;
}

pub trait FieldBasedSignatureScheme {
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
    type Signature: Copy
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
        + Serialize
        + for<'a> Deserialize<'a>;

    fn keygen<R: Rng>(rng: &mut R) -> (Self::PublicKey, Self::SecretKey);

    fn get_public_key(sk: &Self::SecretKey) -> Self::PublicKey;

    fn sign<R: Rng>(
        rng: &mut R,
        pk: &Self::PublicKey,
        sk: &Self::SecretKey,
        message: Self::Data,
    ) -> Result<Self::Signature, Error>;

    fn verify(
        pk: &Self::PublicKey,
        message: Self::Data,
        signature: &Self::Signature,
    ) -> Result<bool, Error>;

    fn keyverify(pk: &Self::PublicKey) -> bool;
}
