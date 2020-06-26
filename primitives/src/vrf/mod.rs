use algebra::{Field, ToBytes, FromBytes, FromBytesChecked, SemanticallyValid, UniformRand};
use rand::Rng;
use std::{hash::Hash, fmt::Debug};
use crate::Error;

pub mod ecvrf;

pub trait FieldBasedVrf {
    type Data: Field;
    type PublicKey: FromBytes + FromBytesChecked + ToBytes + Hash + Eq + Clone +
                    Default + Debug + Send + Sync + UniformRand;
    type SecretKey: ToBytes + Clone + Default;
    type Proof: Copy + Clone + Default + Send + Sync + Debug + Eq + PartialEq + ToBytes
                + FromBytes + FromBytesChecked + SemanticallyValid;
    type GHParams: Clone + Default;

    fn keygen<R: Rng>(
        rng: &mut R,
    ) -> (Self::PublicKey, Self::SecretKey);

    fn get_public_key(
        sk: &Self::SecretKey
    ) -> Self::PublicKey;

    fn prove<R: Rng>
    (
        rng:     &mut R,
        pp:      &Self::GHParams,
        pk:      &Self::PublicKey,
        sk:      &Self::SecretKey,
        message: &[Self::Data],
    ) -> Result<Self::Proof, Error>;

    // Verifies the VRF proof and returns the VRF output
    fn proof_to_hash
    (
        pp:      &Self::GHParams,
        pk:      &Self::PublicKey,
        message: &[Self::Data],
        proof:   &Self::Proof,
    ) -> Result<Self::Data, Error>;

    fn keyverify(pk: &Self::PublicKey) -> bool;
}