use crate::dpc::{delegable_dpc::DelegableDPCComponents, AddressKeyPair};
use algebra::bytes::ToBytes;
use crypto_primitives::{CommitmentScheme, SignatureScheme, PRF};
use std::io::{Result as IoResult, Write};

#[derive(Derivative)]
#[derivative(
    Default(bound = "C: DelegableDPCComponents"),
    Clone(bound = "C: DelegableDPCComponents")
)]
pub struct AddressPublicKey<C: DelegableDPCComponents> {
    pub public_key: <C::AddrC as CommitmentScheme>::Output,
}

impl<C: DelegableDPCComponents> ToBytes for AddressPublicKey<C> {
    fn write<W: Write>(&self, writer: W) -> IoResult<()> {
        self.public_key.write(writer)
    }
}

#[derive(Derivative)]
#[derivative(
    Default(bound = "C: DelegableDPCComponents"),
    Clone(bound = "C: DelegableDPCComponents")
)]
pub struct AddressSecretKey<C: DelegableDPCComponents> {
    pub pk_sig: <C::S as SignatureScheme>::PublicKey,
    pub sk_sig: <C::S as SignatureScheme>::SecretKey,
    pub sk_prf: <C::P as PRF>::Seed,
    pub metadata: [u8; 32],
    pub r_pk: <C::AddrC as CommitmentScheme>::Randomness,
}

#[derive(Derivative)]
#[derivative(Clone(bound = "C: DelegableDPCComponents"))]
pub struct AddressPair<C: DelegableDPCComponents> {
    pub public_key: AddressPublicKey<C>,
    pub secret_key: AddressSecretKey<C>,
}

impl<C: DelegableDPCComponents> AddressKeyPair for AddressPair<C> {
    type AddressSecretKey = AddressSecretKey<C>;
    type AddressPublicKey = AddressPublicKey<C>;
}
