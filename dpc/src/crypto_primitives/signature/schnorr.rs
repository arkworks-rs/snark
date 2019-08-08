use crate::crypto_primitives::SignatureScheme;
use algebra::{
    bytes::ToBytes,
    fields::{Field, PrimeField},
    groups::Group,
    to_bytes,
};
use digest::Digest;
use crate::Error;
use rand::{Rand, Rng};
use std::{
    hash::Hash,
    io::{Result as IoResult, Write},
    marker::PhantomData,
    ops::Sub,
};

pub struct SchnorrSignature<G: Group, D: Digest> {
    _group: PhantomData<G>,
    _hash:  PhantomData<D>,
}

#[derive(Derivative)]
#[derivative(Clone(bound = "G: Group, H: Digest"))]
pub struct SchnorrSigParameters<G: Group, H: Digest> {
    _hash:         PhantomData<H>,
    pub generator: G,
    pub salt:      [u8; 32],
}

pub type SchnorrPublicKey<G> = G;

#[derive(Derivative)]
#[derivative(Clone(bound = "G: Group"), Default(bound = "G: Group"))]
pub struct SchnorrSecretKey<G: Group>(pub G::ScalarField);

impl<G: Group> ToBytes for SchnorrSecretKey<G> {
    #[inline]
    fn write<W: Write>(&self, writer: W) -> IoResult<()> {
        self.0.write(writer)
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = "G: Group"), Default(bound = "G: Group"))]
pub struct SchnorrSig<G: Group> {
    pub prover_response:    G::ScalarField,
    pub verifier_challenge: G::ScalarField,
}

impl<G: Group + Hash, D: Digest + Send + Sync> SignatureScheme for SchnorrSignature<G, D>
where
    G::ScalarField: PrimeField,
{
    type Parameters = SchnorrSigParameters<G, D>;
    type PublicKey = G;
    type SecretKey = SchnorrSecretKey<G>;
    type Signature = SchnorrSig<G>;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Parameters, Error> {
        let setup_time = timer_start!(|| "SchnorrSig::Setup");

        let mut salt = [0u8; 32];
        rng.fill_bytes(&mut salt);
        let generator = G::rand(rng);

        timer_end!(setup_time);
        Ok(SchnorrSigParameters {
            _hash: PhantomData,
            generator,
            salt,
        })
    }

    fn keygen<R: Rng>(
        parameters: &Self::Parameters,
        rng: &mut R,
    ) -> Result<(Self::PublicKey, Self::SecretKey), Error> {
        let keygen_time = timer_start!(|| "SchnorrSig::KeyGen");

        let secret_key = G::ScalarField::rand(rng);
        let public_key = parameters.generator.mul(&secret_key);

        timer_end!(keygen_time);
        Ok((public_key, SchnorrSecretKey(secret_key)))
    }

    fn sign<R: Rng>(
        parameters: &Self::Parameters,
        sk: &Self::SecretKey,
        message: &[u8],
        rng: &mut R,
    ) -> Result<Self::Signature, Error> {
        let sign_time = timer_start!(|| "SchnorrSig::Sign");
        // (k, e);
        let (random_scalar, verifier_challenge) = loop {
            // Sample a random scalar `k` from the prime scalar field.
            let random_scalar: G::ScalarField = G::ScalarField::rand(rng);
            // Commit to the random scalar via r := k · g.
            // This is the prover's first msg in the Sigma protocol.
            let prover_commitment: G = parameters.generator.mul(&random_scalar);

            // Hash everything to get verifier challenge.
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(&parameters.salt);
            hash_input.extend_from_slice(&to_bytes![prover_commitment]?);
            hash_input.extend_from_slice(message);

            // Compute the supposed verifier response: e := H(salt || r || msg);
            if let Some(verifier_challenge) =
                G::ScalarField::from_random_bytes(&D::digest(&hash_input))
            {
                break (random_scalar, verifier_challenge);
            };
        };

        // k - xe;
        let prover_response = random_scalar - &(verifier_challenge * &sk.0);
        let signature = SchnorrSig {
            prover_response,
            verifier_challenge,
        };

        timer_end!(sign_time);
        Ok(signature)
    }

    fn verify(
        parameters: &Self::Parameters,
        pk: &Self::PublicKey,
        message: &[u8],
        signature: &Self::Signature,
    ) -> Result<bool, Error> {
        let verify_time = timer_start!(|| "SchnorrSig::Verify");

        let &SchnorrSig {
            ref prover_response,
            ref verifier_challenge,
        } = signature;
        let mut claimed_prover_commitment = parameters.generator.mul(prover_response);
        let public_key_times_verifier_challenge = pk.mul(verifier_challenge);
        claimed_prover_commitment += &public_key_times_verifier_challenge;

        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(&parameters.salt);
        hash_input.extend_from_slice(&to_bytes![claimed_prover_commitment]?);
        hash_input.extend_from_slice(&message);

        let obtained_verifier_challenge = if let Some(obtained_verifier_challenge) =
            G::ScalarField::from_random_bytes(&D::digest(&hash_input))
        {
            obtained_verifier_challenge
        } else {
            return Ok(false);
        };
        timer_end!(verify_time);
        Ok(verifier_challenge == &obtained_verifier_challenge)
    }

    fn randomize_public_key(
        parameters: &Self::Parameters,
        public_key: &Self::PublicKey,
        randomness: &[u8],
    ) -> Result<Self::PublicKey, Error> {
        let rand_pk_time = timer_start!(|| "SchnorrSig::RandomizePubKey");

        let mut randomized_pk = *public_key;
        let mut base = parameters.generator;
        let mut encoded = G::zero();
        for bit in bytes_to_bits(randomness) {
            if bit {
                encoded += &base;
            }
            base.double_in_place();
        }
        randomized_pk += &encoded;

        timer_end!(rand_pk_time);

        Ok(randomized_pk)
    }

    fn randomize_signature(
        _parameter: &Self::Parameters,
        signature: &Self::Signature,
        randomness: &[u8],
    ) -> Result<Self::Signature, Error> {
        let rand_signature_time = timer_start!(|| "SchnorrSig::RandomizeSig");
        let &SchnorrSig {
            ref prover_response,
            ref verifier_challenge,
        } = signature;
        let mut base = G::ScalarField::one();
        let mut multiplier = G::ScalarField::zero();
        for bit in bytes_to_bits(randomness) {
            if bit {
                multiplier += &base;
            }
            base.double_in_place();
        }

        let new_sig = SchnorrSig {
            prover_response:    prover_response.sub(&(*verifier_challenge * &multiplier)),
            verifier_challenge: *verifier_challenge,
        };
        timer_end!(rand_signature_time);
        Ok(new_sig)
    }
}

pub fn bytes_to_bits(bytes: &[u8]) -> Vec<bool> {
    let mut bits = Vec::with_capacity(bytes.len() * 8);
    for byte in bytes {
        for i in 0..8 {
            let bit = (*byte >> (8 - i - 1)) & 1;
            bits.push(bit == 1);
        }
    }
    bits
}
