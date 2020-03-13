use crate::{Error, SignatureScheme};
use algebra::{
    bytes::ToBytes,
    fields::{Field, PrimeField},
    groups::Group,
    to_bytes, ToConstraintField, UniformRand,
};
use digest::Digest;
use rand::Rng;
use std::{
    hash::Hash,
    io::{Result as IoResult, Write},
    marker::PhantomData,
};

#[cfg(feature = "r1cs")]
pub mod constraints;

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
        let setup_time = start_timer!(|| "SchnorrSig::Setup");

        let mut salt = [0u8; 32];
        rng.fill_bytes(&mut salt);
        let generator = G::rand(rng);

        end_timer!(setup_time);
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
        let keygen_time = start_timer!(|| "SchnorrSig::KeyGen");

        let secret_key = G::ScalarField::rand(rng);
        let public_key = parameters.generator.mul(&secret_key);

        end_timer!(keygen_time);
        Ok((public_key, SchnorrSecretKey(secret_key)))
    }

    fn sign<R: Rng>(
        parameters: &Self::Parameters,
        sk: &Self::SecretKey,
        message: &[u8],
        rng: &mut R,
    ) -> Result<Self::Signature, Error> {
        let sign_time = start_timer!(|| "SchnorrSig::Sign");
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

        end_timer!(sign_time);
        Ok(signature)
    }

    fn verify(
        parameters: &Self::Parameters,
        pk: &Self::PublicKey,
        message: &[u8],
        signature: &Self::Signature,
    ) -> Result<bool, Error> {
        let verify_time = start_timer!(|| "SchnorrSig::Verify");

        let SchnorrSig {
            prover_response,
            verifier_challenge,
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
        end_timer!(verify_time);
        Ok(verifier_challenge == &obtained_verifier_challenge)
    }

    fn randomize_public_key(
        parameters: &Self::Parameters,
        public_key: &Self::PublicKey,
        randomness: &[u8],
    ) -> Result<Self::PublicKey, Error> {
        let rand_pk_time = start_timer!(|| "SchnorrSig::RandomizePubKey");

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

        end_timer!(rand_pk_time);

        Ok(randomized_pk)
    }

    fn randomize_signature(
        _parameter: &Self::Parameters,
        signature: &Self::Signature,
        randomness: &[u8],
    ) -> Result<Self::Signature, Error> {
        let rand_signature_time = start_timer!(|| "SchnorrSig::RandomizeSig");
        let SchnorrSig {
            prover_response,
            verifier_challenge,
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
            prover_response:    *prover_response - &(*verifier_challenge * &multiplier),
            verifier_challenge: *verifier_challenge,
        };
        end_timer!(rand_signature_time);
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

impl<ConstraintF: Field, G: Group + ToConstraintField<ConstraintF>, D: Digest>
ToConstraintField<ConstraintF> for SchnorrSigParameters<G, D>
{
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        self.generator.to_field_elements()
    }
}

pub mod field_impl {

    use crate::{crh::FieldBasedHash, signature::FieldBasedSignatureScheme, CryptoError, Error, compute_truncation_size};
    use algebra::{Field, PrimeField, Group, UniformRand, ProjectiveCurve,
                  convert, leading_zeros, ToBits, ToConstraintField, ToBytes, FromBytes};
    use std::marker::PhantomData;
    use rand::Rng;
    use std::io::{Write, Read, Result as IoResult};

    #[allow(dead_code)]
    pub struct FieldBasedSchnorrSignatureScheme<
        F: PrimeField,
        G: Group,
        H: FieldBasedHash,
    >
    {
        _field:    PhantomData<F>,
        _group:    PhantomData<G>,
        _hash:     PhantomData<H>,
    }

    #[derive(Derivative)]
    #[derivative(
    Copy(bound = "F: PrimeField"),
    Clone(bound = "F: PrimeField"),
    Default(bound = "F: PrimeField"),
    Eq(bound = "F: PrimeField"),
    PartialEq(bound = "F: PrimeField"),
    Debug(bound = "F: PrimeField")
    )]
    pub struct FieldBasedSchnorrSignature<F: PrimeField> {
        pub e:    F,
        pub s:    F,
    }

    impl<F: PrimeField> ToBytes for FieldBasedSchnorrSignature<F> {
        fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
            self.e.write(&mut writer)?;
            self.s.write(&mut writer)
        }
    }

    impl<F: PrimeField> FromBytes for FieldBasedSchnorrSignature<F> {
        fn read<R: Read>(mut reader: R) -> IoResult<Self> {
            let e = F::read(&mut reader)?;
            let s = F::read(&mut reader)?;
            Ok(Self{ e, s })
        }
    }

    impl<F: PrimeField, G: ProjectiveCurve + ToConstraintField<F>, H: FieldBasedHash<Data = F>> FieldBasedSignatureScheme for
    FieldBasedSchnorrSignatureScheme<F, G, H>
    {
        type Data = H::Data;
        type PublicKey = G;
        type SecretKey = G::ScalarField;
        type Signature = FieldBasedSchnorrSignature<F>;

        fn keygen<R: Rng>(rng: &mut R) -> (Self::PublicKey, Self::SecretKey)
        {
            let secret_key = G::ScalarField::rand(rng);
            let public_key = G::prime_subgroup_generator()
                .mul(&secret_key);
            (public_key, secret_key)
        }

        fn sign<R: Rng>(
            rng: &mut R,
            pk: &Self::PublicKey,
            sk: &Self::SecretKey,
            message: &[Self::Data],
        )-> Result<Self::Signature, Error>
        {
            let pk_coords = pk.to_field_elements()?;

            let (e, s) = loop {

                //Sample random element
                let k = G::ScalarField::rand(rng);

                if k.is_zero() {continue};

                //R = k * G
                let r = G::prime_subgroup_generator()
                    .mul(&k);

                let r_coords = r.to_field_elements()?;

                // Compute e = H(m || R || pk.x)
                let mut hash_input = Vec::new();
                hash_input.extend_from_slice(message);
                hash_input.extend_from_slice(r_coords.as_slice());
                hash_input.push(pk_coords[0]);
                let e = H::evaluate(hash_input.as_ref())?;

                let e_leading_zeros = leading_zeros(e.write_bits()) as usize;
                let required_leading_zeros = compute_truncation_size(
                    F::size_in_bits() as i32,
                    G::ScalarField::size_in_bits() as i32,
                );

                //Enforce e bit length is strictly smaller than G::ScalarField modulus bit length
                if e_leading_zeros < required_leading_zeros {continue};

                //We can now safely convert it to the other field
                let e_conv = convert::<F, G::ScalarField>(e)?;

                //Enforce s bit length is strictly smaller than F modulus bit length
                let s = k + &(e_conv * sk);

                let s_leading_zeros = leading_zeros(s.write_bits()) as usize;
                let required_leading_zeros = compute_truncation_size(
                    G::ScalarField::size_in_bits() as i32,
                    F::size_in_bits() as i32,
                );

                if s_leading_zeros < required_leading_zeros {continue};

                let s_conv = convert::<G::ScalarField, F>(s)?;

                break (e, s_conv);
            };

            Ok(FieldBasedSchnorrSignature {e, s})
        }

        fn verify(
            pk: &Self::PublicKey,
            message: &[Self::Data],
            signature: &Self::Signature
        )
            -> Result<bool, Error>
        {

            let pk_coords = pk.to_field_elements()?;

            //Checks
            let e_bits = signature.e.write_bits();
            let e_leading_zeros = leading_zeros(e_bits.clone()) as usize;
            if (F::size_in_bits() - e_leading_zeros) >= G::ScalarField::size_in_bits(){
                return Err(Box::new(CryptoError::IncorrectInputLength("signature.e".to_owned(), e_bits.len() - e_leading_zeros)))
            }

            let s_bits = signature.s.write_bits();
            let s_leading_zeros = leading_zeros(s_bits.clone()) as usize;
            if (G::ScalarField::size_in_bits() - s_leading_zeros) >= F::size_in_bits(){
                return Err(Box::new(CryptoError::IncorrectInputLength("signature.s".to_owned(), s_bits.len() - s_leading_zeros)))
            }

            //Compute R' = s*G - e * pk
            let r_prime = {
                let s_conv = convert::<F, G::ScalarField>(signature.s)?;
                let e_conv =  convert::<F, G::ScalarField>(signature.e)?;
                let s_times_g = G::prime_subgroup_generator().mul(&s_conv);
                let neg_e_times_pk = pk.neg().mul(&e_conv);
                (s_times_g + &neg_e_times_pk)
            };

            let r_prime_coords = r_prime.to_field_elements()?;

            // Compute e' = H(m || R' || pk.x)
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(message);
            hash_input.extend_from_slice(r_prime_coords.as_slice());
            hash_input.push(pk_coords[0]);
            let e_prime = H::evaluate(hash_input.as_ref())?;

            Ok(signature.e == e_prime)
        }

        fn keyverify(pk: &Self::PublicKey) -> bool
        {
            pk.group_membership_test()
        }
    }
}

#[cfg(test)]
mod test {
    use algebra::curves::{
        mnt4753::G1Projective as MNT4G1Projective,
        mnt6753::G1Projective as MNT6G1Projective,
        jubjub::JubJubProjective,
    };
    use algebra::fields::{
        mnt4753::Fr as MNT4Fr,
        mnt6753::Fr as MNT6Fr,
        bls12_381::Fr as BLS12Fr,
    };
    use algebra::{ToBytes, to_bytes, FromBytes};
    use crate::crh::{MNT4PoseidonHash, MNT6PoseidonHash, BLS12PoseidonHash};
    use crate::signature::FieldBasedSignatureScheme;
    use crate::signature::schnorr::field_impl::FieldBasedSchnorrSignatureScheme;
    use rand::{Rng, thread_rng};

    type SchnorrMNT4 = FieldBasedSchnorrSignatureScheme<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash>;
    type SchnorrMNT6 = FieldBasedSchnorrSignatureScheme<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash>;
    type SchnorrBls12 = FieldBasedSchnorrSignatureScheme<BLS12Fr, JubJubProjective, BLS12PoseidonHash>;

    fn sign_and_verify<S: FieldBasedSignatureScheme, R: Rng>(rng: &mut R, message: &[S::Data]) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));
        let sig = S::sign(rng, &pk, &sk, &message).unwrap();
        assert!(S::verify(&pk, &message, &sig).unwrap());

        //Serialization/deserialization test
        let sig_serialized = to_bytes!(sig).unwrap();
        let sig_deserialized = <S as FieldBasedSignatureScheme>::Signature::read(sig_serialized.as_slice()).unwrap();
        assert_eq!(sig, sig_deserialized);
        assert!(S::verify(&pk, &message, &sig_deserialized).unwrap());
    }

    fn failed_verification<S: FieldBasedSignatureScheme, R: Rng>(rng: &mut R, message: &[S::Data], bad_message: &[S::Data]) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));

        //Attempt to verify a signature for a different message
        let sig = S::sign(rng, &pk, &sk, message).unwrap();
        assert!(!S::verify(&pk, bad_message, &sig).unwrap());

        //Attempt to verify a different signature for a message
        let bad_sig = S::sign(rng, &pk, &sk, bad_message).unwrap();
        assert!(!S::verify(&pk, message, &bad_sig).unwrap());

        //Attempt to verify a signature for a message but with different public key
        let (new_pk, _) = S::keygen(rng);
        assert!(!S::verify(&new_pk, message, &sig).unwrap());
    }

    #[test]
    fn mnt4_schnorr_test() {
        let rng = &mut thread_rng();
        let samples = 100;
        for _ in 0..samples {
            let f: MNT4Fr = rng.gen();
            let g: MNT4Fr = rng.gen();
            sign_and_verify::<SchnorrMNT4, _>(rng, &[f, g]);
            failed_verification::<SchnorrMNT4, _>(rng, &[f], &[g]);
        }
    }

    #[test]
    fn mnt6_schnorr_test() {
        let rng = &mut thread_rng();
        let samples = 100;
        for _ in 0..samples{
            let f: MNT6Fr = rng.gen();
            let g: MNT6Fr = rng.gen();
            sign_and_verify::<SchnorrMNT6, _>(rng,&[f, g]);
            failed_verification::<SchnorrMNT6, _>(rng, &[f], &[g]);
        }
    }

    #[test]
    fn bls12_381_schnorr_test() {
        let rng = &mut thread_rng();
        let samples = 100;
        for _ in 0..samples{
            let f: BLS12Fr = rng.gen();
            let g: BLS12Fr = rng.gen();
            sign_and_verify::<SchnorrBls12, _>(rng,&[f, g]);
            failed_verification::<SchnorrBls12, _>(rng, &[f], &[g]);
        }
    }
}