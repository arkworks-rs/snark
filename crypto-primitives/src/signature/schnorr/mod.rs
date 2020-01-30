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

//TODO: Replace to_bytes with to_bits
mod field_impl {

    use crate::{
        crh::FieldBasedHash,
        signature::FieldBasedSignatureScheme,
        Error,
    };
    use algebra::{
        Field, PrimeField, Group, UniformRand, to_bytes, ToBytes, FromBytes,
        AffineCurve, ProjectiveCurve, project,
    };
    use std::marker::PhantomData;
    use rand::Rng;
    use std::{
        ops::Neg,
        io::{Write, Read, Result as IoResult}
    };

    pub struct FieldBasedSchnorrSignatureScheme<
        F: PrimeField,
        G: ProjectiveCurve<BaseField = F>,
        H: FieldBasedHash<Data = F>
    >
    {
        _field:    PhantomData<F>,
        _group:    PhantomData<G>,
        _hash:     PhantomData<H>,
    }

    #[derive(Derivative)]
    #[derivative(
    Clone(bound = "F: PrimeField"),
    Default(bound = "F: PrimeField"),
    Eq(bound = "F: PrimeField"),
    PartialEq(bound = "F: PrimeField"),
    Debug(bound = "F: PrimeField")
    )]
    pub struct FieldBasedSchnorrSignature<F: PrimeField> {
        pub r:    F,
        pub s:    Vec<u8>,
    }

    impl<F: PrimeField> ToBytes for FieldBasedSchnorrSignature<F> {
        #[inline]
        fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
            self.r.write(&mut writer)?;
            self.s.write(&mut writer)
        }
    }

    impl<F: PrimeField> FromBytes for FieldBasedSchnorrSignature<F> {
        fn read<R: Read>(mut reader: R) -> IoResult<Self> {
            let r = F::read(&mut reader)?;
            let s = Vec::<u8>::read(&mut reader)?;
            Ok(Self{r, s})
        }
    }

    impl<F: PrimeField, G: ProjectiveCurve<BaseField = F>, H: FieldBasedHash<Data = F>> FieldBasedSignatureScheme for
    FieldBasedSchnorrSignatureScheme<F, G, H>
        where
            G::Affine: AffineCurve<BaseField = F>,
    {
        type Data = H::Data;
        type PublicKey = G;
        type SecretKey = G::ScalarField;
        type Signature = FieldBasedSchnorrSignature<F>;

        fn keygen<R: Rng>(rng: &mut R) -> Result<(Self::PublicKey, Self::SecretKey), Error>
        {
            let secret_key = G::ScalarField::rand(rng);
            let public_key = G::prime_subgroup_generator()
                .mul(&secret_key);
            Ok((public_key, secret_key))
        }

        fn sign(
            pk: &Self::PublicKey,
            sk: &Self::SecretKey,
            message: &[Self::Data],
        )-> Result<Self::Signature, Error>
        {
            // Compute k' = H(m || pk || sk)
            let pk = pk.into_affine();
            let sk_b = project::<G::ScalarField, F>(*sk)?;
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(message);
            hash_input.push(pk.get_x());
            hash_input.push(pk.get_y());
            hash_input.push(sk_b);
            let hr = H::evaluate(hash_input.as_ref())?;

            let k_prime = project::<F, G::ScalarField>(hr)?;

            // Assert k' != 0
            assert!(!k_prime.is_zero());

            //R = k' * G
            let r = G::prime_subgroup_generator()
                .mul(&k_prime).into_affine();

            //Set k = -k' if r.y is odd otherwise k = k'
            let k = if r.get_y().is_odd() { k_prime.neg() } else { k_prime };

            // Compute e = H(m || R.x || pk)
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(message);
            hash_input.push(r.get_x());
            hash_input.push(pk.get_x());
            hash_input.push(pk.get_y());
            let hr = H::evaluate(hash_input.as_ref())?;
            let e =  project::<F, G::ScalarField>(hr)?;

            let signature = FieldBasedSchnorrSignature {
                r: r.get_x(),
                s: to_bytes!((k + &(e * sk))).unwrap(),
            };

            Ok(signature)
        }

        fn verify(
            pk: &Self::PublicKey,
            message: &[Self::Data],
            signature: &Self::Signature
        )
            -> Result<bool, Error>
        {
            let s = G::ScalarField::read(signature.s.as_slice())?;
            let pk = pk.into_affine();

            debug_assert!(pk.is_in_correct_subgroup_assuming_on_curve());
            debug_assert!(s.pow(&G::ScalarField::characteristic()) == s);
            debug_assert!(signature.r.pow(&G::BaseField::characteristic()) == signature.r);

            // Compute e' = H(m || signature.r.x || pk)
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(message);
            hash_input.push(signature.r);
            hash_input.push(pk.get_x());
            hash_input.push(pk.get_y());
            let hr = H::evaluate(hash_input.as_ref())?;

            let e_prime =  project::<F, G::ScalarField>(hr)?;

            //Compute R' = s*G - e' * pk
            let r_prime = {
                let s_times_g = G::prime_subgroup_generator().mul(&s);
                let neg_e_times_pk = pk.neg().mul(e_prime);
                (s_times_g + &neg_e_times_pk).into_affine()
            };

            Ok((r_prime.get_x() == signature.r) && !r_prime.is_zero() && !r_prime.get_y().is_odd())
        }
    }
}

#[cfg(test)]
mod test {
    use algebra::curves::{
        mnt4753::G1Projective as MNT4G1Projective,
        mnt6753::G1Projective as MNT6G1Projective,
    };
    use algebra::fields::{
        mnt4753::Fr as MNT4Fr,
        mnt6753::Fr as MNT6Fr,
    };
    use crate::crh::{MNT4PoseidonHash, MNT6PoseidonHash};
    use crate::signature::FieldBasedSignatureScheme;
    use crate::signature::schnorr::field_impl::FieldBasedSchnorrSignatureScheme;
    use rand::{Rng, thread_rng};

    type SchnorrMNT4 = FieldBasedSchnorrSignatureScheme<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash>;
    type SchnorrMNT6 = FieldBasedSchnorrSignatureScheme<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash>;

    fn sign_and_verify<S: FieldBasedSignatureScheme>(message: &[S::Data]) {
        let rng = &mut thread_rng();
        let (pk, sk) = S::keygen(rng).unwrap();
        let sig = S::sign(&pk, &sk, &message).unwrap();
        assert!(S::verify(&pk, &message, &sig).unwrap());
    }

    fn failed_verification<S: FieldBasedSignatureScheme>(message: &[S::Data], bad_message: &[S::Data]) {
        let rng = &mut thread_rng();
        let (pk, sk) = S::keygen(rng).unwrap();
        let sig = S::sign(&pk, &sk, message).unwrap();
        assert!(!S::verify(&pk, bad_message, &sig).unwrap());
    }

    #[test]
    fn mnt4_schnorr_test() {
        let samples = 100;
        for _ in 0..samples {
            let rng = &mut thread_rng();
            let f: MNT4Fr = rng.gen();
            let g: MNT4Fr = rng.gen();
            sign_and_verify::<SchnorrMNT4>(&[f, g]);
            failed_verification::<SchnorrMNT4>(&[f], &[g]);
        }
    }

    #[test]
    fn mnt6_schnorr_test() {
        let samples = 100;
        for _ in 0..samples{
            let rng = &mut thread_rng();
            let f: MNT6Fr = rng.gen();
            let g: MNT6Fr = rng.gen();
            sign_and_verify::<SchnorrMNT6>(&[f, g]);
            failed_verification::<SchnorrMNT6>(&[f], &[g]);
        }
    }
}