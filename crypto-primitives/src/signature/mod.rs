use algebra::bytes::ToBytes;
use crate::Error;
use rand::Rng;
use std::hash::Hash;

#[cfg(feature = "r1cs")]
pub mod constraints;
#[cfg(feature = "r1cs")]
pub use constraints::*;


pub mod schnorr;

pub trait SignatureScheme {
    type Parameters: Clone + Send + Sync;
    type PublicKey: ToBytes + Hash + Eq + Clone + Default + Send + Sync;
    type SecretKey: ToBytes + Clone + Default;
    type Signature: Clone + Default + Send + Sync;

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

#[cfg(test)]
mod test {
    use crate::{signature::schnorr::SchnorrSignature, SignatureScheme};
    use algebra::{
        curves::edwards_sw6::EdwardsAffine as Edwards, groups::Group, to_bytes, ToBytes,
    };
    use blake2::Blake2s;
    use rand::thread_rng;
    use algebra::UniformRand;

    fn sign_and_verify<S: SignatureScheme>(message: &[u8]) {
        let rng = &mut thread_rng();
        let parameters = S::setup::<_>(rng).unwrap();
        let (pk, sk) = S::keygen(&parameters, rng).unwrap();
        let sig = S::sign(&parameters, &sk, &message, rng).unwrap();
        assert!(S::verify(&parameters, &pk, &message, &sig).unwrap());
    }

    fn failed_verification<S: SignatureScheme>(message: &[u8], bad_message: &[u8]) {
        let rng = &mut thread_rng();
        let parameters = S::setup::<_>(rng).unwrap();
        let (pk, sk) = S::keygen(&parameters, rng).unwrap();
        let sig = S::sign(&parameters, &sk, message, rng).unwrap();
        assert!(!S::verify(&parameters, &pk, bad_message, &sig).unwrap());
    }

    fn randomize_and_verify<S: SignatureScheme>(message: &[u8], randomness: &[u8]) {
        let rng = &mut thread_rng();
        let parameters = S::setup::<_>(rng).unwrap();
        let (pk, sk) = S::keygen(&parameters, rng).unwrap();
        let sig = S::sign(&parameters, &sk, message, rng).unwrap();
        assert!(S::verify(&parameters, &pk, message, &sig).unwrap());
        let randomized_pk = S::randomize_public_key(&parameters, &pk, randomness).unwrap();
        let randomized_sig = S::randomize_signature(&parameters, &sig, randomness).unwrap();
        assert!(S::verify(&parameters, &randomized_pk, &message, &randomized_sig).unwrap());
    }

    #[test]
    fn schnorr_signature_test() {
        let message = "Hi, I am a Schnorr signature!";
        let rng = &mut thread_rng();
        sign_and_verify::<SchnorrSignature<Edwards, Blake2s>>(message.as_bytes());
        failed_verification::<SchnorrSignature<Edwards, Blake2s>>(
            message.as_bytes(),
            "Bad message".as_bytes(),
        );
        let random_scalar = to_bytes!(<Edwards as Group>::ScalarField::rand(rng)).unwrap();
        randomize_and_verify::<SchnorrSignature<Edwards, Blake2s>>(
            message.as_bytes(),
            &random_scalar.as_slice(),
        );
    }
}
