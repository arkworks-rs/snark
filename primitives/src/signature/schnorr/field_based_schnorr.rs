use crate::{
    compute_truncation_size, crh::FieldBasedHash, signature::FieldBasedSignatureScheme, Error,
};
use algebra::{
    convert, leading_zeros, serialize::*, Field, FromBytes, FromBytesChecked, Group, PrimeField,
    ProjectiveCurve, SemanticallyValid, ToBits, ToBytes, ToConstraintField, UniformRand,
};
use rand::distributions::{Distribution, Standard};
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::io::{Error as IoError, ErrorKind, Read, Result as IoResult, Write};
use std::marker::PhantomData;

#[allow(dead_code)]
pub struct FieldBasedSchnorrSignatureScheme<F: PrimeField, G: Group, H: FieldBasedHash> {
    _field: PhantomData<F>,
    _group: PhantomData<G>,
    _hash: PhantomData<H>,
}

#[derive(Derivative)]
#[derivative(
    Copy(bound = "F: PrimeField, G: Group"),
    Clone(bound = "F: PrimeField, G: Group"),
    Default(bound = "F: PrimeField, G: Group"),
    Eq(bound = "F: PrimeField, G: Group"),
    PartialEq(bound = "F: PrimeField, G: Group"),
    Debug(bound = "F: PrimeField, G: Group")
)]
#[derive(Serialize, Deserialize)]
#[serde(bound(serialize = "F: PrimeField, G: Group"))]
#[serde(bound(deserialize = "F: PrimeField, G: Group"))]
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct FieldBasedSchnorrSignature<F: PrimeField, G: Group> {
    pub e: F,
    pub s: F,
    #[serde(skip)]
    _group: PhantomData<G>,
}

impl<F: PrimeField, G: Group> FieldBasedSchnorrSignature<F, G> {
    #[allow(dead_code)]
    pub fn new(e: F, s: F) -> Self {
        Self {
            e,
            s,
            _group: PhantomData,
        }
    }
}

impl<F: PrimeField, G: Group> ToBytes for FieldBasedSchnorrSignature<F, G> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.e.write(&mut writer)?;
        self.s.write(&mut writer)
    }
}

impl<F: PrimeField, G: Group> FromBytes for FieldBasedSchnorrSignature<F, G> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let e = F::read(&mut reader)?;
        let s = F::read(&mut reader)?;
        Ok(Self {
            e,
            s,
            _group: PhantomData,
        })
    }
}

impl<F: PrimeField, G: Group> FromBytesChecked for FieldBasedSchnorrSignature<F, G> {
    fn read_checked<R: Read>(mut reader: R) -> IoResult<Self> {
        let e = F::read_checked(&mut reader)
            .map_err(|err| IoError::new(ErrorKind::InvalidData, format!("invalid sig.e: {}", err)))
            .and_then(|e| {
                let e_bits = e.write_bits();
                let e_leading_zeros = leading_zeros(e_bits.as_slice()) as usize;
                if (F::size_in_bits() - e_leading_zeros) >= G::ScalarField::size_in_bits() {
                    return Err(IoError::new(
                        ErrorKind::InvalidData,
                        format!(
                            "Invalid bit-length for signature.e: {}",
                            e_bits.len() - e_leading_zeros
                        ),
                    ));
                }
                Ok(e)
            })?;
        let s = F::read_checked(&mut reader)
            .map_err(|err| IoError::new(ErrorKind::InvalidData, format!("invalid sig.e: {}", err)))
            .and_then(|s| {
                let s_bits = s.write_bits();
                let s_leading_zeros = leading_zeros(s_bits.as_slice()) as usize;
                if (G::ScalarField::size_in_bits() - s_leading_zeros) >= F::size_in_bits() {
                    return Err(IoError::new(
                        ErrorKind::InvalidData,
                        format!(
                            "Invalid bit-length for signature.s: {}",
                            s_bits.len() - s_leading_zeros
                        ),
                    ));
                }
                Ok(s)
            })?;
        Ok(Self {
            e,
            s,
            _group: PhantomData,
        })
    }
}

impl<F: PrimeField, G: Group> SemanticallyValid for FieldBasedSchnorrSignature<F, G> {
    fn is_valid(&self) -> bool {
        self.e.is_valid()
            && {
                //Checks e had proper bit-length when converted into a G::ScalarField element
                let e_bits = self.e.write_bits();
                let e_leading_zeros = leading_zeros(e_bits.as_slice()) as usize;
                F::size_in_bits() - e_leading_zeros < G::ScalarField::size_in_bits()
            }
            && self.s.is_valid()
            && {
                //Checks s had proper bit-length when converted into a F element
                let s_bits = self.s.write_bits();
                let s_leading_zeros = leading_zeros(s_bits.as_slice()) as usize;
                G::ScalarField::size_in_bits() - s_leading_zeros < F::size_in_bits()
            }
    }
}

#[derive(Derivative)]
#[derivative(
    Copy(bound = "G: Group"),
    Clone(bound = "G: Group"),
    Default(bound = "G: Group"),
    Hash(bound = "G: Group"),
    Eq(bound = "G: Group"),
    PartialEq(bound = "G: Group"),
    Debug(bound = "G: Group")
)]
#[derive(Serialize, Deserialize)]
#[serde(bound(serialize = "G: Group"))]
#[serde(bound(deserialize = "G: Group"))]
#[serde(transparent)]
pub struct FieldBasedSchnorrPk<G: Group>(pub G);

impl<G: Group> Distribution<FieldBasedSchnorrPk<G>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> FieldBasedSchnorrPk<G> {
        let pk = G::rand(rng);
        FieldBasedSchnorrPk::<G>(pk)
    }
}

impl<G: Group> ToBytes for FieldBasedSchnorrPk<G> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}

impl<G: Group> FromBytes for FieldBasedSchnorrPk<G> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let pk = G::read(&mut reader)?;
        Ok(Self(pk))
    }
}

impl<G: Group> FromBytesChecked for FieldBasedSchnorrPk<G> {
    fn read_checked<R: Read>(mut reader: R) -> IoResult<Self> {
        let pk = G::read_checked(&mut reader)
            .map_err(|e| IoError::new(ErrorKind::InvalidData, format!("invalid schnorr pk: {}", e)))
            .and_then(|p| {
                if p.is_zero() {
                    return Err(IoError::new(
                        ErrorKind::InvalidData,
                        "invalid schnorr pk: point at infinity",
                    ));
                }
                Ok(p)
            })?;
        Ok(Self(pk))
    }
}

impl<G: Group> SemanticallyValid for FieldBasedSchnorrPk<G> {
    #[inline]
    fn is_valid(&self) -> bool {
        self.0.is_valid() &&
        // GingerLib only accepts non-trivial Schnorr public keys. This is usually
        // good practice to avoid using obvious weak keys, and helps preventing
        // exceptional cases if using incomplete arithmetics.
        !self.0.is_zero()
    }
}

// Low-level crypto for the length-restricted Schnorr Signature, does not perform any
// input validity check. It's responsibility of the caller to do so, through keyverify()
// function for the PublicKey, read() or is_valid() functions for FieldBasedSchnorrSignature.
impl<F: PrimeField, G: ProjectiveCurve + ToConstraintField<F>, H: FieldBasedHash<Data = F>>
    FieldBasedSignatureScheme for FieldBasedSchnorrSignatureScheme<F, G, H>
{
    type Data = H::Data;
    type PublicKey = FieldBasedSchnorrPk<G>;
    type SecretKey = G::ScalarField;
    type Signature = FieldBasedSchnorrSignature<F, G>;

    fn keygen<R: Rng>(rng: &mut R) -> (Self::PublicKey, Self::SecretKey) {
        let secret_key = loop {
            let r = G::ScalarField::rand(rng);
            // Reject sk = 0 to avoid generating obviously weak keypair. See keyverify() function
            // for additional explanations.
            if !r.is_zero() {
                break (r);
            }
        };
        let public_key = G::prime_subgroup_generator().mul(&secret_key);
        (FieldBasedSchnorrPk(public_key), secret_key)
    }

    fn get_public_key(sk: &Self::SecretKey) -> Self::PublicKey {
        FieldBasedSchnorrPk(G::prime_subgroup_generator().mul(sk))
    }

    fn sign<R: Rng>(
        rng: &mut R,
        pk: &Self::PublicKey,
        sk: &Self::SecretKey,
        message: Self::Data,
    ) -> Result<Self::Signature, Error> {
        let required_leading_zeros_e = compute_truncation_size(
            F::size_in_bits() as i32,
            G::ScalarField::size_in_bits() as i32,
        );

        let required_leading_zeros_s = compute_truncation_size(
            G::ScalarField::size_in_bits() as i32,
            F::size_in_bits() as i32,
        );

        //Affine coordinates of pk
        let pk_coords = pk.0.to_field_elements()?;

        let (e, s) = loop {
            //Sample random element
            let k = G::ScalarField::rand(rng);

            //R = k * G
            let r = G::prime_subgroup_generator().mul(&k);

            //Affine coordinates of R (even if R is infinity)
            let r_coords = r.to_field_elements()?;

            // Compute e = H(m || R || pk.x)
            let e = {
                let mut digest = H::init_constant_length(4, None);
                digest.update(message);
                r_coords.into_iter().for_each(|coord| {
                    digest.update(coord);
                });
                digest.update(pk_coords[0]);
                digest.finalize()
            }?;

            let e_bits = e.write_bits();
            let e_leading_zeros = leading_zeros(e_bits.as_slice()) as usize;

            //Enforce e bit length is strictly smaller than G::ScalarField modulus bit length
            if e_leading_zeros < required_leading_zeros_e {
                continue;
            };

            //We can now safely convert it to the other field
            let e_conv = convert::<G::ScalarField>(e_bits)?;

            //Enforce s bit length is strictly smaller than F modulus bit length
            let s = k + &(e_conv * sk);
            let s_bits = s.write_bits();
            let s_leading_zeros = leading_zeros(s_bits.as_slice()) as usize;

            if s_leading_zeros < required_leading_zeros_s {
                continue;
            };

            let s_conv = convert::<F>(s_bits)?;

            break (e, s_conv);
        };

        Ok(FieldBasedSchnorrSignature {
            e,
            s,
            _group: PhantomData,
        })
    }

    fn verify(
        pk: &Self::PublicKey,
        message: Self::Data,
        signature: &Self::Signature,
    ) -> Result<bool, Error> {
        let pk_coords = pk.0.to_field_elements()?;

        //Compute R' = s*G - e * pk
        let r_prime = {
            let s_bits = signature.s.write_bits();
            let s_conv = convert::<G::ScalarField>(s_bits)?;

            let e_bits = signature.e.write_bits();
            let e_conv = convert::<G::ScalarField>(e_bits)?;

            let s_times_g = G::prime_subgroup_generator().mul(&s_conv);
            let neg_e_times_pk = pk.0.neg().mul(&e_conv);
            s_times_g + &neg_e_times_pk
        };

        let r_prime_coords = r_prime.to_field_elements()?;

        // Compute e' = H(m || R' || pk.x)
        let e_prime = {
            let mut digest = H::init_constant_length(4, None);
            digest.update(message);
            r_prime_coords.into_iter().for_each(|coord| {
                digest.update(coord);
            });
            digest.update(pk_coords[0]);
            digest.finalize()
        }?;

        Ok(signature.e == e_prime)
    }

    #[inline]
    fn keyverify(pk: &Self::PublicKey) -> bool {
        pk.is_valid()
    }
}

#[cfg(test)]
mod test {
    use crate::crh::{MNT4PoseidonHash, MNT6PoseidonHash};
    use crate::signature::schnorr::field_based_schnorr::FieldBasedSchnorrSignatureScheme;
    use crate::signature::FieldBasedSignatureScheme;
    use algebra::curves::{
        mnt4753::G1Projective as MNT4G1Projective, mnt6753::G1Projective as MNT6G1Projective,
    };
    use algebra::fields::{mnt4753::Fr as MNT4Fr, mnt6753::Fr as MNT6Fr};
    use algebra::{to_bytes, FromBytes, FromBytesChecked, SemanticallyValid, ToBytes};
    use rand::{thread_rng, Rng};

    type SchnorrMNT4 = FieldBasedSchnorrSignatureScheme<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash>;
    type SchnorrMNT6 = FieldBasedSchnorrSignatureScheme<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash>;

    fn sign_and_verify<S: FieldBasedSignatureScheme, R: Rng>(rng: &mut R, message: S::Data) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));
        assert_eq!(pk, S::get_public_key(&sk));
        let sig = S::sign(rng, &pk, &sk, message).unwrap();
        assert!(sig.is_valid());
        assert!(S::verify(&pk, message, &sig).unwrap());

        //Serialization/deserialization test
        let sig_serialized = to_bytes!(sig).unwrap();
        let sig_deserialized =
            <S as FieldBasedSignatureScheme>::Signature::read(sig_serialized.as_slice()).unwrap();
        assert_eq!(sig, sig_deserialized);
        assert!(<S as FieldBasedSignatureScheme>::Signature::read_checked(
            sig_serialized.as_slice()
        )
        .is_ok());
        assert!(S::verify(&pk, message, &sig_deserialized).unwrap());
    }

    fn failed_verification<S: FieldBasedSignatureScheme, R: Rng>(
        rng: &mut R,
        message: S::Data,
        bad_message: S::Data,
    ) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));
        assert_eq!(pk, S::get_public_key(&sk));

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
            sign_and_verify::<SchnorrMNT4, _>(rng, f);
            failed_verification::<SchnorrMNT4, _>(rng, f, g);
        }
    }

    #[test]
    fn mnt6_schnorr_test() {
        let rng = &mut thread_rng();
        let samples = 100;
        for _ in 0..samples {
            let f: MNT6Fr = rng.gen();
            let g: MNT6Fr = rng.gen();
            sign_and_verify::<SchnorrMNT6, _>(rng, f);
            failed_verification::<SchnorrMNT6, _>(rng, f, g);
        }
    }
}
