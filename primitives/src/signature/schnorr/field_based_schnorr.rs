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

    fn get_public_key(sk: &Self::SecretKey) -> Self::PublicKey {
        G::prime_subgroup_generator().mul(sk)
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
            let e_bits = e.write_bits();
            let e_leading_zeros = leading_zeros(e_bits.clone()) as usize;
            let required_leading_zeros = compute_truncation_size(
                F::size_in_bits() as i32,
                G::ScalarField::size_in_bits() as i32,
            );

            //Enforce e bit length is strictly smaller than G::ScalarField modulus bit length
            if e_leading_zeros < required_leading_zeros {continue};

            //We can now safely convert it to the other field
            let e_conv = convert::<G::ScalarField>(e_bits)?;

            //Enforce s bit length is strictly smaller than F modulus bit length
            let s = k + &(e_conv * sk);
            let s_bits = s.write_bits();
            let s_leading_zeros = leading_zeros(s_bits.clone()) as usize;
            let required_leading_zeros = compute_truncation_size(
                G::ScalarField::size_in_bits() as i32,
                F::size_in_bits() as i32,
            );

            if s_leading_zeros < required_leading_zeros {continue};

            let s_conv = convert::<F>(s_bits)?;

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
            let s_conv = convert::<G::ScalarField>(s_bits)?;
            let e_conv =  convert::<G::ScalarField>(e_bits)?;
            let s_times_g = G::prime_subgroup_generator().mul(&s_conv);
            let neg_e_times_pk = pk.neg().mul(&e_conv);
            s_times_g + &neg_e_times_pk
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
    use algebra::{ToBytes, to_bytes, FromBytes};
    use crate::crh::{MNT4PoseidonHash, MNT6PoseidonHash};
    use crate::signature::FieldBasedSignatureScheme;
    use crate::signature::schnorr::field_based_schnorr::FieldBasedSchnorrSignatureScheme;
    use rand::{Rng, thread_rng};

    type SchnorrMNT4 = FieldBasedSchnorrSignatureScheme<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash>;
    type SchnorrMNT6 = FieldBasedSchnorrSignatureScheme<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash>;

    fn sign_and_verify<S: FieldBasedSignatureScheme, R: Rng>(rng: &mut R, message: &[S::Data]) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));
        assert_eq!(pk, S::get_public_key(&sk));
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
}