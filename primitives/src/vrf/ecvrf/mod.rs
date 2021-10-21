use algebra::{Field, PrimeField, FpParameters, convert, leading_zeros, Group, AffineCurve, ProjectiveCurve,
              ToBytes, to_bytes, ToBits, UniformRand, ToConstraintField, FromBytes};
use crate::{crh::{
    FieldBasedHash, FixedLengthCRH,
}, vrf::FieldBasedVrf, Error, CryptoError, compute_truncation_size};
use std::marker::PhantomData;
use rand::Rng;
use std::io::{Write, Read, Result as IoResult};


pub struct FieldBasedEcVrf<
    F: PrimeField,
    G: Group,
    FH: FieldBasedHash,
    GH: FixedLengthCRH,
>
{
    _field:         PhantomData<F>,
    _group:         PhantomData<G>,
    _field_hash:    PhantomData<FH>,
    _group_hash:    PhantomData<GH>,
}

#[derive(Derivative)]
#[derivative(
Copy(bound = "F: PrimeField, G: ProjectiveCurve"),
Clone(bound = "F: PrimeField, G: ProjectiveCurve"),
Default(bound = "F: PrimeField, G: ProjectiveCurve"),
Eq(bound = "F: PrimeField, G: ProjectiveCurve"),
PartialEq(bound = "F: PrimeField, G: ProjectiveCurve"),
Debug(bound = "F: PrimeField, G: ProjectiveCurve")
)]
pub struct FieldBasedEcVrfProof<F: PrimeField, G: ProjectiveCurve> {
    pub gamma:  G,
    pub c:      F,
    pub s:      F,
}

impl<F: PrimeField, G: ProjectiveCurve> ToBytes for FieldBasedEcVrfProof<F, G> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.gamma.into_affine().write(&mut writer)?;
        self.c.write(&mut writer)?;
        self.s.write(&mut writer)
    }
}

impl<F: PrimeField, G: ProjectiveCurve> FromBytes for FieldBasedEcVrfProof<F, G> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let gamma = G::Affine::read(&mut reader)?;
        let c = F::read(&mut reader)?;
        let s = F::read(&mut reader)?;
        Ok(Self{ gamma: gamma.into_projective(), c, s })
    }
}

impl<F, G, FH, GH> FieldBasedVrf for FieldBasedEcVrf<F, G, FH, GH>
    where
        F: PrimeField,
        G: ProjectiveCurve + ToConstraintField<F>,
        FH: FieldBasedHash<Data = F>,
        GH: FixedLengthCRH<Output = G>,
{
    type Data = FH::Data;
    type PublicKey = G;
    type SecretKey = G::ScalarField;
    type Proof = FieldBasedEcVrfProof<F, G>;
    type GHParams = GH::Parameters;

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

    fn prove<R: Rng>(
        rng:               &mut R,
        group_hash_params: &Self::GHParams,
        pk:                &Self::PublicKey,
        sk:                &Self::SecretKey,
        message:           &[Self::Data],
    )-> Result<Self::Proof, Error>
    {
        //Compute mh = hash_to_curve(message)
        let mut message_bytes = Vec::new();
        for (i, field_element) in message.iter().enumerate() {
            // The reason for a secure de-packing is not collision resistance (the non-restricted variant
            // would be still), but that inside the circuit a field element might be proven to hash to
            // one of two possible fingerprints (as there might be two different byte sequences satisfying
            // the depacking constraint mod q). Hence via SNARKs the output of the VRF is not unique and
            // can be chosen between two possible outputs, which is what we definitely do not want in the
            // application of the VRF (the VRF is now rather a verifiable random relation, not function).
            if field_element.into_repr_raw() >= F::Params::MODULUS {
                return Err(Box::new(CryptoError::InvalidElement(format!("message_element_{}", i).to_owned())));
            }
            message_bytes.extend_from_slice(to_bytes!(field_element).unwrap().as_slice())
        }

        let message_on_curve = GH::evaluate(group_hash_params, message_bytes.as_slice())?;

        //Compute gamma = message_on_curve^sk
        let gamma = message_on_curve.mul(sk);

        let (c, s) = loop {

            //Choose random scalar
            let r = G::ScalarField::rand(rng);

            if r.is_zero() {continue};

            //Compute a = g^r
            let a = G::prime_subgroup_generator().mul(&r);

            //Compute b = message_on_curve^r
            let b = message_on_curve.mul(&r);

            //Compute c = H(m||pk.x||a.x||b.x)
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(message);
            hash_input.push(pk.to_field_elements().unwrap()[0]);
            hash_input.push(a.to_field_elements().unwrap()[0]);
            hash_input.push(b.to_field_elements().unwrap()[0]);
            let c = FH::evaluate(hash_input.as_ref())?;
            let c_bits = c.write_bits();
            let c_leading_zeros = leading_zeros(c_bits.clone()) as usize;
            let required_leading_zeros = compute_truncation_size(
                F::size_in_bits() as i32,
                G::ScalarField::size_in_bits() as i32,
            );

            //Enforce c bit length is strictly smaller than G::ScalarField modulus bit length
            if c_leading_zeros < required_leading_zeros {continue};

            let c_conv = convert::<G::ScalarField>(c_bits)?;

            //Compute s = r + sk * c
            let s = r + &((*sk) * &c_conv);
            let s_bits = s.write_bits();
            let s_leading_zeros = leading_zeros(s_bits.clone()) as usize;
            let required_leading_zeros = compute_truncation_size(
                G::ScalarField::size_in_bits() as i32,
                F::size_in_bits() as i32,
            );

            if s_leading_zeros < required_leading_zeros {continue};

            let s_conv = convert::<F>(s_bits)?;

            break (c, s_conv)
        };

        Ok(FieldBasedEcVrfProof {gamma, c, s})
    }

    fn proof_to_hash(
        group_hash_params: &Self::GHParams,
        pk:                &Self::PublicKey,
        message:           &[Self::Data],
        proof:             &Self::Proof
    )
        -> Result<Self::Data, Error>
    {

        //Checks
        let c_bits = proof.c.write_bits();
        let c_leading_zeros = leading_zeros(c_bits.clone()) as usize;
        if (F::size_in_bits() - c_leading_zeros) >= G::ScalarField::size_in_bits(){
            return Err(Box::new(CryptoError::IncorrectInputLength("proof.c".to_owned(), c_bits.len() - c_leading_zeros)))
        }

        let s_bits = proof.s.write_bits();
        let s_leading_zeros = leading_zeros(s_bits.clone()) as usize;
        if (G::ScalarField::size_in_bits() - s_leading_zeros) >= F::size_in_bits(){
            return Err(Box::new(CryptoError::IncorrectInputLength("proof.s".to_owned(), s_bits.len() - s_leading_zeros)))
        }

        if !proof.gamma.group_membership_test() {
            return Err(Box::new(CryptoError::NotPrimeOrder("proof.gamma".to_owned())))
        }

        //Compute mh = hash_to_curve(message)
        let mut message_bytes = Vec::new();
        for (i, field_element) in message.iter().enumerate() {
            // The reason for a secure de-packing is not collision resistance (the non-restricted variant
            // would be still), but that inside the circuit a field element might be proven to hash to
            // one of two possible fingerprints (as there might be two different byte sequences satisfying
            // the depacking constraint mod q). Hence via SNARKs the output of the VRF is not unique and
            // can be chosen between two possible outputs, which is what we definitely do not want in the
            // application of the VRF (the VRF is now rather a verifiable random relation, not function).
            if field_element.into_repr_raw() >= F::Params::MODULUS {
                return Err(Box::new(CryptoError::InvalidElement(format!("message_element_{}", i).to_owned())));
            }
            message_bytes.extend_from_slice(to_bytes!(field_element).unwrap().as_slice())
        }

        let message_on_curve = GH::evaluate(group_hash_params, message_bytes.as_slice())?;

        let c_conv = convert::<G::ScalarField>(c_bits)?;
        let s_conv = convert::<G::ScalarField>(s_bits)?;

        //Compute u = g^s - pk^c
        let u = G::prime_subgroup_generator().mul(&s_conv) - &(pk.mul(&c_conv));

        //Compute v = mh^s - gamma^c
        let v = message_on_curve.mul(&s_conv) - &proof.gamma.mul(&c_conv);

        //Compute c' = H(m||pk.x||u.x||v.x)
        let mut hash_input = Vec::new();
        let pk_coords = pk.to_field_elements()?;
        hash_input.extend_from_slice(message);
        hash_input.push(pk_coords[0]);
        hash_input.push(u.to_field_elements().unwrap()[0]);
        hash_input.push(v.to_field_elements().unwrap()[0]);
        let c_prime = FH::evaluate(hash_input.as_ref())?;

        //Verify valid proof
        match proof.c == c_prime {
            false => Err(Box::new(CryptoError::FailedVerification)),
            true => {
                let gamma_coords = proof.gamma.to_field_elements().unwrap();

                //Compute VRF output
                hash_input = Vec::new();
                hash_input.extend_from_slice(message);
                hash_input.extend_from_slice(gamma_coords.as_slice());
                let output = FH::evaluate(hash_input.as_ref())?;

                //Return VRF output
                Ok(output)
            }
        }
    }

    fn keyverify(
        pk: &Self::PublicKey,
    ) -> bool {
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
    use algebra::{ToBytes, FromBytes, to_bytes};
    use crate::{
        crh::{
            MNT4PoseidonHash, MNT6PoseidonHash,
            bowe_hopwood::BoweHopwoodPedersenCRH,
            pedersen::PedersenWindow,
        },
        vrf::{
            FieldBasedVrf,
            ecvrf::FieldBasedEcVrf,
        },
        FixedLengthCRH
    };
    use rand::{Rng, thread_rng};

    #[derive(Clone)]
    struct TestWindow {}
    impl PedersenWindow for TestWindow {
        const WINDOW_SIZE: usize = 128;
        const NUM_WINDOWS: usize = 2;
    }

    type BHMNT4 = BoweHopwoodPedersenCRH<MNT4G1Projective, TestWindow>;
    type BHMNT6 = BoweHopwoodPedersenCRH<MNT6G1Projective, TestWindow>;

    type EcVrfMNT4 = FieldBasedEcVrf<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash, BHMNT6>;
    type EcVrfMNT6 = FieldBasedEcVrf<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash, BHMNT4>;

    fn prove_and_verify<S: FieldBasedVrf, R: Rng>(rng: &mut R, message: &[S::Data], pp: &S::GHParams) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));
        assert_eq!(pk, S::get_public_key(&sk));

        let proof = S::prove(rng, pp, &pk, &sk, &message).unwrap();
        assert!(S::proof_to_hash(pp, &pk, &message, &proof).is_ok());

        //Serialization/deserialization test
        let proof_serialized = to_bytes!(proof).unwrap();
        let proof_deserialized = <S as FieldBasedVrf>::Proof::read(proof_serialized.as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);
        assert!(S::proof_to_hash(pp, &pk, &message, &proof_deserialized).is_ok());
    }

    fn failed_verification<S: FieldBasedVrf, R: Rng>(rng: &mut R, message: &[S::Data], bad_message: &[S::Data], pp: &S::GHParams) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));
        assert_eq!(pk, S::get_public_key(&sk));

        //Attempt to verify proof for a different message
        let proof = S::prove(rng, pp, &pk, &sk, message).unwrap();
        assert!(S::proof_to_hash(pp, &pk, bad_message, &proof).is_err());

        //Attempt to verify different proof for a message
        let bad_proof = S::prove(rng, pp, &pk, &sk, bad_message).unwrap();
        assert!(S::proof_to_hash(pp, &pk, message, &bad_proof).is_err());

        //Attempt to verify proof for a message with different pk
        let (new_pk, _) = S::keygen(rng);
        assert!(S::proof_to_hash(pp, &new_pk, message, &proof).is_err());
    }

    #[test]
    fn mnt4_ecvrf_test() {
        let rng = &mut thread_rng();
        let pp = <BHMNT6 as FixedLengthCRH>::setup(rng).unwrap();
        let samples = 100;
        for _ in 0..samples {
            let f: MNT4Fr = rng.gen();
            let g: MNT4Fr = rng.gen();
            prove_and_verify::<EcVrfMNT4, _>(rng, &[f], &pp);
            failed_verification::<EcVrfMNT4, _>(rng, &[f], &[g], &pp);
        }
    }

    #[test]
    fn mnt6_ecvrf_test() {
        let rng = &mut thread_rng();
        let pp = <BHMNT4 as FixedLengthCRH>::setup(rng).unwrap();
        let samples = 100;
        for _ in 0..samples {
            let f: MNT6Fr = rng.gen();
            let g: MNT6Fr = rng.gen();
            prove_and_verify::<EcVrfMNT6, _>(rng, &[f], &pp);
            failed_verification::<EcVrfMNT6, _>(rng, &[f], &[g], &pp);
        }
    }
}