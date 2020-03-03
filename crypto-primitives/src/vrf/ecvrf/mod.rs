use algebra::{Field, PrimeField, convert, leading_zeros, Group, AffineCurve, ProjectiveCurve,
              ToBytes, to_bytes, ToBits, UniformRand, ToConstraintField, FromBytes};
use crate::{
    crh::{
    FieldBasedHash, FixedLengthCRH,
    },
    vrf::FieldBasedVrf,
    Error, CryptoError,
    compute_truncation_size,
};
use std::marker::PhantomData;
use rand::Rng;
use std::io::{Write, Read, Result as IoResult};

#[cfg(feature = "r1cs")]
pub mod constraints;

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
Clone(bound = "F: PrimeField, G: ProjectiveCurve"),
Default(bound = "F: PrimeField, G: ProjectiveCurve"),
Eq(bound = "F: PrimeField, G: ProjectiveCurve"),
PartialEq(bound = "F: PrimeField, G: ProjectiveCurve"),
Debug(bound = "F: PrimeField, G: ProjectiveCurve")
)]
pub struct FieldBasedEcVrfProof<F: PrimeField, G: Group> {
    pub gamma:  G,
    pub c:      F,
    pub s:      F,
}

impl<F: PrimeField, G: Group> ToBytes for FieldBasedEcVrfProof<F, G> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.gamma.write(&mut writer)?;
        self.c.write(&mut writer)?;
        self.s.write(&mut writer)
    }
}

impl<F: PrimeField, G: Group> FromBytes for FieldBasedEcVrfProof<F, G> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let gamma = G::read(&mut reader)?;
        let c = F::read(&mut reader)?;
        let s = F::read(&mut reader)?;
        Ok(Self{ gamma, c, s })
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

    fn keygen<R: Rng>(rng: &mut R) -> Result<(Self::PublicKey, Self::SecretKey), Error>
    {
        let secret_key = G::ScalarField::rand(rng);
        let public_key = G::prime_subgroup_generator()
            .mul(&secret_key);
        Ok((public_key, secret_key))
    }

    fn prove<R: Rng>(
        rng:     &mut R,
        pp:      &Self::GHParams,
        pk:      &Self::PublicKey,
        sk:      &Self::SecretKey,
        message: &[Self::Data],
    )-> Result<Self::Proof, Error>
    {
        //Compute mh = hash_to_curve(message)
        let mut message_bytes = Vec::new();
        message.iter().for_each(|f| message_bytes.extend_from_slice(to_bytes!(f).unwrap().as_slice()));
        let mh = GH::evaluate(pp, message_bytes.as_slice())?;

        //Compute gamma = mh^sk
        let gamma = mh.mul(sk);

        let (c, s) = loop {

            //Choose random scalar
            let r = G::ScalarField::rand(rng);

            //Compute a = g^r
            let a = G::prime_subgroup_generator().mul(&r);

            //Compute b = mh^r
            let b = mh.mul(&r);

            //Compute c = H(m||pk||a||b)
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(message);
            hash_input.extend_from_slice(pk.to_field_elements().unwrap().as_slice());
            hash_input.extend_from_slice(a.to_field_elements().unwrap().as_slice());
            hash_input.extend_from_slice(b.to_field_elements().unwrap().as_slice());
            let c = FH::evaluate(hash_input.as_ref())?;

            let c_leading_zeros = leading_zeros(c.write_bits()) as usize;
            let required_leading_zeros = compute_truncation_size(
                F::size_in_bits() as i32,
                G::ScalarField::size_in_bits() as i32,
            );

            //Enforce c bit length is strictly smaller than G::ScalarField modulus bit length
            if c_leading_zeros < required_leading_zeros {continue};

            let c_conv = convert::<F, G::ScalarField>(c)?;

            //Compute s = r + sk * c
            let s = r + &((*sk) * &c_conv);

            let s_leading_zeros = leading_zeros(s.write_bits()) as usize;
            let required_leading_zeros = compute_truncation_size(
                G::ScalarField::size_in_bits() as i32,
                F::size_in_bits() as i32,
            );

            if s_leading_zeros < required_leading_zeros {continue};

            let s_conv = convert::<G::ScalarField, F>(s)?;

            break (c, s_conv)
        };

        Ok(FieldBasedEcVrfProof {gamma, c, s})
    }

    fn verify(
        pp:      &Self::GHParams,
        pk:      &Self::PublicKey,
        message: &[Self::Data],
        proof:   &Self::Proof
    )
        -> Result<Self::Data, Error>
    {

        //Checks
        assert!(
            (F::size_in_bits() - leading_zeros(proof.c.write_bits()) as usize)
                < G::ScalarField::size_in_bits()
        );
        assert!(
            (G::ScalarField::size_in_bits() - leading_zeros(proof.s.write_bits()) as usize)
                < F::size_in_bits()
        );

        //Debug checks: should they be promoted to actual checks ?
        debug_assert!(pk.into_affine().is_in_correct_subgroup_assuming_on_curve());
        debug_assert!(proof.gamma.into_affine().is_in_correct_subgroup_assuming_on_curve());
        debug_assert!(proof.s.pow(&G::BaseField::characteristic()) == proof.s); //Does this check make sense ?
        debug_assert!(proof.c.pow(&G::BaseField::characteristic()) == proof.c);

        //Compute mh = hash_to_curve(message)
        let mut message_bytes = Vec::new();
        message.iter().for_each(|f| message_bytes.extend_from_slice(to_bytes!(f).unwrap().as_slice()));
        let mh = GH::evaluate(pp, message_bytes.as_slice())?;

        let c_conv = convert::<F, G::ScalarField>(proof.c)?;
        let s_conv = convert::<F, G::ScalarField>(proof.s)?;

        //Compute u = g^s - pk^c
        let u = G::prime_subgroup_generator().mul(&s_conv) - &(pk.mul(&c_conv));

        //Compute v = mh^s - gamma^c
        let v = mh.mul(&s_conv) - &proof.gamma.mul(&c_conv);

        //Compute c' = H(m||pk||u||v)
        let mut hash_input = Vec::new();
        let pk_coords = pk.to_field_elements()?;
        hash_input.extend_from_slice(message);
        hash_input.extend_from_slice(pk_coords.as_slice());
        hash_input.extend_from_slice(u.to_field_elements().unwrap().as_slice());
        hash_input.extend_from_slice(v.to_field_elements().unwrap().as_slice());
        let c_prime = FH::evaluate(hash_input.as_ref())?;

        //Verify valid proof
        match proof.c == c_prime {
            false => Err(Box::new(CryptoError::FailedVerification)),
            true => {

                //Compute VRF output
                hash_input = Vec::new();
                hash_input.extend_from_slice(message);
                hash_input.extend_from_slice(pk_coords.as_slice());
                let output = FH::evaluate(hash_input.as_ref())?;

                //Return VRF output
                Ok(output)
            }
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
    use algebra::{ToBytes, FromBytes, to_bytes, ToBits};
    use crate::{crh::{
        MNT4PoseidonHash, MNT6PoseidonHash, BLS12PoseidonHash,
        bowe_hopwood::BoweHopwoodPedersenCRH,
        pedersen::{
            PedersenCRH, PedersenWindow
        },
    }, vrf::{
        FieldBasedVrf,
        ecvrf::FieldBasedEcVrf,
    }, FixedLengthCRH};
    use rand::{Rng, thread_rng};

    #[derive(Clone)]
    struct TestWindow {}
    impl PedersenWindow for TestWindow {
        const WINDOW_SIZE: usize = 128;
        const NUM_WINDOWS: usize = 2;
    }

    #[derive(Clone)]
    struct JubJubTestWindow {}
    impl PedersenWindow for JubJubTestWindow {
        const WINDOW_SIZE: usize = 32;
        const NUM_WINDOWS: usize = 8;
    }

    type BHMNT4 = BoweHopwoodPedersenCRH<MNT4G1Projective, TestWindow>;
    type BHMNT6 = BoweHopwoodPedersenCRH<MNT6G1Projective, TestWindow>;
    type PedersenJubJub = PedersenCRH<JubJubProjective, JubJubTestWindow>;

    type EcVrfMNT4 = FieldBasedEcVrf<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash, BHMNT6>;
    type EcVrfMNT6 = FieldBasedEcVrf<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash, BHMNT4>;
    type EcVrfBLS12 = FieldBasedEcVrf<BLS12Fr, JubJubProjective, BLS12PoseidonHash, PedersenJubJub>;

    fn prove_and_verify<S: FieldBasedVrf, R: Rng>(rng: &mut R, message: &[S::Data], pp: &S::GHParams) {
        let (pk, sk) = S::keygen(rng).unwrap();
        let proof = S::prove(rng, pp, &pk, &sk, &message).unwrap();
        assert!(S::verify(pp, &pk, &message, &proof).is_ok());

        //Serialization/deserialization test
        let proof_serialized = to_bytes!(proof).unwrap();
        let proof_deserialized = <S as FieldBasedVrf>::Proof::read(proof_serialized.as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);
        assert!(S::verify(pp, &pk, &message, &proof_deserialized).is_ok());
    }

    fn failed_verification<S: FieldBasedVrf, R: Rng>(rng: &mut R, message: &[S::Data], bad_message: &[S::Data], pp: &S::GHParams) {
        let (pk, sk) = S::keygen(rng).unwrap();

        //Attempt to verify proof for a different message
        let proof = S::prove(rng, pp, &pk, &sk, message).unwrap();
        assert!(S::verify(pp, &pk, bad_message, &proof).is_err());

        //Attempt to verify different proof for a message
        let bad_proof = S::prove(rng, pp, &pk, &sk, bad_message).unwrap();
        assert!(S::verify(pp, &pk, message, &bad_proof).is_err());

        //Attempt to verify proof for a message with different pk
        let (new_pk, _) = S::keygen(rng).unwrap();
        assert!(S::verify(pp, &new_pk, message, &proof).is_err());
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

    #[test]
    fn bls12_381_ecvrf_test() {
        let rng = &mut thread_rng();
        let pp = <PedersenJubJub as FixedLengthCRH>::setup(rng).unwrap();
        let samples = 100;
        for _ in 0..samples {
            let f: BLS12Fr = rng.gen();
            let g: BLS12Fr = rng.gen();
            prove_and_verify::<EcVrfBLS12, _>(rng, &[f], &pp);
            failed_verification::<EcVrfBLS12, _>(rng, &[f], &[g], &pp);
        }
    }
}