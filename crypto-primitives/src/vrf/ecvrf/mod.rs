use algebra::{
    Field, PrimeField, project,
    Group, AffineCurve, ProjectiveCurve,
    ToBytes, to_bytes, FromBytes,
    UniformRand,
};
use crate::{
    crh::{
    FieldBasedHash, FixedLengthCRH,
    },
    vrf::FieldBasedVrf, Error, CryptoError
};
use std::marker::PhantomData;
use std::io::{Write, Read, Result as IoResult};
use rand::Rng;

#[cfg(feature = "r1cs")]
pub mod constraints;

pub struct FieldBasedEcVrf<
    F: PrimeField,
    G: ProjectiveCurve<BaseField = F>,
    FH: FieldBasedHash<Data = F>,
    GH: FixedLengthCRH<Output = G>,
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
pub struct FieldBasedEcVrfProof<F: PrimeField, G: ProjectiveCurve<BaseField = F>> {
    pub gamma:  G,
    pub c:      F,
    pub s:      Vec<u8>,
}

impl<F: PrimeField, G: ProjectiveCurve<BaseField = F>> ToBytes for FieldBasedEcVrfProof<F, G>
{
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.gamma.write(&mut writer)?;
        self.c.write(&mut writer)?;
        self.s.write(&mut writer)
    }
}

impl <F: PrimeField, G: ProjectiveCurve<BaseField = F>> FromBytes for FieldBasedEcVrfProof<F, G>
{
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let gamma = G::read(&mut reader)?;
        let c = F::read(&mut reader)?;
        let s = Vec::<u8>::read(&mut reader)?;
        Ok(Self{gamma, c, s})
    }
}

impl<F, G, FH, GH> FieldBasedVrf for FieldBasedEcVrf<F, G, FH, GH>
    where
        F: PrimeField,
        G: ProjectiveCurve<BaseField = F>,
        G::Affine: AffineCurve<BaseField = F>,
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

        //Choose random scalar
        let r = G::ScalarField::rand(rng);

        //Compute a = g^r
        let a = G::prime_subgroup_generator().mul(&r).into_affine();

        //Compute b = mh^r
        let b = mh.mul(&r).into_affine();

        //Compute c = H(m||pk||a||b)
        let pk = pk.into_affine();
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(message);
        hash_input.push(pk.get_x());
        hash_input.push(pk.get_y());
        hash_input.push(a.get_x());
        hash_input.push(a.get_y());
        hash_input.push(b.get_x());
        hash_input.push(b.get_y());
        let c = FH::evaluate(hash_input.as_ref())?;

        //Compute s = r + sk * c
        let s = r + &((*sk) * &(project::<F, G::ScalarField>(c)?));

        //Compute gamma = mh^sk
        let gamma = mh.mul(sk);

        Ok(FieldBasedEcVrfProof {gamma, c, s: to_bytes!(s).unwrap()})
    }

    fn verify(
        pp:      &Self::GHParams,
        pk:      &Self::PublicKey,
        message: &[Self::Data],
        proof:   &Self::Proof
    )
        -> Result<Self::Data, Error>
    {
        //Read s from proof
        let s = G::ScalarField::read(proof.s.as_slice())?;

        //Checks;
        let pk = pk.into_affine();
        let gamma = proof.gamma.into_affine();
        debug_assert!(pk.is_in_correct_subgroup_assuming_on_curve());
        debug_assert!(gamma.is_in_correct_subgroup_assuming_on_curve());
        debug_assert!(s.pow(&G::ScalarField::characteristic()) == s);
        debug_assert!(proof.c.pow(&G::BaseField::characteristic()) == proof.c);

        //Compute mh = hash_to_curve(message)
        let mut message_bytes = Vec::new();
        message.iter().for_each(|f| message_bytes.extend_from_slice(to_bytes!(f).unwrap().as_slice()));
        let mh = GH::evaluate(pp, message_bytes.as_slice())?;

        let c_projected = project::<F, G::ScalarField>(proof.c)?;

        //Compute u = g^s - pk^c
        let u = (G::prime_subgroup_generator().mul(&s) - &(pk.mul(c_projected))).into_affine();

        //Compute v = mh^s - gamma^c
        let v = (mh.mul(&s) - &gamma.mul(c_projected)).into_affine();

        //Compute c' = H(m||pk||u||v)
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(message);
        hash_input.push(pk.get_x());
        hash_input.push(pk.get_y());
        hash_input.push(u.get_x());
        hash_input.push(u.get_y());
        hash_input.push(v.get_x());
        hash_input.push(v.get_y());
        let c_prime = FH::evaluate(hash_input.as_ref())?;

        //Verify valid proof
        match proof.c == c_prime {
            false => Err(Box::new(CryptoError::FailedVerification)),
            true => {

                //Compute VRF output
                hash_input = Vec::new();
                hash_input.extend_from_slice(message);
                hash_input.push(pk.get_x());
                hash_input.push(pk.get_y());
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
    };
    use algebra::fields::{
        mnt4753::Fr as MNT4Fr,
        mnt6753::Fr as MNT6Fr,
    };
    use crate::{crh::{
        MNT4PoseidonHash, MNT6PoseidonHash,
        bowe_hopwood::BoweHopwoodPedersenCRH,
        pedersen::PedersenWindow,
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

    type BHMNT4 = BoweHopwoodPedersenCRH<MNT4G1Projective, TestWindow>;
    type BHMNT6 = BoweHopwoodPedersenCRH<MNT6G1Projective, TestWindow>;

    type EcVrfMNT4 = FieldBasedEcVrf<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash, BHMNT6>;
    type EcVrfMNT6 = FieldBasedEcVrf<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash, BHMNT4>;


    fn prove_and_verify<S: FieldBasedVrf, R: Rng>(rng: &mut R, message: &[S::Data], pp: &S::GHParams) {
        let (pk, sk) = S::keygen(rng).unwrap();
        let proof = S::prove(rng, pp, &pk, &sk, &message).unwrap();
        assert!(S::verify(pp, &pk, &message, &proof).is_ok());
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
}