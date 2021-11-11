use crate::{
    compute_truncation_size,
    crh::{FieldBasedHash, FixedLengthCRH},
    vrf::FieldBasedVrf,
    CryptoError, Error,
};
use algebra::{
    convert, leading_zeros, serialize::*, to_bytes, AffineCurve, Field, FromBytes,
    FromBytesChecked, Group, PrimeField, ProjectiveCurve, SemanticallyValid, ToBits, ToBytes,
    ToConstraintField, UniformRand,
};
use rand::distributions::{Distribution, Standard};
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::io::{self, Error as IoError, ErrorKind, Read, Result as IoResult, Write};
use std::marker::PhantomData;

pub struct FieldBasedEcVrf<F: PrimeField, G: Group, FH: FieldBasedHash, GH: FixedLengthCRH> {
    _field: PhantomData<F>,
    _group: PhantomData<G>,
    _field_hash: PhantomData<FH>,
    _group_hash: PhantomData<GH>,
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
#[derive(Serialize, Deserialize)]
#[serde(bound(serialize = "F: PrimeField, G: ProjectiveCurve"))]
#[serde(bound(deserialize = "F: PrimeField, G: ProjectiveCurve"))]
#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct FieldBasedEcVrfProof<F: PrimeField, G: ProjectiveCurve> {
    pub gamma: G,
    pub c: F,
    pub s: F,
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
        Ok(Self {
            gamma: gamma.into_projective(),
            c,
            s,
        })
    }
}

impl<F: PrimeField, G: ProjectiveCurve> FromBytesChecked for FieldBasedEcVrfProof<F, G> {
    fn read_checked<R: Read>(mut reader: R) -> IoResult<Self> {
        let gamma = G::Affine::read_checked(&mut reader)
            .map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid proof.gamma: {}", e),
                )
            })
            .and_then(|p| {
                if p.is_zero() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "invalid proof.gamma: point at infinity",
                    ));
                }
                Ok(p)
            })?;
        let c = F::read_checked(&mut reader)
            .map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid proof.c: {}", err),
                )
            })
            .and_then(|c| {
                let c_bits = c.write_bits();
                let c_leading_zeros = leading_zeros(c_bits.as_slice()) as usize;
                if (F::size_in_bits() - c_leading_zeros) >= G::ScalarField::size_in_bits() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "Invalid bit-length for proof.c: {}",
                            c_bits.len() - c_leading_zeros
                        ),
                    ));
                }
                Ok(c)
            })?;
        let s = F::read_checked(&mut reader)
            .map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid proof.s: {}", err),
                )
            })
            .and_then(|s| {
                let s_bits = s.write_bits();
                let s_leading_zeros = leading_zeros(s_bits.as_slice()) as usize;
                if (G::ScalarField::size_in_bits() - s_leading_zeros) >= F::size_in_bits() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "Invalid bit-length for proof.s: {}",
                            s_bits.len() - s_leading_zeros
                        ),
                    ));
                }
                Ok(s)
            })?;
        Ok(Self {
            gamma: gamma.into_projective(),
            c,
            s,
        })
    }
}

impl<F: PrimeField, G: ProjectiveCurve> SemanticallyValid for FieldBasedEcVrfProof<F, G> {
    fn is_valid(&self) -> bool {
        (self.gamma.is_valid() && !self.gamma.is_zero())
            && self.c.is_valid()
            && {
                //Checks c had proper bit-length when converted into a G::ScalarField element
                let c_bits = self.c.write_bits();
                let c_leading_zeros = leading_zeros(c_bits.as_slice()) as usize;
                F::size_in_bits() - c_leading_zeros < G::ScalarField::size_in_bits()
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
pub struct FieldBasedEcVrfPk<G: Group>(pub G);

impl<G: Group> Distribution<FieldBasedEcVrfPk<G>> for Standard {
    #[inline]
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> FieldBasedEcVrfPk<G> {
        let pk = G::rand(rng);
        FieldBasedEcVrfPk::<G>(pk)
    }
}

impl<G: Group> ToBytes for FieldBasedEcVrfPk<G> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.0.write(&mut writer)
    }
}

impl<G: Group> FromBytes for FieldBasedEcVrfPk<G> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let pk = G::read(&mut reader)?;
        Ok(Self(pk))
    }
}

impl<G: Group> FromBytesChecked for FieldBasedEcVrfPk<G> {
    fn read_checked<R: Read>(mut reader: R) -> IoResult<Self> {
        let pk = G::read_checked(&mut reader)
            .map_err(|e| IoError::new(ErrorKind::InvalidData, format!("invalid ecvrf pk: {}", e)))
            .and_then(|p| {
                if p.is_zero() {
                    return Err(IoError::new(
                        ErrorKind::InvalidData,
                        "invalid ecvrf pk: point at infinity",
                    ));
                }
                Ok(p)
            })?;
        Ok(Self(pk))
    }
}

impl<G: Group> SemanticallyValid for FieldBasedEcVrfPk<G> {
    #[inline]
    fn is_valid(&self) -> bool {
        self.0.is_valid() &&
            // GingerLib only accepts non-trivial ECVRF public keys. This is usually
            // good practice to avoid using obvious weak keys, and helps preventing
            // exceptional cases if using incomplete arithmetics.
            !self.0.is_zero()
    }
}

// Low-level crypto for our length-restricted variant of the DL-based VRF, does not perform any
// input validity check. It's responsibility of the caller to do so, through keyverify()
// function for the PublicKey, read() or is_valid() functions for FieldBasedEcVrfProof.
impl<F, G, FH, GH> FieldBasedVrf for FieldBasedEcVrf<F, G, FH, GH>
where
    F: PrimeField,
    G: ProjectiveCurve + ToConstraintField<F>,
    FH: FieldBasedHash<Data = F>,
    GH: FixedLengthCRH<Output = G>,
{
    type Data = FH::Data;
    type PublicKey = FieldBasedEcVrfPk<G>;
    type SecretKey = G::ScalarField;
    type Proof = FieldBasedEcVrfProof<F, G>;
    type GHParams = GH::Parameters;

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
        (FieldBasedEcVrfPk(public_key), secret_key)
    }

    fn get_public_key(sk: &Self::SecretKey) -> Self::PublicKey {
        FieldBasedEcVrfPk(G::prime_subgroup_generator().mul(sk))
    }

    fn prove<R: Rng>(
        rng: &mut R,
        group_hash_params: &Self::GHParams,
        pk: &Self::PublicKey,
        sk: &Self::SecretKey,
        message: Self::Data,
    ) -> Result<Self::Proof, Error> {
        //Compute mh = hash_to_curve(message)
        let message_on_curve =
            GH::evaluate(group_hash_params, to_bytes!(&message).unwrap().as_slice())?;

        //Compute gamma = message_on_curve^sk
        let gamma = message_on_curve.mul(sk);

        let required_leading_zeros_c = compute_truncation_size(
            F::size_in_bits() as i32,
            G::ScalarField::size_in_bits() as i32,
        );

        let required_leading_zeros_s = compute_truncation_size(
            G::ScalarField::size_in_bits() as i32,
            F::size_in_bits() as i32,
        );

        let (c, s) = loop {
            //Choose random scalar
            let r = G::ScalarField::rand(rng);

            //Compute a = g^r
            let a = G::prime_subgroup_generator().mul(&r);

            //Compute b = message_on_curve^r
            let b = message_on_curve.mul(&r);

            //Compute c = H(m||pk.x||a.x||b.x)
            let c = {
                let mut digest = FH::init_constant_length(4, None);
                digest
                    .update(message)
                    .update(pk.0.to_field_elements()?[0])
                    .update(a.to_field_elements()?[0])
                    .update(b.to_field_elements()?[0])
                    .finalize()
            }?;

            let c_bits = c.write_bits();
            let c_leading_zeros = leading_zeros(c_bits.as_slice()) as usize;

            //Enforce c bit length is strictly smaller than G::ScalarField modulus bit length
            if c_leading_zeros < required_leading_zeros_c {
                continue;
            };

            let c_conv = convert::<G::ScalarField>(c_bits)?;

            //Compute s = r + sk * c
            let s = r + &((*sk) * &c_conv);
            let s_bits = s.write_bits();
            let s_leading_zeros = leading_zeros(s_bits.as_slice()) as usize;

            if s_leading_zeros < required_leading_zeros_s {
                continue;
            };

            let s_conv = convert::<F>(s_bits)?;

            break (c, s_conv);
        };

        Ok(FieldBasedEcVrfProof { gamma, c, s })
    }

    fn proof_to_hash(
        group_hash_params: &Self::GHParams,
        pk: &Self::PublicKey,
        message: Self::Data,
        proof: &Self::Proof,
    ) -> Result<Self::Data, Error> {
        //Compute mh = hash_to_curve(message)
        let message_on_curve =
            GH::evaluate(group_hash_params, to_bytes!(&message).unwrap().as_slice())?;

        let c_bits = proof.c.write_bits();
        let s_bits = proof.s.write_bits();
        let c_conv = convert::<G::ScalarField>(c_bits)?;
        let s_conv = convert::<G::ScalarField>(s_bits)?;

        //Compute u = g^s - pk^c
        let u = G::prime_subgroup_generator().mul(&s_conv) - &(pk.0.mul(&c_conv));

        //Compute v = mh^s - gamma^c
        let v = message_on_curve.mul(&s_conv) - &proof.gamma.mul(&c_conv);

        //Compute c' = H(m||pk.x||u.x||v.x)
        let c_prime = {
            let mut digest = FH::init_constant_length(4, None);

            digest
                .update(message.clone())
                .update(pk.0.to_field_elements()?[0])
                .update(u.to_field_elements()?[0])
                .update(v.to_field_elements()?[0])
                .finalize()
        }?;

        //Verify valid proof
        match proof.c == c_prime {
            false => Err(Box::new(CryptoError::FailedVerification)),
            true => {
                let gamma_coords = proof.gamma.to_field_elements()?;

                //Compute VRF output
                let output = {
                    let mut digest = FH::init_constant_length(3, None);
                    digest.update(message);
                    gamma_coords.into_iter().for_each(|c| {
                        digest.update(c);
                    });
                    digest.finalize()
                }?;

                //Return VRF output
                Ok(output)
            }
        }
    }

    fn keyverify(pk: &Self::PublicKey) -> bool {
        pk.is_valid()
    }
}

#[cfg(test)]
mod test {
    use crate::{
        crh::{
            bowe_hopwood::BoweHopwoodPedersenCRH, pedersen::PedersenWindow, MNT4PoseidonHash,
            MNT6PoseidonHash,
        },
        vrf::{ecvrf::FieldBasedEcVrf, FieldBasedVrf},
        FixedLengthCRH,
    };
    use algebra::curves::{
        mnt4753::G1Projective as MNT4G1Projective, mnt6753::G1Projective as MNT6G1Projective,
    };
    use algebra::fields::{mnt4753::Fr as MNT4Fr, mnt6753::Fr as MNT6Fr};
    use algebra::{to_bytes, FromBytes, FromBytesChecked, SemanticallyValid, ToBytes};
    use rand::{thread_rng, Rng};

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

    fn prove_and_verify<S: FieldBasedVrf, R: Rng>(rng: &mut R, message: S::Data, pp: &S::GHParams) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));
        assert_eq!(pk, S::get_public_key(&sk));

        let proof = S::prove(rng, pp, &pk, &sk, message).unwrap();
        assert!(proof.is_valid());
        assert!(S::proof_to_hash(pp, &pk, message, &proof).is_ok());
    }

    fn failed_verification<S: FieldBasedVrf, R: Rng>(
        rng: &mut R,
        message: S::Data,
        bad_message: S::Data,
        pp: &S::GHParams,
    ) {
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

    fn serialize_deserialize<S: FieldBasedVrf, R: Rng>(
        rng: &mut R,
        message: S::Data,
        pp: &S::GHParams,
    ) {
        let (pk, sk) = S::keygen(rng);
        let proof = S::prove(rng, pp, &pk, &sk, message).unwrap();

        let proof_serialized = to_bytes!(proof).unwrap();

        let proof_deserialized =
            <S as FieldBasedVrf>::Proof::read(proof_serialized.as_slice()).unwrap();
        assert_eq!(proof, proof_deserialized);
        assert!(<S as FieldBasedVrf>::Proof::read_checked(proof_serialized.as_slice()).is_ok());
        assert!(S::proof_to_hash(pp, &pk, message, &proof_deserialized).is_ok());
    }

    #[test]
    fn mnt4_ecvrf_test() {
        let rng = &mut thread_rng();
        let pp = <BHMNT6 as FixedLengthCRH>::setup(rng).unwrap();
        let samples = 100;
        for _ in 0..samples {
            let f: MNT4Fr = rng.gen();
            let g: MNT4Fr = rng.gen();
            prove_and_verify::<EcVrfMNT4, _>(rng, f, &pp);
            failed_verification::<EcVrfMNT4, _>(rng, f, g, &pp);
            serialize_deserialize::<EcVrfMNT4, _>(rng, f, &pp);
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
            prove_and_verify::<EcVrfMNT6, _>(rng, f, &pp);
            failed_verification::<EcVrfMNT6, _>(rng, f, g, &pp);
            serialize_deserialize::<EcVrfMNT6, _>(rng, f, &pp);
        }
    }
}
