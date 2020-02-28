use algebra::{PrimeField, FpParameters, ProjectiveCurve, Group};
use r1cs_std::{
    fields::fp::FpGadget,
    alloc::AllocGadget,
    groups::GroupGadget,
    eq::EqGadget,
    ToBytesGadget,
};
use crate::{
    vrf::{
        ecvrf::{FieldBasedEcVrfProof, FieldBasedEcVrf},
        FieldBasedVrfGadget,
    },
    crh::{
        FixedLengthCRH, FixedLengthCRHGadget,
        FieldBasedHash, FieldBasedHashGadget,
    },
};
use r1cs_core::{ConstraintSystem, SynthesisError, ToConstraintField};
use std::{
    marker::PhantomData,
    borrow::Borrow,
};
use r1cs_std::to_field_gadget_vec::ToConstraintFieldGadget;

#[derive(Derivative)]
#[derivative(
Debug(bound = "ConstraintF: PrimeField, G: Group, GG: GroupGadget<G, ConstraintF>"),
Clone(bound = "ConstraintF: PrimeField, G: Group, GG: GroupGadget<G, ConstraintF>"),
PartialEq(bound = "ConstraintF: PrimeField, G: Group, GG: GroupGadget<G, ConstraintF>"),
Eq(bound = "ConstraintF: PrimeField, G: Group, GG: GroupGadget<G, ConstraintF>")
)]
pub struct FieldBasedEcVrfProofGadget<ConstraintF, G, GG>
where
    ConstraintF: PrimeField,
    G:           Group,
    GG:          GroupGadget<G, ConstraintF>,
{
    pub gamma:   GG,
    pub c:       FpGadget<ConstraintF>,
    pub s:       FpGadget<ConstraintF>,
    _field:      PhantomData<ConstraintF>,
    _group:      PhantomData<G>,
}

impl<ConstraintF, G, GG> AllocGadget<FieldBasedEcVrfProof<ConstraintF, G>, ConstraintF>
for FieldBasedEcVrfProofGadget<ConstraintF, G, GG>
    where
        ConstraintF: PrimeField,
        G:           ProjectiveCurve,
        GG:          GroupGadget<G, ConstraintF>,
{
    fn alloc<FN, T, CS: ConstraintSystem<ConstraintF>>(mut cs: CS, f: FN) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<FieldBasedEcVrfProof<ConstraintF, G>>,
    {
        f().and_then(|proof| {
            let FieldBasedEcVrfProof {
                gamma,
                c,
                s
            } = proof.borrow().clone();
            let gamma = GG::alloc(cs.ns(|| "alloc gamma"), || Ok(gamma))?;
            let c = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc r"), || Ok(c))?;
            let s = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc s"), || Ok(s))?;
            Ok(Self{gamma, c, s, _field: PhantomData, _group: PhantomData})
        })
    }

    fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(mut cs: CS, f: FN) -> Result<Self, SynthesisError>
        where
            FN: FnOnce() -> Result<T, SynthesisError>,
            T: Borrow<FieldBasedEcVrfProof<ConstraintF, G>>,
    {
        f().and_then(|proof| {
            let FieldBasedEcVrfProof {
                gamma,
                c,
                s
            } = proof.borrow().clone();
            let gamma = GG::alloc_input(cs.ns(|| "alloc gamma"), || Ok(gamma))?;
            let c = FpGadget::<ConstraintF>::alloc_input(cs.ns(|| "alloc r"), || Ok(c))?;
            let s = FpGadget::<ConstraintF>::alloc_input(cs.ns(|| "alloc s"), || Ok(s))?;
            Ok(Self{gamma, c, s, _field: PhantomData, _group: PhantomData})
        })
    }
}

pub struct FieldBasedEcVrfProofVerificationGadget<
    ConstraintF: PrimeField,
    G:  ProjectiveCurve,
    GG: GroupGadget<G, ConstraintF>,
    FH:  FieldBasedHash<Data = ConstraintF>,
    FHG: FieldBasedHashGadget<FH, ConstraintF>,
    GH:  FixedLengthCRH<Output = G>,
    GHG: FixedLengthCRHGadget<GH, ConstraintF, OutputGadget = GG>,
>
where
{
    _field:               PhantomData<ConstraintF>,
    _group:               PhantomData<G>,
    _group_gadget:        PhantomData<GG>,
    _field_hash:          PhantomData<FH>,
    _field_hash_gadget:   PhantomData<FHG>,
    _group_hash:          PhantomData<GH>,
    _group_hash_gadget:   PhantomData<GHG>,
}

impl<ConstraintF, G, GG, FH, FHG, GH, GHG> FieldBasedVrfGadget<FieldBasedEcVrf<ConstraintF, G, FH, GH>, ConstraintF>
for FieldBasedEcVrfProofVerificationGadget<ConstraintF, G, GG, FH, FHG, GH, GHG>
    where
        ConstraintF: PrimeField,
        G:           ProjectiveCurve + ToConstraintField<ConstraintF>,
        GG:          GroupGadget<G, ConstraintF, Value = G> + ToConstraintFieldGadget<ConstraintF, FieldGadget = FHG::DataGadget>,
        FH:          FieldBasedHash<Data = ConstraintF>,
        FHG:         FieldBasedHashGadget<FH, ConstraintF, DataGadget = FpGadget<ConstraintF>>,
        GH:          FixedLengthCRH<Output = G>,
        GHG:         FixedLengthCRHGadget<GH, ConstraintF, OutputGadget = GG>,
{
    type DataGadget = FpGadget<ConstraintF>;
    type ProofGadget = FieldBasedEcVrfProofGadget<ConstraintF, G, GG>;
    type PublicKeyGadget = GG;
    type GHParametersGadget = GHG::ParametersGadget;

    fn check_verify_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs:     CS,
        pp:         &Self::GHParametersGadget,
        public_key: &Self::PublicKeyGadget,
        proof:      &Self::ProofGadget,
        message:    &[Self::DataGadget]
    ) -> Result<Self::DataGadget, SynthesisError> {

        //Check mh = hash_to_curve(message)
        let mut message_bytes = Vec::new();
        for (i, fg) in message.iter().enumerate() {
            let fg_bytes = fg.to_bytes(
                cs.ns(|| format!("message_{}_to_bytes", i))
            )?;
            message_bytes.extend_from_slice(fg_bytes.as_slice())
        }

        let mh = GHG::check_evaluation_gadget(
            cs.ns(|| "check mh"),
            pp,
            message_bytes.as_slice()
        )?;

        //Hardcode g, serialize c and s
        let g = GG::alloc_hardcoded(
            cs.ns(|| "hardcode generator"),
            || Ok(G::prime_subgroup_generator())
        )?;

        let c_bits = {

            //Serialize e taking into account the length restriction
            let mut to_skip = 0usize;
            let moduli_diff = ConstraintF::Params::MODULUS_BITS as i32 -
                <G::ScalarField as PrimeField>::Params::MODULUS_BITS as i32;
            if moduli_diff >= 0 {
                to_skip = moduli_diff as usize + 1;
            }

            let c_bits = proof.c
                .to_bits_with_length_restriction(cs.ns(|| "proof.c_to_bits"), to_skip)?;
            debug_assert!(c_bits.len() as u32 == ConstraintF::Params::MODULUS_BITS - to_skip as u32);
            c_bits
        };

        let mut s_bits = {

            //Serialize s taking into account the length restriction
            let mut to_skip = 0usize;
            let moduli_diff = <G::ScalarField as PrimeField>::Params::MODULUS_BITS as i32 -
                ConstraintF::Params::MODULUS_BITS as i32;
            if moduli_diff >= 0 {
                to_skip = moduli_diff as usize + 1;
            }

            let s_bits = proof.s
                .to_bits_with_length_restriction(cs.ns(|| "proof.s_to_bits"), to_skip)?;
            debug_assert!(s_bits.len() as u32 == <G::ScalarField as PrimeField>::Params::MODULUS_BITS - to_skip as u32);
            s_bits
        };
        s_bits.reverse();

        //Check u = g^s - pk^c
        let u =
        {
            let neg_c_times_pk = public_key
                .mul_bits(cs.ns(|| "pk * c + g"), &g, c_bits.as_slice().iter().rev())?
                .sub(cs.ns(|| "c * pk"), &g)?
                .negate(cs.ns(|| "- (c * pk)"))?;
            GG::mul_bits_precomputed(&(g.get_value().unwrap()),
                                     cs.ns(|| "(s * G) - (c * pk)"),
                                     &neg_c_times_pk,
                                     s_bits.as_slice()
            )?
        };

        //Check v = mh^s - gamma^c
        let v =
        {
            let neg_c_times_gamma = proof.gamma
                .mul_bits(cs.ns(|| "c * gamma + g"), &g, c_bits.as_slice().iter().rev())?
                .sub(cs.ns(|| "c * gamma"), &g)?
                .negate(cs.ns(|| "- (c * gamma)"))?;
            mh.mul_bits(cs.ns(|| "(s * mh) - (c * gamma)"), &neg_c_times_gamma, s_bits.as_slice().iter())?
        };

        //Check c' = H(m||pk||u||v)
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(message);
        hash_input.extend_from_slice(public_key.to_field_gadget_elements().unwrap().as_slice());
        hash_input.extend_from_slice(u.to_field_gadget_elements().unwrap().as_slice());
        hash_input.extend_from_slice(v.to_field_gadget_elements().unwrap().as_slice());

        let c_prime = FHG::check_evaluation_gadget(
            cs.ns(|| "check c_prime"),
            hash_input.as_slice()
        )?;

        //Enforce c = c'
        proof.c.enforce_equal(cs.ns(|| "check c == c'"), &c_prime)?;

        //Check and return VRF output
        hash_input = Vec::new();
        hash_input.extend_from_slice(message);
        hash_input.extend_from_slice(public_key.to_field_gadget_elements().unwrap().as_slice());

        let vrf_output = FHG::check_evaluation_gadget(
            cs.ns(|| "check vrf_output"),
            hash_input.as_slice()
        )?;

        Ok(vrf_output)
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
    use crate::{vrf::{
        FieldBasedVrf, FieldBasedVrfGadget,
        ecvrf::{FieldBasedEcVrf, constraints::FieldBasedEcVrfProofVerificationGadget},
    }, crh::{
        MNT4PoseidonHash, MNT6PoseidonHash, MNT4PoseidonHashGadget, MNT6PoseidonHashGadget,
        bowe_hopwood::{
            BoweHopwoodPedersenCRH,
            constraints::BoweHopwoodPedersenCRHGadget,
        },
        pedersen::PedersenWindow,
    }, FixedLengthCRH};

    use r1cs_core::ConstraintSystem;
    use r1cs_std::alloc::AllocGadget;

    use r1cs_std::groups::curves::short_weierstrass::mnt::{
        mnt4::mnt4753::MNT4G1Gadget,
        mnt6::mnt6753::MNT6G1Gadget,
    };

    use rand::{Rng, thread_rng};
    use r1cs_std::test_constraint_system::TestConstraintSystem;

    #[derive(Clone)]
    struct TestWindow {}
    impl PedersenWindow for TestWindow {
        const WINDOW_SIZE: usize = 128;
        const NUM_WINDOWS: usize = 2;
    }

    type BHMNT4 = BoweHopwoodPedersenCRH<MNT4G1Projective, TestWindow>;
    type BHMNT6 = BoweHopwoodPedersenCRH<MNT6G1Projective, TestWindow>;

    type BHMNT4Gadget = BoweHopwoodPedersenCRHGadget<MNT6G1Projective, MNT4Fr, MNT6G1Gadget>;
    type BHMNT6Gadget = BoweHopwoodPedersenCRHGadget<MNT4G1Projective, MNT6Fr, MNT4G1Gadget>;

    type EcVrfMNT4 = FieldBasedEcVrf<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash, BHMNT6>;
    type EcVrfMNT6 = FieldBasedEcVrf<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash, BHMNT4>;

    type EcVrfMNT4Gadget = FieldBasedEcVrfProofVerificationGadget<
        MNT4Fr,
        MNT6G1Projective,
        MNT6G1Gadget,
        MNT4PoseidonHash,
        MNT4PoseidonHashGadget,
        BHMNT6,
        BHMNT4Gadget>;

    type EcVrfMNT6Gadget = FieldBasedEcVrfProofVerificationGadget<
        MNT6Fr,
        MNT4G1Projective,
        MNT4G1Gadget,
        MNT6PoseidonHash,
        MNT6PoseidonHashGadget,
        BHMNT4,
        BHMNT6Gadget>;

    fn prove<S: FieldBasedVrf, R: Rng>(rng: &mut R, pp: &S::GHParams, message: &[S::Data])
        -> (S::Proof, S::PublicKey, S::SecretKey)
    {
        let (pk, sk) = S::keygen(rng).unwrap();
        let proof = S::prove(rng, pp, &pk, &sk, &message).unwrap();
        (proof, pk, sk)
    }

    #[test]
    fn mnt4_ecvrf_gadget_test() {

        let mut cs = TestConstraintSystem::<MNT4Fr>::new();

        //Generate VRF proof for a random field element f and get the proof and the keypair too
        let rng = &mut thread_rng();
        let message: MNT4Fr = rng.gen();
        let pp = <BHMNT6 as FixedLengthCRH>::setup(rng).unwrap();
        let (proof, pk, sk) = prove::<EcVrfMNT4, _>(rng, &pp, &[message]);

        //Alloc proof, pk and message
        let proof_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::ProofGadget::alloc(
            cs.ns(|| "alloc proof"),
            || Ok(proof)
        ).unwrap();

        let pk_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::PublicKeyGadget::alloc(cs.ns(|| "alloc pk"), || Ok(pk)).unwrap();

        let pp_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::GHParametersGadget::alloc(cs.ns(|| "alloc gh params"), || Ok(&pp)).unwrap();

        let message_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::DataGadget::alloc(
            cs.ns(|| "alloc message"),
            || Ok(message)
        ).unwrap();

        //Verify proof
        EcVrfMNT4Gadget::check_verify_gadget(
            cs.ns(|| "verify proof1"),
            &pp_g,
            &pk_g,
            &proof_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(cs.is_satisfied());

        //Wrong proof: generate a proof for a different message and check constraints fail
        let new_message: MNT4Fr = rng.gen();
        let new_proof = EcVrfMNT4::prove(rng, &pp, &pk, &sk, &[new_message]).unwrap();
        let new_proof_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::ProofGadget::alloc(
            cs.ns(|| "alloc new proof"),
            || Ok(new_proof)
        ).unwrap();

        //Verify new proof: expected to fail
        EcVrfMNT4Gadget::check_verify_gadget(
            cs.ns(|| "verify proof2"),
            &pp_g,
            &pk_g,
            &new_proof_g,
            &[message_g]
        ).unwrap();

        assert!(!cs.is_satisfied());
        println!("{:?}", cs.which_is_unsatisfied());
    }

    #[test]
    fn mnt6_ecvrf_gadget_test() {

        let mut cs = TestConstraintSystem::<MNT6Fr>::new();

        //Generate VRF proof for a random field element f and get the proof and the keypair too
        let rng = &mut thread_rng();
        let message: MNT6Fr = rng.gen();
        let pp = <BHMNT4 as FixedLengthCRH>::setup(rng).unwrap();
        let (proof, pk, sk) = prove::<EcVrfMNT6, _>(rng, &pp, &[message]);

        //Alloc proof, pk and message
        let proof_g = <EcVrfMNT6Gadget as FieldBasedVrfGadget<EcVrfMNT6, MNT6Fr>>::ProofGadget::alloc(
            cs.ns(|| "alloc proof"),
            || Ok(proof)
        ).unwrap();
        let pk_g = <EcVrfMNT6Gadget as FieldBasedVrfGadget<EcVrfMNT6, MNT6Fr>>::PublicKeyGadget::alloc(cs.ns(|| "alloc pk"), || Ok(pk)).unwrap();
        let pp_g = <EcVrfMNT6Gadget as FieldBasedVrfGadget<EcVrfMNT6, MNT6Fr>>::GHParametersGadget::alloc(cs.ns(|| "alloc gh params"), || Ok(&pp)).unwrap();
        let message_g = <EcVrfMNT6Gadget as FieldBasedVrfGadget<EcVrfMNT6, MNT6Fr>>::DataGadget::alloc(
            cs.ns(|| "alloc message"),
            || Ok(message)
        ).unwrap();

        //Verify proof
        EcVrfMNT6Gadget::check_verify_gadget(
            cs.ns(|| "verify proof1"),
            &pp_g,
            &pk_g,
            &proof_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(cs.is_satisfied());

        //Wrong proof: generate a proof for a different message and check constraints fail
        let new_message: MNT6Fr = rng.gen();
        let new_proof = EcVrfMNT6::prove(rng, &pp, &pk, &sk, &[new_message]).unwrap();
        let new_proof_g = <EcVrfMNT6Gadget as FieldBasedVrfGadget<EcVrfMNT6, MNT6Fr>>::ProofGadget::alloc(
            cs.ns(|| "alloc new proof"),
            || Ok(new_proof)
        ).unwrap();

        //Verify new proof: expected to fail
        EcVrfMNT6Gadget::check_verify_gadget(
            cs.ns(|| "verify proof2"),
            &pp_g,
            &pk_g,
            &new_proof_g,
            &[message_g]
        ).unwrap();

        assert!(!cs.is_satisfied());
        println!("{:?}", cs.which_is_unsatisfied());
    }

    #[test]
    fn random_ecvrf_gadget_test() {

        //Generate VRF proof for a random field element f and get the proof and the public key
        let rng = &mut thread_rng();
        let pp = <BHMNT6 as FixedLengthCRH>::setup(rng).unwrap();

        let samples = 10;
        for _ in 0..samples {
            let message: MNT4Fr = rng.gen();
            let (sig, pk, _) = prove::<EcVrfMNT4, _>(rng, &pp, &[message]);
            let mut cs = TestConstraintSystem::<MNT4Fr>::new();

            //Alloc proof, pk, hash params and message
            let proof_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::ProofGadget::alloc(
                cs.ns(|| "alloc proof"),
                || Ok(sig)
            ).unwrap();

            let pk_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::PublicKeyGadget::alloc(
                cs.ns(|| "alloc pk"),
                || Ok(pk)
            ).unwrap();

            let pp_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::GHParametersGadget::alloc(
                cs.ns(|| "alloc gh params"),
                || Ok(&pp)
            ).unwrap();

            let message_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::DataGadget::alloc(
                cs.ns(|| "alloc message"),
                || Ok(message)
            ).unwrap();

            //Verify sig
            EcVrfMNT4Gadget::check_verify_gadget(
                cs.ns(|| "verify proof"),
                &pp_g,
                &pk_g,
                &proof_g,
                &[message_g.clone()]
            ).unwrap();
            assert!(cs.is_satisfied());
        }
    }
}