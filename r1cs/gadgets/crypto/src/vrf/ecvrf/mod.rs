use algebra::{Group, PrimeField, ProjectiveCurve};
use r1cs_std::{
    alloc::AllocGadget, bits::ToBytesGadget, eq::EqGadget, fields::fp::FpGadget,
    groups::GroupGadget, to_field_gadget_vec::ToConstraintFieldGadget,
};

use primitives::{
    compute_truncation_size,
    crh::{FieldBasedHash, FixedLengthCRH},
    vrf::ecvrf::{FieldBasedEcVrf, FieldBasedEcVrfProof},
};

use crate::{
    crh::{FieldBasedHashGadget, FixedLengthCRHGadget},
    vrf::FieldBasedVrfGadget,
};
use primitives::vrf::ecvrf::FieldBasedEcVrfPk;
use r1cs_core::{ConstraintSystem, SynthesisError, ToConstraintField};
use r1cs_std::bits::boolean::Boolean;
use rand::rngs::OsRng;
use std::{borrow::Borrow, marker::PhantomData};

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
    G: Group,
    GG: GroupGadget<G, ConstraintF>,
{
    pub gamma: GG,
    pub c: FpGadget<ConstraintF>,
    pub s: FpGadget<ConstraintF>,
    _field: PhantomData<ConstraintF>,
    _group: PhantomData<G>,
}

impl<ConstraintF, G, GG> FieldBasedEcVrfProofGadget<ConstraintF, G, GG>
where
    ConstraintF: PrimeField,
    G: ProjectiveCurve,
    GG: GroupGadget<G, ConstraintF>,
{
    fn alloc_internal<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: FN,
        gamma_on_curve: bool,
        gamma_prime_order: bool,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedEcVrfProof<ConstraintF, G>>,
    {
        let (gamma, c, s) = match f() {
            Ok(proof) => {
                let proof = *proof.borrow();
                (Ok(proof.gamma), Ok(proof.c), Ok(proof.s))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let gamma = match (gamma_on_curve, gamma_prime_order) {
            (false, false) => GG::alloc_without_check(cs.ns(|| "alloc gamma unchecked"), || gamma)?,
            (true, false) => GG::alloc(cs.ns(|| "alloc gamma"), || gamma).and_then(|gamma_g| {
                gamma_g.is_zero(cs.ns(|| "is gamma zero"))?.enforce_equal(
                    cs.ns(|| "gamma must not be zero"),
                    &Boolean::constant(false),
                )?;
                Ok(gamma_g)
            })?,
            (true, true) => GG::alloc_checked(cs.ns(|| "alloc gamma checked"), || gamma).and_then(
                |gamma_g| {
                    gamma_g.is_zero(cs.ns(|| "is gamma zero"))?.enforce_equal(
                        cs.ns(|| "gamma must not be zero"),
                        &Boolean::constant(false),
                    )?;
                    Ok(gamma_g)
                },
            )?,
            _ => unreachable!(),
        };
        let c = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc c"), || c)?;
        let s = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc s"), || s)?;
        Ok(Self {
            gamma,
            c,
            s,
            _field: PhantomData,
            _group: PhantomData,
        })
    }
}

impl<ConstraintF, G, GG> AllocGadget<FieldBasedEcVrfProof<ConstraintF, G>, ConstraintF>
    for FieldBasedEcVrfProofGadget<ConstraintF, G, GG>
where
    ConstraintF: PrimeField,
    G: ProjectiveCurve,
    GG: GroupGadget<G, ConstraintF>,
{
    fn alloc_without_check<FN, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedEcVrfProof<ConstraintF, G>>,
    {
        Self::alloc_internal(cs, f, false, false)
    }

    fn alloc<FN, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedEcVrfProof<ConstraintF, G>>,
    {
        Self::alloc_internal(cs, f, true, false)
    }

    fn alloc_checked<FN, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedEcVrfProof<ConstraintF, G>>,
    {
        Self::alloc_internal(cs, f, true, true)
    }

    fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedEcVrfProof<ConstraintF, G>>,
    {
        let (gamma, c, s) = match f() {
            Ok(proof) => {
                let proof = *proof.borrow();
                (Ok(proof.gamma), Ok(proof.c), Ok(proof.s))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let gamma = GG::alloc_input(cs.ns(|| "alloc gamma"), || gamma)?;
        let c = FpGadget::<ConstraintF>::alloc_input(cs.ns(|| "alloc c"), || c)?;
        let s = FpGadget::<ConstraintF>::alloc_input(cs.ns(|| "alloc s"), || s)?;
        Ok(Self {
            gamma,
            c,
            s,
            _field: PhantomData,
            _group: PhantomData,
        })
    }
}

pub struct FieldBasedEcVrfPkGadget<
    ConstraintF: PrimeField,
    G: Group,
    GG: GroupGadget<G, ConstraintF>,
> {
    pub pk: GG,
    _field: PhantomData<ConstraintF>,
    _group: PhantomData<G>,
}

impl<ConstraintF, G, GG> AllocGadget<FieldBasedEcVrfPk<G>, ConstraintF>
    for FieldBasedEcVrfPkGadget<ConstraintF, G, GG>
where
    ConstraintF: PrimeField,
    G: Group,
    GG: GroupGadget<G, ConstraintF>,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedEcVrfPk<G>>,
    {
        let pk =
            GG::alloc(cs.ns(|| "alloc pk"), || f().map(|pk| pk.borrow().0)).and_then(|pk_g| {
                pk_g.is_zero(cs.ns(|| "is pk zero"))?
                    .enforce_equal(cs.ns(|| "pk must not be zero"), &Boolean::constant(false))?;
                Ok(pk_g)
            })?;
        Ok(Self {
            pk,
            _field: PhantomData,
            _group: PhantomData,
        })
    }

    fn alloc_without_check<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedEcVrfPk<G>>,
    {
        let pk = GG::alloc_without_check(cs.ns(|| "alloc pk"), || f().map(|pk| pk.borrow().0))?;
        Ok(Self {
            pk,
            _field: PhantomData,
            _group: PhantomData,
        })
    }

    fn alloc_checked<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedEcVrfPk<G>>,
    {
        let pk = GG::alloc_checked(cs.ns(|| "alloc pk checked"), || f().map(|pk| pk.borrow().0))
            .and_then(|pk_g| {
                pk_g.is_zero(cs.ns(|| "is pk zero"))?
                    .enforce_equal(cs.ns(|| "pk must not be zero"), &Boolean::constant(false))?;
                Ok(pk_g)
            })?;
        Ok(Self {
            pk,
            _field: PhantomData,
            _group: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedEcVrfPk<G>>,
    {
        let pk = GG::alloc_input(cs.ns(|| "alloc pk"), || f().map(|pk| pk.borrow().0))?;
        Ok(Self {
            pk,
            _field: PhantomData,
            _group: PhantomData,
        })
    }
}

pub struct FieldBasedEcVrfProofVerificationGadget<
    ConstraintF: PrimeField,
    G: ProjectiveCurve,
    GG: GroupGadget<G, ConstraintF>,
    FH: FieldBasedHash<Data = ConstraintF>,
    FHG: FieldBasedHashGadget<FH, ConstraintF>,
    GH: FixedLengthCRH<Output = G>,
    GHG: FixedLengthCRHGadget<GH, ConstraintF, OutputGadget = GG>,
> {
    _field: PhantomData<ConstraintF>,
    _group: PhantomData<G>,
    _group_gadget: PhantomData<GG>,
    _field_hash: PhantomData<FH>,
    _field_hash_gadget: PhantomData<FHG>,
    _group_hash: PhantomData<GH>,
    _group_hash_gadget: PhantomData<GHG>,
}

// This implementation supports both complete and incomplete (safe) point addition.
// Assumes provided key material to be already checked.
//
// In case of incomplete point addition, with negligible probability, the
// proof creation might fail at first attempt and must be re-run (in order to sample
// fresh randomnesses).
// Furthermore, two exceptional cases (gamma, c, s) have to be treated outside the circuit:
// - if c * pk = s * G, i.e. when u is trivial (therefore leaking the sk), OR
// - if c * gamma = s * mh, i.e. when v is trivial (therefore also leaking the sk), THEN
// the circuit is not satisfiable.
impl<ConstraintF, G, GG, FH, FHG, GH, GHG>
    FieldBasedVrfGadget<FieldBasedEcVrf<ConstraintF, G, FH, GH>, ConstraintF>
    for FieldBasedEcVrfProofVerificationGadget<ConstraintF, G, GG, FH, FHG, GH, GHG>
where
    ConstraintF: PrimeField,
    G: ProjectiveCurve + ToConstraintField<ConstraintF>,
    GG: GroupGadget<G, ConstraintF, Value = G>
        + ToConstraintFieldGadget<ConstraintF, FieldGadget = FHG::DataGadget>,
    FH: FieldBasedHash<Data = ConstraintF>,
    FHG: FieldBasedHashGadget<FH, ConstraintF, DataGadget = FpGadget<ConstraintF>>,
    GH: FixedLengthCRH<Output = G>,
    GHG: FixedLengthCRHGadget<GH, ConstraintF, OutputGadget = GG>,
{
    type DataGadget = FpGadget<ConstraintF>;
    type ProofGadget = FieldBasedEcVrfProofGadget<ConstraintF, G, GG>;
    type PublicKeyGadget = FieldBasedEcVrfPkGadget<ConstraintF, G, GG>;
    type GHParametersGadget = GHG::ParametersGadget;

    fn enforce_proof_to_hash_verification<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        group_hash_params: &Self::GHParametersGadget,
        public_key: &Self::PublicKeyGadget,
        proof: &Self::ProofGadget,
        message: Self::DataGadget,
    ) -> Result<Self::DataGadget, SynthesisError> {
        //Check mh = hash_to_curve(message)
        let message_bytes = message.to_bytes_strict(cs.ns(|| "message_to_bytes_restricted"))?;

        let message_on_curve = GHG::check_evaluation_gadget(
            cs.ns(|| "check message_on_curve"),
            group_hash_params,
            message_bytes.as_slice(),
        )?;

        //Serialize c and s

        let c_bits = {
            //Serialize e taking into account the length restriction
            let to_skip = compute_truncation_size(
                ConstraintF::size_in_bits() as i32,
                G::ScalarField::size_in_bits() as i32,
            );

            let c_bits = proof
                .c
                .to_bits_with_length_restriction(cs.ns(|| "c_to_bits"), to_skip)?;

            debug_assert!(c_bits.len() == ConstraintF::size_in_bits() - to_skip);
            c_bits
        };

        let mut s_bits = {
            //Serialize s taking into account the length restriction

            //Before computing the number of bits to truncate from s, we first have to normalize
            //it, i.e. considering its number of bits equals to G::ScalarField::MODULUS_BITS;
            let moduli_diff =
                ConstraintF::size_in_bits() as i32 - G::ScalarField::size_in_bits() as i32;
            let to_skip_init = (if moduli_diff > 0 { moduli_diff } else { 0 }) as usize;

            //Now we can compare the two moduli and decide the bits to truncate
            let to_skip = to_skip_init
                + compute_truncation_size(
                    G::ScalarField::size_in_bits() as i32,
                    ConstraintF::size_in_bits() as i32,
                );

            let s_bits = proof
                .s
                .to_bits_with_length_restriction(cs.ns(|| "s_to_bits"), to_skip as usize)?;

            debug_assert!(s_bits.len() == G::ScalarField::size_in_bits() + to_skip_init - to_skip);
            s_bits
        };
        s_bits.reverse();

        //Hardcode g
        let g = GG::from_value(
            cs.ns(|| "hardcode generator"),
            &G::prime_subgroup_generator(),
        );

        // Random shift to avoid exceptional cases if add is incomplete.
        // With overwhelming probability the circuit will be satisfiable,
        // otherwise the prover can sample another shift by re-running
        // the proof creation.
        let shift = GG::alloc(cs.ns(|| "alloc random shift"), || {
            let mut rng = OsRng::default();
            Ok(loop {
                let r = G::rand(&mut rng);
                if !r.is_zero() {
                    break (r);
                }
            })
        })?;

        //Check u = g^s - pk^c
        let u = {
            let neg_c_times_pk = public_key
                .pk
                .mul_bits(
                    cs.ns(|| "pk * c + shift"),
                    &shift,
                    c_bits.as_slice().iter().rev(),
                )?
                .negate(cs.ns(|| "- (c * pk + shift)"))?;
            GG::mul_bits_fixed_base(
                &g.get_constant(),
                cs.ns(|| "(s * G + shift)"),
                &shift,
                s_bits.as_slice(),
            )?
            // If add is incomplete, and s * G - c * pk = 0, the circuit of the add won't be satisfiable
            .add(cs.ns(|| "(s * G) - (c * pk)"), &neg_c_times_pk)?
        };

        //Check v = mh^s - gamma^c
        let v = {
            let neg_c_times_gamma = proof
                .gamma
                .mul_bits(
                    cs.ns(|| "c * gamma + shift"),
                    &shift,
                    c_bits.as_slice().iter().rev(),
                )?
                .negate(cs.ns(|| "- (c * gamma + shift)"))?;
            message_on_curve
                .mul_bits(
                    cs.ns(|| "(s * mh + shift)"),
                    &shift,
                    s_bits.as_slice().iter(),
                )?
                // If add is incomplete, and s * mh - c * gamma = 0, the circuit of the add won't be satisfiable
                .add(cs.ns(|| "(s * mh) - (c * gamma"), &neg_c_times_gamma)?
        };

        // Check c' = H(m||pk.x||u.x||v.x)
        // Best constraints-efficiency is achieved when m is one field element
        // (or an odd number of field elements).
        let mut hash_input = Vec::new();
        hash_input.push(message.clone());
        hash_input.push(
            public_key
                .pk
                .to_field_gadget_elements(cs.ns(|| "pk to fes"))
                .unwrap()[0]
                .clone(),
        );
        hash_input.push(u.to_field_gadget_elements(cs.ns(|| "u to fes")).unwrap()[0].clone());
        hash_input.push(v.to_field_gadget_elements(cs.ns(|| "v to fes")).unwrap()[0].clone());

        let c_prime =
            FHG::enforce_hash_constant_length(cs.ns(|| "check c_prime"), hash_input.as_slice())?;

        //Enforce c = c'
        proof.c.enforce_equal(cs.ns(|| "check c == c'"), &c_prime)?;

        //Check and return VRF output
        hash_input = Vec::new();
        hash_input.push(message);
        hash_input.extend_from_slice(
            proof
                .gamma
                .to_field_gadget_elements(cs.ns(|| "gamma to fes"))
                .unwrap()
                .as_slice(),
        );

        let vrf_output =
            FHG::enforce_hash_constant_length(cs.ns(|| "check vrf_output"), hash_input.as_slice())?;

        Ok(vrf_output)
    }
}

#[cfg(test)]
mod test {
    use crate::{
        crh::{
            bowe_hopwood::BoweHopwoodPedersenCRHGadget, MNT4PoseidonHashGadget,
            MNT6PoseidonHashGadget,
        },
        vrf::{ecvrf::FieldBasedEcVrfProofVerificationGadget, FieldBasedVrfGadget},
    };
    use algebra::curves::{
        mnt4753::G1Projective as MNT4G1Projective, mnt6753::G1Projective as MNT6G1Projective,
    };
    use algebra::fields::{mnt4753::Fr as MNT4Fr, mnt6753::Fr as MNT6Fr};
    use primitives::{
        crh::{
            bowe_hopwood::{BoweHopwoodPedersenCRH, BoweHopwoodPedersenParameters},
            pedersen::PedersenWindow,
            FixedLengthCRH, MNT4PoseidonHash, MNT6PoseidonHash,
        },
        vrf::{
            ecvrf::{FieldBasedEcVrf, FieldBasedEcVrfProof},
            FieldBasedVrf,
        },
    };

    use r1cs_core::ConstraintSystem;
    use r1cs_std::alloc::AllocGadget;
    use r1cs_std::instantiated::{
        mnt4_753::G1Gadget as MNT4G1Gadget, mnt6_753::G1Gadget as MNT6G1Gadget,
    };
    use r1cs_std::test_constraint_system::TestConstraintSystem;

    use primitives::vrf::ecvrf::FieldBasedEcVrfPk;
    use rand::{thread_rng, Rng};

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

    type BHMNT4Parameters = BoweHopwoodPedersenParameters<MNT6G1Projective>;
    type BHMNT6Parameters = BoweHopwoodPedersenParameters<MNT4G1Projective>;

    type EcVrfMNT4 = FieldBasedEcVrf<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash, BHMNT6>;
    type EcVrfMNT6 = FieldBasedEcVrf<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash, BHMNT4>;

    type EcVrfMNT4Pk = FieldBasedEcVrfPk<MNT6G1Projective>;
    type EcVrfMNT6Pk = FieldBasedEcVrfPk<MNT4G1Projective>;

    type EcVrfMNT4Proof = FieldBasedEcVrfProof<MNT4Fr, MNT6G1Projective>;
    type EcVrfMNT6Proof = FieldBasedEcVrfProof<MNT6Fr, MNT4G1Projective>;

    type EcVrfMNT4Gadget = FieldBasedEcVrfProofVerificationGadget<
        MNT4Fr,
        MNT6G1Projective,
        MNT6G1Gadget,
        MNT4PoseidonHash,
        MNT4PoseidonHashGadget,
        BHMNT6,
        BHMNT4Gadget,
    >;

    type EcVrfMNT6Gadget = FieldBasedEcVrfProofVerificationGadget<
        MNT6Fr,
        MNT4G1Projective,
        MNT4G1Gadget,
        MNT6PoseidonHash,
        MNT6PoseidonHashGadget,
        BHMNT4,
        BHMNT6Gadget,
    >;

    fn prove<S: FieldBasedVrf, R: Rng>(
        rng: &mut R,
        pp: &S::GHParams,
        message: S::Data,
    ) -> (S::Proof, S::PublicKey) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));
        let proof = S::prove(rng, pp, &pk, &sk, message).unwrap();
        (proof, pk)
    }

    fn mnt4_ecvrf_gadget_generate_constraints(
        message: MNT4Fr,
        pk: &EcVrfMNT4Pk,
        proof: EcVrfMNT4Proof,
        pp: &BHMNT4Parameters,
    ) -> bool {
        let mut cs = TestConstraintSystem::<MNT4Fr>::new();

        //Alloc proof, pk and message
        let proof_g =
            <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::ProofGadget::alloc(
                cs.ns(|| "alloc proof"),
                || Ok(proof),
            )
            .unwrap();

        let pk_g =
            <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::PublicKeyGadget::alloc(
                cs.ns(|| "alloc pk"),
                || Ok(pk),
            )
            .unwrap();

        let pp_g =
            <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::GHParametersGadget::alloc(
                cs.ns(|| "alloc gh params"),
                || Ok(pp),
            )
            .unwrap();

        let message_g =
            <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::DataGadget::alloc(
                cs.ns(|| "alloc message"),
                || Ok(message),
            )
            .unwrap();

        //Verify proof
        EcVrfMNT4Gadget::enforce_proof_to_hash_verification(
            cs.ns(|| "verify proof1"),
            &pp_g,
            &pk_g,
            &proof_g,
            message_g,
        )
        .unwrap();

        if !cs.is_satisfied() {
            println!("**********Unsatisfied constraints***********");
            println!("{:?}", cs.which_is_unsatisfied());
        }

        cs.is_satisfied()
    }

    #[test]
    fn mnt4_ecvrf_gadget_test() {
        let rng = &mut thread_rng();
        let message: MNT4Fr = rng.gen();
        let pp = <BHMNT6 as FixedLengthCRH>::setup(rng).unwrap();
        let (proof, pk) = prove::<EcVrfMNT4, _>(rng, &pp, message);

        //Positive case
        assert!(mnt4_ecvrf_gadget_generate_constraints(
            message, &pk, proof, &pp
        ));

        //Change message
        let wrong_message: MNT4Fr = rng.gen();
        assert!(!mnt4_ecvrf_gadget_generate_constraints(
            wrong_message,
            &pk,
            proof,
            &pp
        ));

        //Change pk
        let wrong_pk: EcVrfMNT4Pk = rng.gen();
        assert!(!mnt4_ecvrf_gadget_generate_constraints(
            message, &wrong_pk, proof, &pp
        ));

        //Change proof
        let (wrong_proof, _) = prove::<EcVrfMNT4, _>(rng, &pp, wrong_message);
        assert!(!mnt4_ecvrf_gadget_generate_constraints(
            message,
            &pk,
            wrong_proof,
            &pp
        ));
    }

    fn mnt6_ecvrf_gadget_generate_constraints(
        message: MNT6Fr,
        pk: &EcVrfMNT6Pk,
        proof: EcVrfMNT6Proof,
        pp: &BHMNT6Parameters,
    ) -> bool {
        let mut cs = TestConstraintSystem::<MNT6Fr>::new();

        //Alloc proof, pk and message
        let proof_g =
            <EcVrfMNT6Gadget as FieldBasedVrfGadget<EcVrfMNT6, MNT6Fr>>::ProofGadget::alloc(
                cs.ns(|| "alloc proof"),
                || Ok(proof),
            )
            .unwrap();
        let pk_g =
            <EcVrfMNT6Gadget as FieldBasedVrfGadget<EcVrfMNT6, MNT6Fr>>::PublicKeyGadget::alloc(
                cs.ns(|| "alloc pk"),
                || Ok(pk),
            )
            .unwrap();
        let pp_g =
            <EcVrfMNT6Gadget as FieldBasedVrfGadget<EcVrfMNT6, MNT6Fr>>::GHParametersGadget::alloc(
                cs.ns(|| "alloc gh params"),
                || Ok(pp),
            )
            .unwrap();
        let message_g =
            <EcVrfMNT6Gadget as FieldBasedVrfGadget<EcVrfMNT6, MNT6Fr>>::DataGadget::alloc(
                cs.ns(|| "alloc message"),
                || Ok(message),
            )
            .unwrap();

        //Verify proof
        EcVrfMNT6Gadget::enforce_proof_to_hash_verification(
            cs.ns(|| "verify proof1"),
            &pp_g,
            &pk_g,
            &proof_g,
            message_g,
        )
        .unwrap();

        if !cs.is_satisfied() {
            println!("**********Unsatisfied constraints***********");
            println!("{:?}", cs.which_is_unsatisfied());
        }

        cs.is_satisfied()
    }

    #[ignore]
    #[test]
    fn mnt6_ecvrf_gadget_test() {
        let rng = &mut thread_rng();
        let message: MNT6Fr = rng.gen();
        let pp = <BHMNT4 as FixedLengthCRH>::setup(rng).unwrap();
        let (proof, pk) = prove::<EcVrfMNT6, _>(rng, &pp, message);

        //Positive case
        assert!(mnt6_ecvrf_gadget_generate_constraints(
            message, &pk, proof, &pp
        ));

        //Change message
        let wrong_message: MNT6Fr = rng.gen();
        assert!(!mnt6_ecvrf_gadget_generate_constraints(
            wrong_message,
            &pk,
            proof,
            &pp
        ));

        //Change pk
        let wrong_pk: EcVrfMNT6Pk = rng.gen();
        assert!(!mnt6_ecvrf_gadget_generate_constraints(
            message, &wrong_pk, proof, &pp
        ));

        //Change proof
        let (wrong_proof, _) = prove::<EcVrfMNT6, _>(rng, &pp, wrong_message);
        assert!(!mnt6_ecvrf_gadget_generate_constraints(
            message,
            &pk,
            wrong_proof,
            &pp
        ));
    }

    #[ignore]
    #[test]
    fn random_ecvrf_gadget_test() {
        //Generate VRF proof for a random field element f and get the proof and the public key
        let rng = &mut thread_rng();
        let pp = <BHMNT6 as FixedLengthCRH>::setup(rng).unwrap();

        let samples = 10;
        for _ in 0..samples {
            let message: MNT4Fr = rng.gen();
            let (sig, pk) = prove::<EcVrfMNT4, _>(rng, &pp, message);
            let mut cs = TestConstraintSystem::<MNT4Fr>::new();

            //Alloc proof, pk, hash params and message
            let proof_g =
                <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::ProofGadget::alloc(
                    cs.ns(|| "alloc proof"),
                    || Ok(sig),
                )
                .unwrap();

            let pk_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::PublicKeyGadget::alloc(
                cs.ns(|| "alloc pk"),
                || Ok(pk)
            ).unwrap();

            let pp_g = <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::GHParametersGadget::alloc(
                cs.ns(|| "alloc gh params"),
                || Ok(&pp)
            ).unwrap();

            let message_g =
                <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::DataGadget::alloc(
                    cs.ns(|| "alloc message"),
                    || Ok(message),
                )
                .unwrap();

            //Verify proof
            EcVrfMNT4Gadget::enforce_proof_to_hash_verification(
                cs.ns(|| "verify proof"),
                &pp_g,
                &pk_g,
                &proof_g,
                message_g,
            )
            .unwrap();

            if !cs.is_satisfied() {
                println!("**********Unsatisfied constraints***********");
                println!("{:?}", cs.which_is_unsatisfied());
            }

            assert!(cs.is_satisfied());

            //Negative case: wrong message (or wrong proof for another message)
            let new_message: MNT4Fr = rng.gen();

            let new_message_g =
                <EcVrfMNT4Gadget as FieldBasedVrfGadget<EcVrfMNT4, MNT4Fr>>::DataGadget::alloc(
                    cs.ns(|| "alloc new_message"),
                    || Ok(new_message),
                )
                .unwrap();

            EcVrfMNT4Gadget::enforce_proof_to_hash_verification(
                cs.ns(|| "verify new proof"),
                &pp_g,
                &pk_g,
                &proof_g,
                new_message_g,
            )
            .unwrap();

            assert!(!cs.is_satisfied());
        }
    }
}
