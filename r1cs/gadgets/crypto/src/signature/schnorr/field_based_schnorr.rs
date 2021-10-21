use crate::{crh::FieldBasedHashGadget, signature::FieldBasedSigGadget};
use algebra::{Group, PrimeField, ProjectiveCurve, ToConstraintField};
use primitives::signature::schnorr::field_based_schnorr::FieldBasedSchnorrPk;
use primitives::{
    compute_truncation_size,
    crh::FieldBasedHash,
    signature::schnorr::field_based_schnorr::{
        FieldBasedSchnorrSignature, FieldBasedSchnorrSignatureScheme,
    },
};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::alloc::ConstantGadget;
use r1cs_std::{
    alloc::AllocGadget, bits::boolean::Boolean, eq::EqGadget, fields::fp::FpGadget,
    groups::GroupGadget, to_field_gadget_vec::ToConstraintFieldGadget,
};
use rand::rngs::OsRng;
use std::{borrow::Borrow, marker::PhantomData};

#[derive(Derivative)]
#[derivative(
    Debug(bound = "ConstraintF: PrimeField, G: Group"),
    Clone(bound = "ConstraintF: PrimeField, G: Group"),
    PartialEq(bound = "ConstraintF: PrimeField, G: Group"),
    Eq(bound = "ConstraintF: PrimeField, G: Group")
)]
pub struct FieldBasedSchnorrSigGadget<ConstraintF: PrimeField, G: Group> {
    pub e: FpGadget<ConstraintF>,
    pub s: FpGadget<ConstraintF>,
    _field: PhantomData<ConstraintF>,
    _group: PhantomData<G>,
}

impl<ConstraintF, G> AllocGadget<FieldBasedSchnorrSignature<ConstraintF, G>, ConstraintF>
    for FieldBasedSchnorrSigGadget<ConstraintF, G>
where
    ConstraintF: PrimeField,
    G: Group,
{
    fn alloc<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedSchnorrSignature<ConstraintF, G>>,
    {
        let (e, s) = match f() {
            Ok(sig) => {
                let sig = *sig.borrow();
                (Ok(sig.e), Ok(sig.s))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let e = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc e"), || e)?;
        let s = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc s"), || s)?;
        Ok(Self {
            e,
            s,
            _field: PhantomData,
            _group: PhantomData,
        })
    }

    fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        f: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<FieldBasedSchnorrSignature<ConstraintF, G>>,
    {
        let (e, s) = match f() {
            Ok(sig) => {
                let sig = *sig.borrow();
                (Ok(sig.e), Ok(sig.s))
            }
            _ => (
                Err(SynthesisError::AssignmentMissing),
                Err(SynthesisError::AssignmentMissing),
            ),
        };

        let e = FpGadget::<ConstraintF>::alloc_input(cs.ns(|| "alloc e"), || e)?;
        let s = FpGadget::<ConstraintF>::alloc_input(cs.ns(|| "alloc s"), || s)?;
        Ok(Self {
            e,
            s,
            _field: PhantomData,
            _group: PhantomData,
        })
    }
}

impl<ConstraintF, G> ConstantGadget<FieldBasedSchnorrSignature<ConstraintF, G>, ConstraintF>
    for FieldBasedSchnorrSigGadget<ConstraintF, G>
where
    ConstraintF: PrimeField,
    G: Group,
{
    fn from_value<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value: &FieldBasedSchnorrSignature<ConstraintF, G>,
    ) -> Self {
        let e = FpGadget::<ConstraintF>::from_value(cs.ns(|| "hardcode e"), &value.e);
        let s = FpGadget::<ConstraintF>::from_value(cs.ns(|| "hardcode s"), &value.s);
        Self {
            e,
            s,
            _field: PhantomData,
            _group: PhantomData,
        }
    }

    fn get_constant(&self) -> FieldBasedSchnorrSignature<ConstraintF, G> {
        let e = self.e.value.unwrap();
        let s = self.s.value.unwrap();
        FieldBasedSchnorrSignature::<ConstraintF, G>::new(e, s)
    }
}

impl<ConstraintF, G> EqGadget<ConstraintF> for FieldBasedSchnorrSigGadget<ConstraintF, G>
where
    ConstraintF: PrimeField,
    G: Group,
{
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        let b1 = self.e.is_eq(cs.ns(|| "b1"), &other.e)?;
        let b2 = self.s.is_eq(cs.ns(|| "b2"), &other.s)?;
        Boolean::and(cs.ns(|| "b1 && b2"), &b1, &b2)
    }

    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.e.conditional_enforce_equal(
            cs.ns(|| "self.e =? other.e"),
            &other.e,
            should_enforce,
        )?;
        self.s.conditional_enforce_equal(
            cs.ns(|| "self.s =? other.s"),
            &other.s,
            should_enforce,
        )?;
        Ok(())
    }

    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.e.conditional_enforce_not_equal(
            cs.ns(|| "self.e !=? other.e"),
            &other.e,
            should_enforce,
        )?;
        self.s.conditional_enforce_not_equal(
            cs.ns(|| "self.s !=? other.s"),
            &other.s,
            should_enforce,
        )?;
        Ok(())
    }
}

impl<ConstraintF, G> ToConstraintFieldGadget<ConstraintF>
    for FieldBasedSchnorrSigGadget<ConstraintF, G>
where
    ConstraintF: PrimeField,
    G: Group,
{
    type FieldGadget = FpGadget<ConstraintF>;

    fn to_field_gadget_elements<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Self::FieldGadget>, SynthesisError> {
        Ok(vec![self.e.clone(), self.s.clone()])
    }
}

#[derive(Clone, Eq, PartialEq)]
pub struct FieldBasedSchnorrPkGadget<
    ConstraintF: PrimeField,
    G: Group,
    GG: GroupGadget<G, ConstraintF>,
> {
    pub pk: GG,
    _field: PhantomData<ConstraintF>,
    _group: PhantomData<G>,
}

impl<ConstraintF, G, GG> AllocGadget<FieldBasedSchnorrPk<G>, ConstraintF>
    for FieldBasedSchnorrPkGadget<ConstraintF, G, GG>
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
        T: Borrow<FieldBasedSchnorrPk<G>>,
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
        T: Borrow<FieldBasedSchnorrPk<G>>,
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
        T: Borrow<FieldBasedSchnorrPk<G>>,
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
        T: Borrow<FieldBasedSchnorrPk<G>>,
    {
        let pk = GG::alloc_input(cs.ns(|| "alloc pk"), || f().map(|pk| pk.borrow().0))?;
        Ok(Self {
            pk,
            _field: PhantomData,
            _group: PhantomData,
        })
    }
}

impl<ConstraintF, G, GG> ConstantGadget<FieldBasedSchnorrPk<G>, ConstraintF>
    for FieldBasedSchnorrPkGadget<ConstraintF, G, GG>
where
    ConstraintF: PrimeField,
    G: Group,
    GG: GroupGadget<G, ConstraintF, Value = G>,
{
    fn from_value<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value: &FieldBasedSchnorrPk<G>,
    ) -> Self {
        let pk = GG::from_value(cs.ns(|| "hardcode pk"), &value.0);
        Self {
            pk,
            _field: PhantomData,
            _group: PhantomData,
        }
    }

    fn get_constant(&self) -> FieldBasedSchnorrPk<G> {
        FieldBasedSchnorrPk::<G>(self.pk.get_value().unwrap())
    }
}

impl<ConstraintF, G, GG> EqGadget<ConstraintF> for FieldBasedSchnorrPkGadget<ConstraintF, G, GG>
where
    ConstraintF: PrimeField,
    G: Group,
    GG: GroupGadget<G, ConstraintF, Value = G>,
{
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        self.pk.is_eq(cs, &other.pk)
    }

    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.pk
            .conditional_enforce_equal(cs, &other.pk, should_enforce)
    }

    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.pk
            .conditional_enforce_not_equal(cs, &other.pk, should_enforce)
    }
}

impl<ConstraintF, G, GG> ToConstraintFieldGadget<ConstraintF>
    for FieldBasedSchnorrPkGadget<ConstraintF, G, GG>
where
    ConstraintF: PrimeField,
    G: Group,
    GG: GroupGadget<G, ConstraintF, Value = G>
        + ToConstraintFieldGadget<ConstraintF, FieldGadget = FpGadget<ConstraintF>>,
{
    type FieldGadget = FpGadget<ConstraintF>;

    fn to_field_gadget_elements<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Self::FieldGadget>, SynthesisError> {
        self.pk.to_field_gadget_elements(cs)
    }
}

pub struct FieldBasedSchnorrSigVerificationGadget<
    ConstraintF: PrimeField,
    G: Group,
    GG: GroupGadget<G, ConstraintF>,
    H: FieldBasedHash<Data = ConstraintF>,
    HG: FieldBasedHashGadget<H, ConstraintF>,
> {
    _field: PhantomData<ConstraintF>,
    _group: PhantomData<G>,
    _group_gadget: PhantomData<GG>,
    _hash: PhantomData<H>,
    _hash_gadget: PhantomData<HG>,
}

// This implementation supports both complete and incomplete (safe) point addition.
// Assumes provided key material to be already checked.
//
// In case of incomplete point addition, with negligible probability, the
// proof creation might fail at first attempt and must be re-run (in order to sample
// fresh randomnesses).
// Furthermore, one exceptional case (e, s) has to be treated outside the circuit:
// if e * pk = s * G, i.e. when R' is trivial (therefore leaking the sk), then
// the circuit is not satisfiable.
impl<ConstraintF, G, GG, H, HG> FieldBasedSchnorrSigVerificationGadget<ConstraintF, G, GG, H, HG>
where
    ConstraintF: PrimeField,
    G: ProjectiveCurve + ToConstraintField<ConstraintF>,
    GG: GroupGadget<G, ConstraintF, Value = G>
        + ToConstraintFieldGadget<ConstraintF, FieldGadget = HG::DataGadget>,
    H: FieldBasedHash<Data = ConstraintF>,
    HG: FieldBasedHashGadget<H, ConstraintF, DataGadget = FpGadget<ConstraintF>>,
{
    fn enforce_signature_computation<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        public_key: &GG,
        signature: &FieldBasedSchnorrSigGadget<ConstraintF, G>,
        message: FpGadget<ConstraintF>,
    ) -> Result<FpGadget<ConstraintF>, SynthesisError> {
        //Enforce e' * pk
        let e_bits = {
            //Serialize e taking into account the length restriction
            let to_skip = compute_truncation_size(
                ConstraintF::size_in_bits() as i32,
                G::ScalarField::size_in_bits() as i32,
            );

            let e_bits = signature
                .e
                .to_bits_with_length_restriction(cs.ns(|| "e_to_bits"), to_skip)?;

            debug_assert!(e_bits.len() == ConstraintF::size_in_bits() - to_skip);
            e_bits
        };

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

        let neg_e_times_pk = public_key
            .mul_bits(
                cs.ns(|| "pk * e + shift"),
                &shift,
                e_bits.as_slice().iter().rev(),
            )?
            .negate(cs.ns(|| "- (pk * e + shift)"))?;

        //Enforce s * G and R' = s*G - e*pk
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

            let s_bits = signature
                .s
                .to_bits_with_length_restriction(cs.ns(|| "s_to_bits"), to_skip as usize)?;

            debug_assert!(s_bits.len() == G::ScalarField::size_in_bits() + to_skip_init - to_skip);
            s_bits
        };

        s_bits.reverse();
        let g = GG::from_value(
            cs.ns(|| "hardcode generator"),
            &G::prime_subgroup_generator(),
        );
        let r_prime = GG::mul_bits_fixed_base(
            &g.get_constant(),
            cs.ns(|| "(s * G + shift)"),
            &shift,
            s_bits.as_slice(),
        )?
        // If add is incomplete, and s * G - e * pk = 0, the circuit of the add won't be satisfiable
        .add(cs.ns(|| "s * G - e * pk "), &neg_e_times_pk)?;

        let r_prime_coords = r_prime.to_field_gadget_elements(cs.ns(|| "r_prime to fes"))?;

        // Check e' = H(m || R' || pk.x)
        // Best constraints-efficiency is achieved when m is one field element
        // (or an odd number of field elements).
        let mut hash_input = Vec::new();
        hash_input.push(message);
        hash_input.extend_from_slice(r_prime_coords.as_slice());
        hash_input.push(
            public_key
                .to_field_gadget_elements(cs.ns(|| "pk to fes"))
                .unwrap()[0]
                .clone(),
        );

        HG::enforce_hash_constant_length(cs.ns(|| "check e_prime"), hash_input.as_slice())
    }
}

impl<ConstraintF, G, GG, H, HG>
    FieldBasedSigGadget<FieldBasedSchnorrSignatureScheme<ConstraintF, G, H>, ConstraintF>
    for FieldBasedSchnorrSigVerificationGadget<ConstraintF, G, GG, H, HG>
where
    ConstraintF: PrimeField,
    G: ProjectiveCurve + ToConstraintField<ConstraintF>,
    GG: GroupGadget<G, ConstraintF, Value = G>
        + ToConstraintFieldGadget<ConstraintF, FieldGadget = HG::DataGadget>,
    H: FieldBasedHash<Data = ConstraintF>,
    HG: FieldBasedHashGadget<H, ConstraintF, DataGadget = FpGadget<ConstraintF>>,
{
    type DataGadget = FpGadget<ConstraintF>;
    type SignatureGadget = FieldBasedSchnorrSigGadget<ConstraintF, G>;
    type PublicKeyGadget = FieldBasedSchnorrPkGadget<ConstraintF, G, GG>;

    fn enforce_signature_verdict<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        public_key: &Self::PublicKeyGadget,
        signature: &Self::SignatureGadget,
        message: Self::DataGadget,
    ) -> Result<Boolean, SynthesisError> {
        let e_prime = Self::enforce_signature_computation(
            cs.ns(|| "is sig verified"),
            &public_key.pk,
            signature,
            message,
        )?;

        //Enforce result of signature verification
        let is_verified = signature.e.is_eq(cs.ns(|| "is e == e_prime"), &e_prime)?;

        Ok(is_verified)
    }

    fn conditionally_enforce_signature_verification<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        public_key: &Self::PublicKeyGadget,
        signature: &Self::SignatureGadget,
        message: Self::DataGadget,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        let e_prime = Self::enforce_signature_computation(
            cs.ns(|| "is sig verified"),
            &public_key.pk,
            signature,
            message,
        )?;
        signature.e.conditional_enforce_equal(
            cs.ns(|| "conditional verify signature"),
            &e_prime,
            should_enforce,
        )?;
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use algebra::curves::{
        mnt4753::G1Projective as MNT4G1Projective, mnt6753::G1Projective as MNT6G1Projective,
    };
    use algebra::fields::{mnt4753::Fr as MNT4Fr, mnt6753::Fr as MNT6Fr};

    use primitives::{
        crh::{MNT4PoseidonHash, MNT6PoseidonHash},
        signature::{schnorr::field_based_schnorr::*, FieldBasedSignatureScheme},
    };

    use crate::{
        crh::{MNT4PoseidonHashGadget, MNT6PoseidonHashGadget},
        signature::{schnorr::field_based_schnorr::*, FieldBasedSigGadget},
    };

    use r1cs_core::ConstraintSystem;
    use r1cs_std::alloc::AllocGadget;

    use r1cs_std::instantiated::{
        mnt4_753::G1Gadget as MNT4G1Gadget, mnt6_753::G1Gadget as MNT6G1Gadget,
    };

    use r1cs_std::test_constraint_system::TestConstraintSystem;
    use rand::{thread_rng, Rng};

    type SchnorrMNT4 = FieldBasedSchnorrSignatureScheme<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash>;
    type SchnorrMNT6 = FieldBasedSchnorrSignatureScheme<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash>;

    type SchnorrMNT4Sig = FieldBasedSchnorrSignature<MNT4Fr, MNT6G1Projective>;
    type SchnorrMNT6Sig = FieldBasedSchnorrSignature<MNT6Fr, MNT4G1Projective>;

    type SchnorrMNT4Pk = FieldBasedSchnorrPk<MNT6G1Projective>;
    type SchnorrMNT6Pk = FieldBasedSchnorrPk<MNT4G1Projective>;

    type SchnorrMNT4Gadget = FieldBasedSchnorrSigVerificationGadget<
        MNT4Fr,
        MNT6G1Projective,
        MNT6G1Gadget,
        MNT4PoseidonHash,
        MNT4PoseidonHashGadget,
    >;
    type SchnorrMNT6Gadget = FieldBasedSchnorrSigVerificationGadget<
        MNT6Fr,
        MNT4G1Projective,
        MNT4G1Gadget,
        MNT6PoseidonHash,
        MNT6PoseidonHashGadget,
    >;

    fn sign<S: FieldBasedSignatureScheme, R: Rng>(
        rng: &mut R,
        message: S::Data,
    ) -> (S::Signature, S::PublicKey) {
        let (pk, sk) = S::keygen(rng);
        assert!(S::keyverify(&pk));
        let sig = S::sign(rng, &pk, &sk, message).unwrap();
        (sig, pk)
    }

    fn mnt4_schnorr_gadget_generate_constraints(
        message: MNT4Fr,
        pk: &SchnorrMNT4Pk,
        sig: SchnorrMNT4Sig,
    ) -> bool {
        let mut cs = TestConstraintSystem::<MNT4Fr>::new();

        //Alloc signature, pk and message
        let sig_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc sig"),
            || Ok(sig)
        ).unwrap();
        let pk_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::PublicKeyGadget::alloc(cs.ns(|| "alloc pk"), || Ok(pk)).unwrap();
        let message_g =
            <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::DataGadget::alloc(
                cs.ns(|| "alloc message"),
                || Ok(message),
            )
            .unwrap();

        //Verify sig
        SchnorrMNT4Gadget::enforce_signature_verification(
            cs.ns(|| "verify sig1"),
            &pk_g,
            &sig_g,
            message_g.clone(),
        )
        .unwrap();

        let is_cs_satisfied = cs.is_satisfied();

        //Verify sig
        let is_verified = SchnorrMNT4Gadget::enforce_signature_verdict(
            cs.ns(|| "sig1 result"),
            &pk_g,
            &sig_g,
            message_g,
        )
        .unwrap();

        assert_eq!(is_verified.get_value().unwrap(), is_cs_satisfied);

        if !is_cs_satisfied {
            println!("**********Unsatisfied constraints***********");
            println!("{:?}", cs.which_is_unsatisfied());
        }

        is_cs_satisfied
    }

    #[test]
    fn mnt4_schnorr_gadget_test() {
        //Sign a random field element f and get the signature and the public key
        let rng = &mut thread_rng();
        let message: MNT4Fr = rng.gen();
        let (sig, pk) = sign::<SchnorrMNT4, _>(rng, message);

        //Positive case
        assert!(mnt4_schnorr_gadget_generate_constraints(message, &pk, sig));

        //Change message
        let wrong_message: MNT4Fr = rng.gen();
        assert!(!mnt4_schnorr_gadget_generate_constraints(
            wrong_message,
            &pk,
            sig
        ));

        //Change pk
        let wrong_pk: SchnorrMNT4Pk = rng.gen();
        assert!(!mnt4_schnorr_gadget_generate_constraints(
            message, &wrong_pk, sig
        ));

        //Change sig
        let (wrong_sig, _) = sign::<SchnorrMNT4, _>(rng, wrong_message);
        assert!(!mnt4_schnorr_gadget_generate_constraints(
            message, &pk, wrong_sig
        ));
    }

    fn mnt6_schnorr_gadget_generate_constraints(
        message: MNT6Fr,
        pk: &SchnorrMNT6Pk,
        sig: SchnorrMNT6Sig,
    ) -> bool {
        let mut cs = TestConstraintSystem::<MNT6Fr>::new();

        //Alloc signature, pk and message
        let sig_g = <SchnorrMNT6Gadget as FieldBasedSigGadget<SchnorrMNT6, MNT6Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc sig"),
            || Ok(sig)
        ).unwrap();
        let pk_g = <SchnorrMNT6Gadget as FieldBasedSigGadget<SchnorrMNT6, MNT6Fr>>::PublicKeyGadget::alloc(cs.ns(|| "alloc pk"), || Ok(pk)).unwrap();
        let message_g =
            <SchnorrMNT6Gadget as FieldBasedSigGadget<SchnorrMNT6, MNT6Fr>>::DataGadget::alloc(
                cs.ns(|| "alloc message"),
                || Ok(message),
            )
            .unwrap();

        //Verify sig
        SchnorrMNT6Gadget::enforce_signature_verification(
            cs.ns(|| "verify sig1"),
            &pk_g,
            &sig_g,
            message_g.clone(),
        )
        .unwrap();

        let is_cs_satisfied = cs.is_satisfied();

        let is_verified = SchnorrMNT6Gadget::enforce_signature_verdict(
            cs.ns(|| "sig1 result"),
            &pk_g,
            &sig_g,
            message_g,
        )
        .unwrap();

        assert_eq!(is_verified.get_value().unwrap(), is_cs_satisfied);

        if !is_cs_satisfied {
            println!("**********Unsatisfied constraints***********");
            println!("{:?}", cs.which_is_unsatisfied());
        }

        is_cs_satisfied
    }

    #[ignore]
    #[test]
    fn mnt6_schnorr_gadget_test() {
        //Sign a random field element f and get the signature and the public key
        let rng = &mut thread_rng();
        let message: MNT6Fr = rng.gen();
        let (sig, pk) = sign::<SchnorrMNT6, _>(rng, message);

        //Positive case
        assert!(mnt6_schnorr_gadget_generate_constraints(message, &pk, sig));

        //Change message
        let wrong_message: MNT6Fr = rng.gen();
        assert!(!mnt6_schnorr_gadget_generate_constraints(
            wrong_message,
            &pk,
            sig
        ));

        //Change pk
        let wrong_pk: SchnorrMNT6Pk = rng.gen();
        assert!(!mnt6_schnorr_gadget_generate_constraints(
            message, &wrong_pk, sig
        ));

        //Change sig
        let (wrong_sig, _) = sign::<SchnorrMNT6, _>(rng, wrong_message);
        assert!(!mnt6_schnorr_gadget_generate_constraints(
            message, &pk, wrong_sig
        ));
    }

    #[ignore]
    #[test]
    fn random_schnorr_gadget_test() {
        let rng = &mut thread_rng();

        let samples = 10;
        for _ in 0..samples {
            let message: MNT4Fr = rng.gen();
            let (sig, pk) = sign::<SchnorrMNT4, _>(rng, message);
            let mut cs = TestConstraintSystem::<MNT4Fr>::new();

            //Alloc signature, pk and message
            let sig_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::SignatureGadget::alloc(
                cs.ns(|| "alloc sig"),
                || Ok(sig)
            ).unwrap();

            let pk_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::PublicKeyGadget::alloc(
                cs.ns(|| "alloc pk"),
                || Ok(pk)
            ).unwrap();

            let message_g =
                <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::DataGadget::alloc(
                    cs.ns(|| "alloc message"),
                    || Ok(message),
                )
                .unwrap();

            //Verify sig
            let is_verified = SchnorrMNT4Gadget::enforce_signature_verdict(
                cs.ns(|| "sig result"),
                &pk_g,
                &sig_g,
                message_g.clone(),
            )
            .unwrap();

            assert!(is_verified.get_value().unwrap());

            SchnorrMNT4Gadget::enforce_signature_verification(
                cs.ns(|| "verify sig"),
                &pk_g,
                &sig_g,
                message_g,
            )
            .unwrap();

            assert!(cs.is_satisfied());

            //Negative case: wrong message (or wrong sig for another message)
            let new_message: MNT4Fr = rng.gen();
            let new_message_g =
                <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::DataGadget::alloc(
                    cs.ns(|| "alloc new_message"),
                    || Ok(new_message),
                )
                .unwrap();

            let is_verified = SchnorrMNT4Gadget::enforce_signature_verdict(
                cs.ns(|| "new sig result"),
                &pk_g,
                &sig_g,
                new_message_g.clone(),
            )
            .unwrap();

            if !cs.is_satisfied() {
                println!("**********Unsatisfied constraints***********");
                println!("{:?}", cs.which_is_unsatisfied());
            }

            assert!(!is_verified.get_value().unwrap());
            assert!(cs.is_satisfied());

            SchnorrMNT4Gadget::enforce_signature_verification(
                cs.ns(|| "verify new sig"),
                &pk_g,
                &sig_g,
                new_message_g,
            )
            .unwrap();

            assert!(!cs.is_satisfied());
        }
    }
}
