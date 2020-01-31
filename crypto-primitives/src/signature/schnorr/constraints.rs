use algebra::{groups::Group, Field};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use crate::signature::SigRandomizePkGadget;

use std::{borrow::Borrow, marker::PhantomData};

use crate::signature::schnorr::{SchnorrPublicKey, SchnorrSigParameters, SchnorrSignature};
use digest::Digest;

pub struct SchnorrSigGadgetParameters<G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>>
{
    generator: GG,
    _group:    PhantomData<*const G>,
    _engine:   PhantomData<*const ConstraintF>,
}

impl<G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>> Clone
    for SchnorrSigGadgetParameters<G, ConstraintF, GG>
{
    fn clone(&self) -> Self {
        Self {
            generator: self.generator.clone(),
            _group:    PhantomData,
            _engine:   PhantomData,
        }
    }
}

#[derive(Derivative)]
#[derivative(
    Debug(bound = "G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>"),
    Clone(bound = "G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>"),
    PartialEq(bound = "G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>"),
    Eq(bound = "G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>")
)]
pub struct SchnorrSigGadgetPk<G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>> {
    pub_key: GG,
    #[doc(hidden)]
    _group: PhantomData<*const G>,
    #[doc(hidden)]
    _engine: PhantomData<*const ConstraintF>,
}

pub struct SchnorrRandomizePkGadget<G: Group, ConstraintF: Field, GG: GroupGadget<G, ConstraintF>> {
    #[doc(hidden)]
    _group: PhantomData<*const G>,
    #[doc(hidden)]
    _group_gadget: PhantomData<*const GG>,
    #[doc(hidden)]
    _engine: PhantomData<*const ConstraintF>,
}

impl<G, GG, D, ConstraintF> SigRandomizePkGadget<SchnorrSignature<G, D>, ConstraintF>
    for SchnorrRandomizePkGadget<G, ConstraintF, GG>
where
    G: Group,
    GG: GroupGadget<G, ConstraintF>,
    D: Digest + Send + Sync,
    ConstraintF: Field,
{
    type ParametersGadget = SchnorrSigGadgetParameters<G, ConstraintF, GG>;
    type PublicKeyGadget = SchnorrSigGadgetPk<G, ConstraintF, GG>;

    fn check_randomization_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        parameters: &Self::ParametersGadget,
        public_key: &Self::PublicKeyGadget,
        randomness: &[UInt8],
    ) -> Result<Self::PublicKeyGadget, SynthesisError> {
        let base = parameters.generator.clone();
        let randomness = randomness
            .iter()
            .flat_map(|b| b.into_bits_le())
            .collect::<Vec<_>>();
        let rand_pk = base.mul_bits(
            &mut cs.ns(|| "Compute Randomizer"),
            &public_key.pub_key,
            randomness.iter(),
        )?;
        Ok(SchnorrSigGadgetPk {
            pub_key: rand_pk,
            _group:  PhantomData,
            _engine: PhantomData,
        })
    }
}

impl<G, ConstraintF, GG, D> AllocGadget<SchnorrSigParameters<G, D>, ConstraintF>
    for SchnorrSigGadgetParameters<G, ConstraintF, GG>
where
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
    D: Digest,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SchnorrSigParameters<G, D>>,
    {
        let generator = GG::alloc_checked(cs, || f().map(|pp| pp.borrow().generator))?;
        Ok(Self {
            generator,
            _engine: PhantomData,
            _group: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SchnorrSigParameters<G, D>>,
    {
        let generator = GG::alloc_input(cs, || f().map(|pp| pp.borrow().generator))?;
        Ok(Self {
            generator,
            _engine: PhantomData,
            _group: PhantomData,
        })
    }
}

impl<G, ConstraintF, GG> AllocGadget<SchnorrPublicKey<G>, ConstraintF>
    for SchnorrSigGadgetPk<G, ConstraintF, GG>
where
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
{
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(cs: CS, f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SchnorrPublicKey<G>>,
    {
        let pub_key = GG::alloc_input(cs, || f().map(|pk| *pk.borrow()))?;
        Ok(Self {
            pub_key,
            _engine: PhantomData,
            _group: PhantomData,
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<SchnorrPublicKey<G>>,
    {
        let pub_key = GG::alloc_input(cs, || f().map(|pk| *pk.borrow()))?;
        Ok(Self {
            pub_key,
            _engine: PhantomData,
            _group: PhantomData,
        })
    }
}

impl<G, ConstraintF, GG> ConditionalEqGadget<ConstraintF> for SchnorrSigGadgetPk<G, ConstraintF, GG>
where
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
{
    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.pub_key.conditional_enforce_equal(
            &mut cs.ns(|| "PubKey equality"),
            &other.pub_key,
            condition,
        )?;
        Ok(())
    }

    fn cost() -> usize {
        <GG as ConditionalEqGadget<ConstraintF>>::cost()
    }
}

impl<G, ConstraintF, GG> EqGadget<ConstraintF> for SchnorrSigGadgetPk<G, ConstraintF, GG>
where
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
{
}

impl<G, ConstraintF, GG> ToBytesGadget<ConstraintF> for SchnorrSigGadgetPk<G, ConstraintF, GG>
where
    G: Group,
    ConstraintF: Field,
    GG: GroupGadget<G, ConstraintF>,
{
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.pub_key.to_bytes(&mut cs.ns(|| "PubKey To Bytes"))
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.pub_key
            .to_bytes_strict(&mut cs.ns(|| "PubKey To Bytes"))
    }
}

//TODO: Replace to_bytes with to_bits
mod field_impl
{
    use algebra::{AffineCurve, PrimeField, ProjectiveCurve};
    use crate::{
        signature::{
            schnorr::field_impl::{FieldBasedSchnorrSignature, FieldBasedSchnorrSignatureScheme},
            FieldBasedSigGadget,
        },
        crh::{FieldBasedHashGadget, FieldBasedHash},
    };
    use r1cs_std::{
        fields::fp::FpGadget,
        alloc::AllocGadget,
        eq::EqGadget,
        groups::GroupGadget,
        bits::ToBitsGadget,
        ToBytesGadget,
        boolean::Boolean
    };
    use r1cs_core::{ConstraintSystem, SynthesisError};
    use std::{
        borrow::Borrow,
        marker::PhantomData,
    };
    use r1cs_std::bits::uint8::UInt8;
    use r1cs_std::groups::AffineGroupGadget;

    #[derive(Derivative)]
    #[derivative(
    Debug(bound = "ConstraintF: PrimeField"),
    Clone(bound = "ConstraintF: PrimeField"),
    PartialEq(bound = "ConstraintF: PrimeField"),
    Eq(bound = "ConstraintF: PrimeField")
    )]
    pub struct FieldBasedSchnorrSigGadget<
        ConstraintF: PrimeField,
    >
    {
        pub r:       FpGadget<ConstraintF>,
        pub s:       Vec<UInt8>,
        _field:      PhantomData<ConstraintF>,
    }

    impl<ConstraintF> AllocGadget<FieldBasedSchnorrSignature<ConstraintF>, ConstraintF>
    for FieldBasedSchnorrSigGadget<ConstraintF>
        where
            ConstraintF: PrimeField,
    {
        fn alloc<FN, T, CS: ConstraintSystem<ConstraintF>>(mut cs: CS, f: FN) -> Result<Self, SynthesisError>
            where
                FN: FnOnce() -> Result<T, SynthesisError>,
                T: Borrow<FieldBasedSchnorrSignature<ConstraintF>>,
        {
            f().and_then(|sig| {
                let FieldBasedSchnorrSignature {
                    r,
                    s
                } = sig.borrow().clone();
                let r = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc r"), || Ok(r))?;
                let s = UInt8::alloc_vec(cs.ns(|| "alloc s"), s.as_slice())?;
                Ok(Self{r, s, _field: PhantomData})
            })
        }

        fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(mut cs: CS, f: FN) -> Result<Self, SynthesisError>
            where
                FN: FnOnce() -> Result<T, SynthesisError>,
                T: Borrow<FieldBasedSchnorrSignature<ConstraintF>>,
        {
            f().and_then(|sig| {
                let FieldBasedSchnorrSignature {
                    r,
                    s
                } = sig.borrow().clone();
                let r = FpGadget::<ConstraintF>::alloc_input(cs.ns(|| "alloc r"), || Ok(r))?;
                let s = UInt8::alloc_input_vec(cs.ns(|| "alloc s"), s.as_slice())?;
                Ok(Self{r, s, _field: PhantomData})
            })
        }
    }

    pub struct FieldBasedSchnorrSigVerificationGadget<
        ConstraintF: PrimeField,
        G:  ProjectiveCurve<BaseField = ConstraintF>,
        GG: GroupGadget<G, ConstraintF>,
        H:  FieldBasedHash<Data = ConstraintF>,
        HG: FieldBasedHashGadget<H, ConstraintF>,
    >
    {
        _field:         PhantomData<ConstraintF>,
        _group:         PhantomData<G>,
        _group_gadget:  PhantomData<GG>,
        _hash:          PhantomData<H>,
        _hash_gadget:   PhantomData<HG>,
    }

    impl<ConstraintF, G, GG, H, HG> FieldBasedSigGadget<FieldBasedSchnorrSignatureScheme<ConstraintF, G, H>, ConstraintF>
    for FieldBasedSchnorrSigVerificationGadget<ConstraintF, G, GG, H, HG>
        where
            ConstraintF: PrimeField,
            G:           ProjectiveCurve<BaseField = ConstraintF>,
            G::Affine:   AffineCurve<BaseField = ConstraintF>,
            GG:          AffineGroupGadget<G, ConstraintF, FpGadget<ConstraintF>>,
            H:           FieldBasedHash<Data = ConstraintF>,
            HG:          FieldBasedHashGadget<H, ConstraintF, DataGadget = FpGadget<ConstraintF>>,
    {
        type DataGadget = FpGadget<ConstraintF>;
        type SignatureGadget = FieldBasedSchnorrSigGadget<ConstraintF>;
        type PublicKeyGadget = GG;

        fn check_verify_gadget<CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            public_key: &Self::PublicKeyGadget,
            signature: &Self::SignatureGadget,
            message: &[Self::DataGadget]
        ) -> Result<(), SynthesisError> {

            //Hardcode generator
            let g = GG::alloc_hardcoded(cs.ns(|| "hardcode generator"), || Ok(G::prime_subgroup_generator()))?;

            // Check e' = H(m || signature.r || pk)
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(message);
            hash_input.push(signature.r.clone());
            hash_input.push(public_key.get_x());
            hash_input.push(public_key.get_y());

            let e_prime = HG::check_evaluation_gadget(
                cs.ns(|| "check e_prime"),
                hash_input.as_slice()
            )?;

            //Enforce R' = s*G - e'*pk
            let e_prime_bits = e_prime
                .to_bytes(cs.ns(|| "e_prime_to_bytes"))?
                .to_bits(cs.ns(|| "e_prime_to_bits"))?;

            //We exploit hardcoded generator as `result` param here to avoid edge cases in addition
            let neg_e_prime_times_pk = public_key
                .mul_bits(cs.ns(|| "pk * e_prime + g"), &g, e_prime_bits.as_slice().iter())?
                .sub(cs.ns(|| "subtract g"), &g)?
                .negate(cs.ns(|| "- (e_prime * pk)"))?;

            //Enforce signature.s * G
            let s_bits = signature.s.to_bits(cs.ns(|| "s to bits"))?;

            //We must subtract the g previously added too.
            let r_prime = g
                .mul_bits_precomputed(cs.ns(|| "(s * G) - (e_prime * pk)"), &neg_e_prime_times_pk, s_bits.as_slice())?;

            //Enforce R'.x == r
            signature.r.enforce_equal(cs.ns(|| "sig.r == R'.x"), &r_prime.get_x())?;

            //Enforce R'.y is even
            let y_odd = r_prime.get_y().is_odd(cs.ns(|| "R'.y is_odd"))?;
            y_odd.enforce_equal(cs.ns(||"enforce R'.y not odd"), &Boolean::constant(false))?;
            Ok(())
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
    use crate::{
        signature::{
            FieldBasedSignatureScheme, FieldBasedSigGadget,
            schnorr::field_impl::*,
            schnorr::constraints::field_impl::*,
        },
        crh::{
            MNT4PoseidonHash, MNT4PoseidonHashGadget,
            MNT6PoseidonHash, MNT6PoseidonHashGadget,
        },
    };

    use r1cs_core::ConstraintSystem;
    use r1cs_std::alloc::AllocGadget;

    use r1cs_std::groups::curves::short_weierstrass::mnt::{
        mnt4::mnt4753::MNT4G1Gadget,
        mnt6::mnt6753::MNT6G1Gadget,
    };

    use rand::{Rng, thread_rng};
    use r1cs_std::test_constraint_system::TestConstraintSystem;

    type SchnorrMNT4 = FieldBasedSchnorrSignatureScheme<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash>;
    type SchnorrMNT6 = FieldBasedSchnorrSignatureScheme<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash>;

    type SchnorrMNT4Gadget = FieldBasedSchnorrSigVerificationGadget<
        MNT4Fr, MNT6G1Projective, MNT6G1Gadget, MNT4PoseidonHash, MNT4PoseidonHashGadget
    >;
    type SchnorrMNT6Gadget = FieldBasedSchnorrSigVerificationGadget<
        MNT6Fr, MNT4G1Projective, MNT4G1Gadget, MNT6PoseidonHash, MNT6PoseidonHashGadget
    >;

    fn sign<S: FieldBasedSignatureScheme, R: Rng>(rng: &mut R, message: &[S::Data]) -> (S::Signature, S::PublicKey, S::SecretKey)
    {
        let (pk, sk) = S::keygen(rng).unwrap();
        let sig = S::sign(&pk, &sk, &message).unwrap();
        (sig, pk, sk)
    }

    #[test]
    fn mnt4_schnorr_gadget_test() {

        let mut cs = TestConstraintSystem::<MNT4Fr>::new();

        //Sign a random field element f and get the signature and the public key
        let rng = &mut thread_rng();
        let message: MNT4Fr = rng.gen();
        let (sig, pk, sk) = sign::<SchnorrMNT4, _>(rng, &[message]);

        //Alloc signature, pk and message
        let sig_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc sig"),
            || Ok(sig)
        ).unwrap();
        let pk_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::PublicKeyGadget::alloc(cs.ns(|| "alloc pk"), || Ok(pk)).unwrap();
        let message_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::DataGadget::alloc(
            cs.ns(|| "alloc message"),
            || Ok(message)
        ).unwrap();

        //Verify sig
        SchnorrMNT4Gadget::check_verify_gadget(
            cs.ns(|| "verify sig1"),
            &pk_g,
            &sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(cs.is_satisfied());

        /*
        //Wrong message
        let message_new: MNT4Fr = rng.gen();
        let message_new_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::DataGadget::alloc(
            cs.ns(|| "alloc message_new"),
            || Ok(message_new)
        ).unwrap();

        SchnorrMNT4Gadget::check_verify_gadget(
            cs.ns(|| "verify sig2"),
            &pk_g,
            &sig_g,
            &[message_new_g]
        ).unwrap();*/

        //Wrong sig: generate a signature for a different message and check constraints fail
        let new_message: MNT4Fr = rng.gen();
        let new_sig = SchnorrMNT4::sign(&pk, &sk, &[new_message]).unwrap();
        let new_sig_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc new sig"),
            || Ok(new_sig)
        ).unwrap();

        //Verify new sig: expected to fail
        SchnorrMNT4Gadget::check_verify_gadget(
            cs.ns(|| "verify sig2"),
            &pk_g,
            &new_sig_g,
            &[message_g]
        ).unwrap();

        assert!(!cs.is_satisfied());
        println!("{:?}", cs.which_is_unsatisfied());
    }

    #[test]
    fn mnt6_schnorr_gadget_test() {

        let mut cs = TestConstraintSystem::<MNT6Fr>::new();

        //Sign a random field element f and get the signature and the public key
        let rng = &mut thread_rng();
        let message: MNT6Fr = rng.gen();
        let (sig, pk, sk) = sign::<SchnorrMNT6, _>(rng, &[message]);

        //Alloc signature, pk and message
        let sig_g = <SchnorrMNT6Gadget as FieldBasedSigGadget<SchnorrMNT6, MNT6Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc sig"),
            || Ok(sig)
        ).unwrap();
        let pk_g = <SchnorrMNT6Gadget as FieldBasedSigGadget<SchnorrMNT6, MNT6Fr>>::PublicKeyGadget::alloc(cs.ns(|| "alloc pk"), || Ok(pk)).unwrap();
        let message_g = <SchnorrMNT6Gadget as FieldBasedSigGadget<SchnorrMNT6, MNT6Fr>>::DataGadget::alloc(
            cs.ns(|| "alloc message"),
            || Ok(message)
        ).unwrap();

        //Verify sig
        SchnorrMNT6Gadget::check_verify_gadget(
            cs.ns(|| "verify sig1"),
            &pk_g,
            &sig_g,
            &[message_g.clone()]
        ).unwrap();

        println!("{:?}", cs.which_is_unsatisfied());

        assert!(cs.is_satisfied());

        /*
        //Wrong message
        let message_new: MNT6Fr = rng.gen();
        let message_new_g = <SchnorrMNT6Gadget as FieldBasedSigGadget<SchnorrMNT6, MNT6Fr>>::DataGadget::alloc(
            cs.ns(|| "alloc message_new"),
            || Ok(message_new)
        ).unwrap();

        SchnorrMNT6Gadget::check_verify_gadget(
            cs.ns(|| "verify sig2"),
            &pk_g,
            &sig_g,
            &[message_new_g]
        ).unwrap();*/

        //Wrong sig: let's generate a signature for a different message
        let new_message: MNT6Fr = rng.gen();
        let new_sig = SchnorrMNT6::sign(&pk, &sk, &[new_message]).unwrap();
        let new_sig_g = <SchnorrMNT6Gadget as FieldBasedSigGadget<SchnorrMNT6, MNT6Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc new sig"),
            || Ok(new_sig)
        ).unwrap();

        //Verify new sig: expected to fail
        SchnorrMNT6Gadget::check_verify_gadget(
            cs.ns(|| "verify sig2"),
            &pk_g,
            &new_sig_g,
            &[message_g]
        ).unwrap();

        assert!(!cs.is_satisfied());

        println!("{:?}", cs.which_is_unsatisfied());
    }

    #[test]
    fn random_schnorr_gadget_test() {

        //Sign a random field element f and get the signature and the public key
        let rng = &mut thread_rng();

        let samples = 10;
        for _ in 0..samples {
            let message: MNT4Fr = rng.gen();
            let (sig, pk, _) = sign::<SchnorrMNT4, _>(rng, &[message]);
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

            let message_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::DataGadget::alloc(
                cs.ns(|| "alloc message"),
                || Ok(message)
            ).unwrap();

            //Verify sig
            SchnorrMNT4Gadget::check_verify_gadget(
                cs.ns(|| "verify sig"),
                &pk_g,
                &sig_g,
                &[message_g.clone()]
            ).unwrap();
            assert!(cs.is_satisfied());
        }
    }
}