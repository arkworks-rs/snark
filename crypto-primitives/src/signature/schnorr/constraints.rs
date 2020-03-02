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

mod field_impl
{
    use algebra::{PrimeField, FpParameters, ProjectiveCurve, Group, ToConstraintField};
    use crate::{
        signature::{
            schnorr::field_impl::{FieldBasedSchnorrSignature, FieldBasedSchnorrSignatureScheme},
            FieldBasedSigGadget,
        },
        crh::{FieldBasedHashGadget, FieldBasedHash},
        compute_truncation_size
    };
    use r1cs_std::{
        fields::fp::FpGadget,
        alloc::AllocGadget,
        eq::EqGadget,
        groups::GroupGadget,
    };
    use r1cs_core::{ConstraintSystem, SynthesisError};
    use std::{
        borrow::Borrow,
        marker::PhantomData,
    };
    use r1cs_std::to_field_gadget_vec::ToConstraintFieldGadget;
    use r1cs_std::bits::boolean::Boolean;
    use r1cs_std::eq::EquVerdictGadget;

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
        pub e:       FpGadget<ConstraintF>,
        pub s:       FpGadget<ConstraintF>,
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
                    e,
                    s
                } = sig.borrow().clone();
                let e = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc e"), || Ok(e))?;
                let s = FpGadget::<ConstraintF>::alloc(cs.ns(|| "alloc s"), || Ok(s))?;
                Ok(Self{e, s, _field: PhantomData})
            })
        }

        fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(mut cs: CS, f: FN) -> Result<Self, SynthesisError>
            where
                FN: FnOnce() -> Result<T, SynthesisError>,
                T: Borrow<FieldBasedSchnorrSignature<ConstraintF>>,
        {
            f().and_then(|sig| {
                let FieldBasedSchnorrSignature {
                    e,
                    s
                } = sig.borrow().clone();
                let e = FpGadget::<ConstraintF>::alloc_input(cs.ns(|| "alloc e"), || Ok(e))?;
                let s = FpGadget::<ConstraintF>::alloc_input(cs.ns(|| "alloc s"), || Ok(s))?;
                Ok(Self{e, s, _field: PhantomData})
            })
        }
    }

    #[allow(dead_code)]
    pub struct FieldBasedSchnorrSigVerificationGadget<
        ConstraintF: PrimeField,
        G:  Group,
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
            G:           ProjectiveCurve + ToConstraintField<ConstraintF>,
            GG:          GroupGadget<G, ConstraintF, Value = G> + ToConstraintFieldGadget<ConstraintF, FieldGadget = HG::DataGadget>,
            H:           FieldBasedHash<Data = ConstraintF>,
            HG:          FieldBasedHashGadget<H, ConstraintF, DataGadget = FpGadget<ConstraintF>>,
    {
        type DataGadget = FpGadget<ConstraintF>;
        type SignatureGadget = FieldBasedSchnorrSigGadget<ConstraintF>;
        type PublicKeyGadget = GG;

        fn check_gadget<CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            public_key: &Self::PublicKeyGadget,
            signature: &Self::SignatureGadget,
            message: &[Self::DataGadget]
        ) -> Result<Boolean, SynthesisError> {

            //Enforce e' * pk
            let e_bits = {

                //Serialize e taking into account the length restriction
                let to_skip = compute_truncation_size(
                    ConstraintF::Params::MODULUS_BITS as i32,
                    <G::ScalarField as PrimeField>::Params::MODULUS_BITS as i32,
                );

                let e_bits = signature.e
                    .to_bits_with_length_restriction(cs.ns(|| "e_to_bits"), to_skip)?;

                debug_assert!(e_bits.len() as u32 == ConstraintF::Params::MODULUS_BITS - to_skip as u32);
                e_bits
            };

            //Let's hardcode generator and use it as `result` param here to avoid edge cases in addition
            let g = GG::alloc_hardcoded(cs.ns(|| "hardcode generator"), || Ok(G::prime_subgroup_generator()))?;
            let neg_e_times_pk = public_key
                .mul_bits(cs.ns(|| "pk * e + g"), &g, e_bits.as_slice().iter().rev())?
                .sub(cs.ns(|| "subtract g"), &g)?
                .negate(cs.ns(|| "- (e * pk)"))?;

            //Enforce s * G and R' = s*G - e*pk
            let mut s_bits = {

                //Serialize s taking into account the length restriction

                //Before computing the number of bits to truncate from s, we first have to normalize
                //it, i.e. considering its number of bits equals to G::ScalarField::MODULUS_BITS;
                let moduli_diff = ConstraintF::Params::MODULUS_BITS as i32 -
                    <G::ScalarField as PrimeField>::Params::MODULUS_BITS as i32;
                let to_skip_init = (if moduli_diff > 0 {moduli_diff} else {0}) as usize;

                //Now we can compare the two modulus and decide the bits to truncate
                let to_skip = to_skip_init + compute_truncation_size(
                    <G::ScalarField as PrimeField>::Params::MODULUS_BITS as i32,
                    ConstraintF::Params::MODULUS_BITS as i32,
                );

                let s_bits = signature.s
                    .to_bits_with_length_restriction(cs.ns(|| "s_to_bits"), to_skip as usize)?;

                debug_assert!(s_bits.len() as u32 == <G::ScalarField as PrimeField>::Params::MODULUS_BITS + to_skip_init as u32 - to_skip as u32);
                s_bits
            };

            s_bits.reverse();
            let r_prime = GG::mul_bits_precomputed(
                    &(g.get_value().unwrap()),
                    cs.ns(|| "(s * G) - (e * pk)"),
                    &neg_e_times_pk,
                    s_bits.as_slice()
                )?;

            let r_prime_x = {
                let r_prime_coords = r_prime.to_field_gadget_elements()?;
                r_prime_coords[0].clone()
            };

            // Check e' = H(m || signature.r || pk)
            let mut hash_input = Vec::new();
            hash_input.extend_from_slice(message);
            hash_input.push(r_prime_x);
            hash_input.extend_from_slice(public_key.to_field_gadget_elements().unwrap().as_slice());

            let e_prime = HG::check_evaluation_gadget(
                cs.ns(|| "check e_prime"),
                hash_input.as_slice()
            )?;

            //Enforce result of signature verification
            let is_verified = signature.e.enforce_verdict(cs.ns(|| "is e == e_prime"), &e_prime)?;

            Ok(is_verified)
        }

        fn check_verify_gadget<CS: ConstraintSystem<ConstraintF>>(
            mut cs: CS,
            public_key: &Self::PublicKeyGadget,
            signature: &Self::SignatureGadget,
            message: &[Self::DataGadget]
        ) -> Result<(), SynthesisError> {
            let is_verified = Self::check_gadget(cs.ns(|| "is sig verified"), public_key, signature, message)?;
            is_verified.enforce_equal(cs.ns(|| "signature must be verified"), &Boolean::Constant(true))?;
            Ok(())
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
    use crate::{
        signature::{
            FieldBasedSignatureScheme, FieldBasedSigGadget,
            schnorr::field_impl::*,
            schnorr::constraints::field_impl::*,
        },
        crh::{
            MNT4PoseidonHash, MNT4PoseidonHashGadget,
            MNT6PoseidonHash, MNT6PoseidonHashGadget,
            BLS12PoseidonHash, BLS12PoseidonHashGadget,
        },
    };

    use r1cs_core::ConstraintSystem;
    use r1cs_std::alloc::AllocGadget;

    use r1cs_std::groups::curves::short_weierstrass::mnt::{
        mnt4::mnt4753::MNT4G1Gadget,
        mnt6::mnt6753::MNT6G1Gadget,
    };

    use r1cs_std::groups::curves::twisted_edwards::jubjub::JubJubGadget;

    use rand::{Rng, thread_rng};
    use r1cs_std::test_constraint_system::TestConstraintSystem;

    type SchnorrMNT4 = FieldBasedSchnorrSignatureScheme<MNT4Fr, MNT6G1Projective, MNT4PoseidonHash>;
    type SchnorrMNT6 = FieldBasedSchnorrSignatureScheme<MNT6Fr, MNT4G1Projective, MNT6PoseidonHash>;
    type SchnorrBLS12 = FieldBasedSchnorrSignatureScheme<BLS12Fr, JubJubProjective, BLS12PoseidonHash>;

    type SchnorrMNT4Gadget = FieldBasedSchnorrSigVerificationGadget<
        MNT4Fr, MNT6G1Projective, MNT6G1Gadget, MNT4PoseidonHash, MNT4PoseidonHashGadget
    >;
    type SchnorrMNT6Gadget = FieldBasedSchnorrSigVerificationGadget<
        MNT6Fr, MNT4G1Projective, MNT4G1Gadget, MNT6PoseidonHash, MNT6PoseidonHashGadget
    >;
    type SchnorrBLS12Gadget = FieldBasedSchnorrSigVerificationGadget<
        BLS12Fr, JubJubProjective, JubJubGadget, BLS12PoseidonHash, BLS12PoseidonHashGadget
    >;

    fn sign<S: FieldBasedSignatureScheme, R: Rng>(rng: &mut R, message: &[S::Data]) -> (S::Signature, S::PublicKey, S::SecretKey)
    {
        let (pk, sk) = S::keygen(rng).unwrap();
        let sig = S::sign(rng, &pk, &sk, &message).unwrap();
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
        let is_verified = SchnorrMNT4Gadget::check_gadget(
            cs.ns(|| "sig1 result"),
            &pk_g,
            &sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(is_verified.get_value().unwrap());

        SchnorrMNT4Gadget::check_verify_gadget(
            cs.ns(|| "verify sig1"),
            &pk_g,
            &sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(cs.is_satisfied());

        //Wrong sig: check constraints fail
        let new_message: MNT4Fr = rng.gen();
        let new_sig = SchnorrMNT4::sign(rng, &pk, &sk, &[new_message]).unwrap();
        let new_sig_g = <SchnorrMNT4Gadget as FieldBasedSigGadget<SchnorrMNT4, MNT4Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc new sig"),
            || Ok(new_sig)
        ).unwrap();

        //Verify new sig: expected to fail
        let is_verified = SchnorrMNT4Gadget::check_gadget(
            cs.ns(|| "sig2 result"),
            &pk_g,
            &new_sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(!is_verified.get_value().unwrap());
        assert!(cs.is_satisfied());

        SchnorrMNT4Gadget::check_verify_gadget(
            cs.ns(|| "verify sig2"),
            &pk_g,
            &new_sig_g,
            &[message_g]
        ).unwrap();

        println!("{:?}", cs.which_is_unsatisfied());
        assert!(!cs.is_satisfied());
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
        let is_verified = SchnorrMNT6Gadget::check_gadget(
            cs.ns(|| "sig1 result"),
            &pk_g,
            &sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(is_verified.get_value().unwrap());

        SchnorrMNT6Gadget::check_verify_gadget(
            cs.ns(|| "verify sig1"),
            &pk_g,
            &sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(cs.is_satisfied());

        //Wrong sig: check constraints fail
        let new_message: MNT6Fr = rng.gen();
        let new_sig = SchnorrMNT6::sign(rng, &pk, &sk, &[new_message]).unwrap();
        let new_sig_g = <SchnorrMNT6Gadget as FieldBasedSigGadget<SchnorrMNT6, MNT6Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc new sig"),
            || Ok(new_sig)
        ).unwrap();

        //Verify new sig: expected to fail
        let is_verified = SchnorrMNT6Gadget::check_gadget(
            cs.ns(|| "sig2 result"),
            &pk_g,
            &new_sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(!is_verified.get_value().unwrap());
        assert!(cs.is_satisfied());

        SchnorrMNT6Gadget::check_verify_gadget(
            cs.ns(|| "verify sig2"),
            &pk_g,
            &new_sig_g,
            &[message_g]
        ).unwrap();

        println!("{:?}", cs.which_is_unsatisfied());
        assert!(!cs.is_satisfied());
    }

    #[test]
    fn bls12_381_schnorr_gadget_test() {

        let mut cs = TestConstraintSystem::<BLS12Fr>::new();

        //Sign a random field element f and get the signature and the public key
        let rng = &mut thread_rng();
        let message: BLS12Fr = rng.gen();
        let (sig, pk, sk) = sign::<SchnorrBLS12, _>(rng, &[message]);

        //Alloc signature, pk and message
        let sig_g = <SchnorrBLS12Gadget as FieldBasedSigGadget<SchnorrBLS12, BLS12Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc sig"),
            || Ok(sig)
        ).unwrap();
        let pk_g = <SchnorrBLS12Gadget as FieldBasedSigGadget<SchnorrBLS12, BLS12Fr>>::PublicKeyGadget::alloc(cs.ns(|| "alloc pk"), || Ok(pk)).unwrap();
        let message_g = <SchnorrBLS12Gadget as FieldBasedSigGadget<SchnorrBLS12, BLS12Fr>>::DataGadget::alloc(
            cs.ns(|| "alloc message"),
            || Ok(message)
        ).unwrap();

        //Verify sig
        let is_verified = SchnorrBLS12Gadget::check_gadget(
            cs.ns(|| "sig1 result"),
            &pk_g,
            &sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(is_verified.get_value().unwrap());

        SchnorrBLS12Gadget::check_verify_gadget(
            cs.ns(|| "verify sig1"),
            &pk_g,
            &sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(cs.is_satisfied());

        //Wrong sig: check constraints fail
        let new_message: BLS12Fr = rng.gen();
        let new_sig = SchnorrBLS12::sign(rng, &pk, &sk, &[new_message]).unwrap();
        let new_sig_g = <SchnorrBLS12Gadget as FieldBasedSigGadget<SchnorrBLS12, BLS12Fr>>::SignatureGadget::alloc(
            cs.ns(|| "alloc new sig"),
            || Ok(new_sig)
        ).unwrap();

        //Verify new sig: expected to fail
        let is_verified = SchnorrBLS12Gadget::check_gadget(
            cs.ns(|| "sig2 result"),
            &pk_g,
            &new_sig_g,
            &[message_g.clone()]
        ).unwrap();

        assert!(!is_verified.get_value().unwrap());
        assert!(cs.is_satisfied());

        SchnorrBLS12Gadget::check_verify_gadget(
            cs.ns(|| "verify sig2"),
            &pk_g,
            &new_sig_g,
            &[message_g]
        ).unwrap();

        println!("{:?}", cs.which_is_unsatisfied());
        assert!(!cs.is_satisfied());
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
            let is_verified = SchnorrMNT4Gadget::check_gadget(
                cs.ns(|| "sig result"),
                &pk_g,
                &sig_g,
                &[message_g.clone()]
            ).unwrap();

            assert!(is_verified.get_value().unwrap());

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