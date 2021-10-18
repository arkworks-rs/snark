use crate::{Groth16, PreparedVerifyingKey, Proof, VerifyingKey};
use ark_crypto_primitives::snark::constraints::{CircuitSpecificSetupSNARKGadget, SNARKGadget};
use ark_crypto_primitives::snark::{BooleanInputVar, SNARK};
use ark_ec::{AffineCurve, PairingEngine};
use ark_r1cs_std::groups::CurveVar;
use ark_r1cs_std::{
    alloc::{AllocVar, AllocationMode},
    bits::boolean::Boolean,
    bits::uint8::UInt8,
    eq::EqGadget,
    pairing::PairingVar,
    ToBitsGadget, ToBytesGadget,
};
use ark_relations::r1cs::{Namespace, SynthesisError};
use ark_std::{borrow::Borrow, marker::PhantomData, vec::Vec};

/// The proof variable for the Groth16 construction
#[derive(Derivative)]
#[derivative(Clone(bound = "P::G1Var: Clone, P::G2Var: Clone"))]
pub struct ProofVar<E: PairingEngine, P: PairingVar<E>> {
    /// The `A` element in `G1`.
    pub a: P::G1Var,
    /// The `B` element in `G2`.
    pub b: P::G2Var,
    /// The `C` element in `G1`.
    pub c: P::G1Var,
}

/// A variable representing the Groth16 verifying key in the constraint system.
#[derive(Derivative)]
#[derivative(
    Clone(bound = "P::G1Var: Clone, P::GTVar: Clone, P::G1PreparedVar: Clone, \
    P::G2PreparedVar: Clone, ")
)]
pub struct VerifyingKeyVar<E: PairingEngine, P: PairingVar<E>> {
    #[doc(hidden)]
    pub alpha_g1: P::G1Var,
    #[doc(hidden)]
    pub beta_g2: P::G2Var,
    #[doc(hidden)]
    pub gamma_g2: P::G2Var,
    #[doc(hidden)]
    pub delta_g2: P::G2Var,
    #[doc(hidden)]
    pub gamma_abc_g1: Vec<P::G1Var>,
}

impl<E: PairingEngine, P: PairingVar<E>> VerifyingKeyVar<E, P> {
    /// Prepare `self` for use in proof verification.
    pub fn prepare(&self) -> Result<PreparedVerifyingKeyVar<E, P>, SynthesisError> {
        let alpha_g1_pc = P::prepare_g1(&self.alpha_g1)?;
        let beta_g2_pc = P::prepare_g2(&self.beta_g2)?;

        let alpha_g1_beta_g2 = P::pairing(alpha_g1_pc, beta_g2_pc)?;
        let gamma_g2_neg_pc = P::prepare_g2(&self.gamma_g2.negate()?)?;
        let delta_g2_neg_pc = P::prepare_g2(&self.delta_g2.negate()?)?;

        Ok(PreparedVerifyingKeyVar {
            alpha_g1_beta_g2,
            gamma_g2_neg_pc,
            delta_g2_neg_pc,
            gamma_abc_g1: self.gamma_abc_g1.clone(),
        })
    }
}

/// Preprocessed verification key parameters variable for the Groth16 construction
#[derive(Derivative)]
#[derivative(
    Clone(bound = "P::G1Var: Clone, P::GTVar: Clone, P::G1PreparedVar: Clone, \
    P::G2PreparedVar: Clone, ")
)]
pub struct PreparedVerifyingKeyVar<E: PairingEngine, P: PairingVar<E>> {
    #[doc(hidden)]
    pub alpha_g1_beta_g2: P::GTVar,
    #[doc(hidden)]
    pub gamma_g2_neg_pc: P::G2PreparedVar,
    #[doc(hidden)]
    pub delta_g2_neg_pc: P::G2PreparedVar,
    #[doc(hidden)]
    pub gamma_abc_g1: Vec<P::G1Var>,
}

/// Constraints for the verifier of the SNARK of [[Groth16]](https://eprint.iacr.org/2016/260.pdf).
pub struct Groth16VerifierGadget<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    _pairing_engine: PhantomData<E>,
    _pairing_gadget: PhantomData<P>,
}

impl<E: PairingEngine, P: PairingVar<E, E::Fq>> SNARKGadget<E::Fr, E::Fq, Groth16<E>>
    for Groth16VerifierGadget<E, P>
{
    type ProcessedVerifyingKeyVar = PreparedVerifyingKeyVar<E, P>;
    type VerifyingKeyVar = VerifyingKeyVar<E, P>;
    type InputVar = BooleanInputVar<E::Fr, E::Fq>;
    type ProofVar = ProofVar<E, P>;

    type VerifierSize = usize;

    fn verifier_size(
        circuit_vk: &<Groth16<E> as SNARK<E::Fr>>::VerifyingKey,
    ) -> Self::VerifierSize {
        circuit_vk.gamma_abc_g1.len()
    }

    /// Allocates `N::Proof` in `cs` without performing
    /// subgroup checks.
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_proof_unchecked<T: Borrow<Proof<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self::ProofVar, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        f().and_then(|proof| {
            let proof = proof.borrow();
            let a = CurveVar::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "Proof.a"),
                || Ok(proof.a.into_projective()),
                mode,
            )?;
            let b = CurveVar::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "Proof.b"),
                || Ok(proof.b.into_projective()),
                mode,
            )?;
            let c = CurveVar::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "Proof.c"),
                || Ok(proof.c.into_projective()),
                mode,
            )?;
            Ok(ProofVar { a, b, c })
        })
    }

    /// Allocates `N::Proof` in `cs` without performing
    /// subgroup checks.
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_verification_key_unchecked<T: Borrow<VerifyingKey<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self::VerifyingKeyVar, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        f().and_then(|vk| {
            let vk = vk.borrow();
            let alpha_g1 = P::G1Var::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "alpha_g1"),
                || Ok(vk.alpha_g1.into_projective()),
                mode,
            )?;
            let beta_g2 = P::G2Var::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "beta_g2"),
                || Ok(vk.beta_g2.into_projective()),
                mode,
            )?;
            let gamma_g2 = P::G2Var::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "gamma_g2"),
                || Ok(vk.gamma_g2.into_projective()),
                mode,
            )?;
            let delta_g2 = P::G2Var::new_variable_omit_prime_order_check(
                ark_relations::ns!(cs, "delta_g2"),
                || Ok(vk.delta_g2.into_projective()),
                mode,
            )?;
            let gamma_abc_g1 = vk
                .gamma_abc_g1
                .iter()
                .map(|g| {
                    P::G1Var::new_variable_omit_prime_order_check(
                        ark_relations::ns!(cs, "gamma_abc_g1"),
                        || Ok(g.into_projective()),
                        mode,
                    )
                })
                .collect::<Result<Vec<_>, _>>()?;

            Ok(VerifyingKeyVar {
                alpha_g1,
                beta_g2,
                gamma_g2,
                delta_g2,
                gamma_abc_g1,
            })
        })
    }

    #[tracing::instrument(target = "r1cs", skip(circuit_pvk, x, proof))]
    fn verify_with_processed_vk(
        circuit_pvk: &Self::ProcessedVerifyingKeyVar,
        x: &Self::InputVar,
        proof: &Self::ProofVar,
    ) -> Result<Boolean<E::Fq>, SynthesisError> {
        let circuit_pvk = circuit_pvk.clone();

        let g_ic = {
            let mut g_ic: P::G1Var = circuit_pvk.gamma_abc_g1[0].clone();
            let mut input_len = 1;
            let mut public_inputs = x.clone().into_iter();
            for (input, b) in public_inputs
                .by_ref()
                .zip(circuit_pvk.gamma_abc_g1.iter().skip(1))
            {
                let encoded_input_i: P::G1Var = b.scalar_mul_le(input.to_bits_le()?.iter())?;
                g_ic += encoded_input_i;
                input_len += 1;
            }
            // Check that the input and the query in the verification are of the
            // same length.
            assert!(input_len == circuit_pvk.gamma_abc_g1.len() && public_inputs.next().is_none());
            g_ic
        };

        let test_exp = {
            let proof_a_prep = P::prepare_g1(&proof.a)?;
            let proof_b_prep = P::prepare_g2(&proof.b)?;
            let proof_c_prep = P::prepare_g1(&proof.c)?;

            let g_ic_prep = P::prepare_g1(&g_ic)?;

            P::miller_loop(
                &[proof_a_prep, g_ic_prep, proof_c_prep],
                &[
                    proof_b_prep,
                    circuit_pvk.gamma_g2_neg_pc.clone(),
                    circuit_pvk.delta_g2_neg_pc.clone(),
                ],
            )?
        };

        let test = P::final_exponentiation(&test_exp)?;
        test.is_eq(&circuit_pvk.alpha_g1_beta_g2)
    }

    #[tracing::instrument(target = "r1cs", skip(circuit_vk, x, proof))]
    fn verify(
        circuit_vk: &Self::VerifyingKeyVar,
        x: &Self::InputVar,
        proof: &Self::ProofVar,
    ) -> Result<Boolean<E::Fq>, SynthesisError> {
        let pvk = circuit_vk.prepare()?;
        Self::verify_with_processed_vk(&pvk, x, proof)
    }
}

impl<E, P> CircuitSpecificSetupSNARKGadget<E::Fr, E::Fq, Groth16<E>> for Groth16VerifierGadget<E, P>
where
    E: PairingEngine,
    P: PairingVar<E, E::Fq>,
{
}

impl<E, P> AllocVar<PreparedVerifyingKey<E>, E::Fq> for PreparedVerifyingKeyVar<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T: Borrow<PreparedVerifyingKey<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        f().and_then(|pvk| {
            let pvk = pvk.borrow();
            let alpha_g1_beta_g2 = P::GTVar::new_variable(
                ark_relations::ns!(cs, "alpha_g1_beta_g2"),
                || Ok(pvk.alpha_g1_beta_g2.clone()),
                mode,
            )?;

            let gamma_g2_neg_pc = P::G2PreparedVar::new_variable(
                ark_relations::ns!(cs, "gamma_g2_neg_pc"),
                || Ok(pvk.gamma_g2_neg_pc.clone()),
                mode,
            )?;

            let delta_g2_neg_pc = P::G2PreparedVar::new_variable(
                ark_relations::ns!(cs, "delta_g2_neg_pc"),
                || Ok(pvk.delta_g2_neg_pc.clone()),
                mode,
            )?;

            let gamma_abc_g1 = Vec::new_variable(
                ark_relations::ns!(cs, "gamma_abc_g1"),
                || Ok(pvk.vk.gamma_abc_g1.clone()),
                mode,
            )?;

            Ok(Self {
                alpha_g1_beta_g2,
                gamma_g2_neg_pc,
                delta_g2_neg_pc,
                gamma_abc_g1,
            })
        })
    }
}

impl<E, P> AllocVar<VerifyingKey<E>, E::Fq> for VerifyingKeyVar<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T: Borrow<VerifyingKey<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        f().and_then(|vk| {
            let VerifyingKey {
                alpha_g1,
                beta_g2,
                gamma_g2,
                delta_g2,
                gamma_abc_g1,
            } = vk.borrow().clone();
            let alpha_g1 =
                P::G1Var::new_variable(ark_relations::ns!(cs, "alpha_g1"), || Ok(alpha_g1), mode)?;
            let beta_g2 =
                P::G2Var::new_variable(ark_relations::ns!(cs, "beta_g2"), || Ok(beta_g2), mode)?;
            let gamma_g2 =
                P::G2Var::new_variable(ark_relations::ns!(cs, "gamma_g2"), || Ok(gamma_g2), mode)?;
            let delta_g2 =
                P::G2Var::new_variable(ark_relations::ns!(cs, "delta_g2"), || Ok(delta_g2), mode)?;

            let gamma_abc_g1 = Vec::new_variable(cs.clone(), || Ok(gamma_abc_g1), mode)?;
            Ok(Self {
                alpha_g1,
                beta_g2,
                gamma_g2,
                delta_g2,
                gamma_abc_g1,
            })
        })
    }
}

impl<E, P> AllocVar<Proof<E>, E::Fq> for ProofVar<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    #[tracing::instrument(target = "r1cs", skip(cs, f))]
    fn new_variable<T: Borrow<Proof<E>>>(
        cs: impl Into<Namespace<E::Fq>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        f().and_then(|proof| {
            let Proof { a, b, c } = proof.borrow().clone();
            let a = P::G1Var::new_variable(ark_relations::ns!(cs, "a"), || Ok(a), mode)?;
            let b = P::G2Var::new_variable(ark_relations::ns!(cs, "b"), || Ok(b), mode)?;
            let c = P::G1Var::new_variable(ark_relations::ns!(cs, "c"), || Ok(c), mode)?;
            Ok(Self { a, b, c })
        })
    }
}

impl<E, P> ToBytesGadget<E::Fq> for VerifyingKeyVar<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    #[inline]
    #[tracing::instrument(target = "r1cs", skip(self))]
    fn to_bytes(&self) -> Result<Vec<UInt8<E::Fq>>, SynthesisError> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.alpha_g1.to_bytes()?);
        bytes.extend_from_slice(&self.beta_g2.to_bytes()?);
        bytes.extend_from_slice(&self.gamma_g2.to_bytes()?);
        bytes.extend_from_slice(&self.delta_g2.to_bytes()?);
        for g in &self.gamma_abc_g1 {
            bytes.extend_from_slice(&g.to_bytes()?);
        }
        Ok(bytes)
    }
}

#[cfg(test)]
mod test {
    use crate::{constraints::Groth16VerifierGadget, Groth16};
    use ark_crypto_primitives::snark::constraints::SNARKGadget;
    use ark_crypto_primitives::snark::{CircuitSpecificSetupSNARK, SNARK};
    use ark_ec::PairingEngine;
    use ark_ff::{Field, UniformRand};
    use ark_mnt4_298::{
        constraints::PairingVar as MNT4PairingVar, Fr as MNT4Fr, MNT4_298 as MNT4PairingEngine,
    };
    use ark_mnt6_298::Fr as MNT6Fr;
    use ark_r1cs_std::bits::boolean::Boolean;
    use ark_r1cs_std::{alloc::AllocVar, eq::EqGadget};
    use ark_relations::{
        lc, ns,
        r1cs::{ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError},
    };
    use ark_std::ops::MulAssign;
    use ark_std::test_rng;

    #[derive(Copy, Clone)]
    struct Circuit<F: Field> {
        a: Option<F>,
        b: Option<F>,
        num_constraints: usize,
        num_variables: usize,
    }

    impl<ConstraintF: Field> ConstraintSynthesizer<ConstraintF> for Circuit<ConstraintF> {
        fn generate_constraints(
            self,
            cs: ConstraintSystemRef<ConstraintF>,
        ) -> Result<(), SynthesisError> {
            let a = cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
            let b = cs.new_witness_variable(|| self.b.ok_or(SynthesisError::AssignmentMissing))?;
            let c = cs.new_input_variable(|| {
                let mut a = self.a.ok_or(SynthesisError::AssignmentMissing)?;
                let b = self.b.ok_or(SynthesisError::AssignmentMissing)?;

                a.mul_assign(&b);
                Ok(a)
            })?;

            for _ in 0..(self.num_variables - 3) {
                cs.new_witness_variable(|| self.a.ok_or(SynthesisError::AssignmentMissing))?;
            }

            for _ in 0..self.num_constraints {
                cs.enforce_constraint(lc!() + a, lc!() + b, lc!() + c)
                    .unwrap();
            }
            Ok(())
        }
    }

    type TestSNARK = Groth16<MNT4PairingEngine>;
    type TestSNARKGadget = Groth16VerifierGadget<MNT4PairingEngine, MNT4PairingVar>;

    #[test]
    fn groth16_snark_test() {
        let mut rng = test_rng();
        let a = MNT4Fr::rand(&mut rng);
        let b = MNT4Fr::rand(&mut rng);
        let mut c = a;
        c.mul_assign(&b);

        let circ = Circuit {
            a: Some(a.clone()),
            b: Some(b.clone()),
            num_constraints: 100,
            num_variables: 25,
        };

        let (pk, vk) = TestSNARK::setup(circ, &mut rng).unwrap();

        let proof = TestSNARK::prove(&pk, circ.clone(), &mut rng).unwrap();

        assert!(
            TestSNARK::verify(&vk, &vec![c], &proof).unwrap(),
            "The native verification check fails."
        );

        let cs_sys = ConstraintSystem::<MNT6Fr>::new();
        let cs = ConstraintSystemRef::new(cs_sys);

        let input_gadget = <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::InputVar::new_input(ns!(cs, "new_input"), || Ok(vec![c]))
        .unwrap();
        let proof_gadget = <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::ProofVar::new_witness(ns!(cs, "alloc_proof"), || Ok(proof))
        .unwrap();
        let vk_gadget = <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::VerifyingKeyVar::new_constant(ns!(cs, "alloc_vk"), vk.clone())
        .unwrap();
        <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::verify(&vk_gadget, &input_gadget, &proof_gadget)
        .unwrap()
        .enforce_equal(&Boolean::constant(true))
        .unwrap();

        assert!(
            cs.is_satisfied().unwrap(),
            "Constraints not satisfied: {}",
            cs.which_is_unsatisfied().unwrap().unwrap_or_default()
        );

        let pvk = TestSNARK::process_vk(&vk).unwrap();
        let pvk_gadget = <TestSNARKGadget as SNARKGadget<
            <MNT4PairingEngine as PairingEngine>::Fr,
            <MNT4PairingEngine as PairingEngine>::Fq,
            TestSNARK,
        >>::ProcessedVerifyingKeyVar::new_constant(
            ns!(cs, "alloc_pvk"), pvk.clone()
        )
        .unwrap();
        TestSNARKGadget::verify_with_processed_vk(&pvk_gadget, &input_gadget, &proof_gadget)
            .unwrap()
            .enforce_equal(&Boolean::constant(true))
            .unwrap();

        assert!(
            cs.is_satisfied().unwrap(),
            "Constraints not satisfied: {}",
            cs.which_is_unsatisfied().unwrap().unwrap_or_default()
        );
    }
}
