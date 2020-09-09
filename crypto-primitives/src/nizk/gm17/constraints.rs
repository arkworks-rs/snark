use crate::{
    nizk::{gm17::Gm17, NIZKVerifierGadget},
    Vec,
};
use algebra_core::{AffineCurve, PairingEngine, ToConstraintField};
use r1cs_core::{ConstraintSynthesizer, Namespace, SynthesisError};
use r1cs_std::prelude::*;

use core::{borrow::Borrow, marker::PhantomData};
use gm17::{PreparedVerifyingKey, Proof, VerifyingKey};

#[derive(Derivative)]
#[derivative(Clone(bound = "P::G1Var: Clone, P::G2Var: Clone"))]
pub struct ProofVar<E: PairingEngine, P: PairingVar<E>> {
    pub a: P::G1Var,
    pub b: P::G2Var,
    pub c: P::G1Var,
}

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P::G1Var: Clone, P::GTVar: Clone, P::G1PreparedVar: Clone, \
    P::G2PreparedVar: Clone, ")
)]
pub struct VerifyingKeyVar<E: PairingEngine, P: PairingVar<E>> {
    pub h_g2: P::G2Var,
    pub g_alpha_g1: P::G1Var,
    pub h_beta_g2: P::G2Var,
    pub g_gamma_g1: P::G1Var,
    pub h_gamma_g2: P::G2Var,
    pub query: Vec<P::G1Var>,
}

impl<E: PairingEngine, P: PairingVar<E>> VerifyingKeyVar<E, P> {
    #[tracing::instrument(target = "r1cs", skip(self))]
    pub fn prepare(&self) -> Result<PreparedVerifyingKeyVar<E, P>, SynthesisError> {
        let g_alpha_pc = P::prepare_g1(&self.g_alpha_g1)?;
        let h_beta_pc = P::prepare_g2(&self.h_beta_g2)?;
        let g_gamma_pc = P::prepare_g1(&self.g_gamma_g1)?;
        let h_gamma_pc = P::prepare_g2(&self.h_gamma_g2)?;
        let h_pc = P::prepare_g2(&self.h_g2)?;
        Ok(PreparedVerifyingKeyVar {
            g_alpha: self.g_alpha_g1.clone(),
            h_beta: self.h_beta_g2.clone(),
            g_alpha_pc,
            h_beta_pc,
            g_gamma_pc,
            h_gamma_pc,
            h_pc,
            query: self.query.clone(),
        })
    }
}

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P::G1Var: Clone, P::GTVar: Clone, P::G1PreparedVar: Clone, \
    P::G2PreparedVar: Clone, ")
)]
pub struct PreparedVerifyingKeyVar<E: PairingEngine, P: PairingVar<E>> {
    pub g_alpha: P::G1Var,
    pub h_beta: P::G2Var,
    pub g_alpha_pc: P::G1PreparedVar,
    pub h_beta_pc: P::G2PreparedVar,
    pub g_gamma_pc: P::G1PreparedVar,
    pub h_gamma_pc: P::G2PreparedVar,
    pub h_pc: P::G2PreparedVar,
    pub query: Vec<P::G1Var>,
}

pub struct Gm17VerifierGadget<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
{
    _pairing_engine: PhantomData<E>,
    _pairing_gadget: PhantomData<P>,
}

impl<E, P, C, V> NIZKVerifierGadget<Gm17<E, C, V>, E::Fq> for Gm17VerifierGadget<E, P>
where
    E: PairingEngine,
    C: ConstraintSynthesizer<E::Fr>,
    V: ToConstraintField<E::Fr>,
    P: PairingVar<E>,
{
    type PreparedVerificationKeyVar = PreparedVerifyingKeyVar<E, P>;
    type VerificationKeyVar = VerifyingKeyVar<E, P>;
    type ProofVar = ProofVar<E, P>;

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
                r1cs_core::ns!(cs, "Proof.a"),
                || Ok(proof.a.into_projective()),
                mode,
            )?;
            let b = CurveVar::new_variable_omit_prime_order_check(
                r1cs_core::ns!(cs, "Proof.b"),
                || Ok(proof.b.into_projective()),
                mode,
            )?;
            let c = CurveVar::new_variable_omit_prime_order_check(
                r1cs_core::ns!(cs, "Proof.c"),
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
    ) -> Result<Self::VerificationKeyVar, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        f().and_then(|vk| {
            let vk = vk.borrow();
            let g_alpha_g1 = P::G1Var::new_variable_omit_prime_order_check(
                r1cs_core::ns!(cs, "g_alpha"),
                || Ok(vk.g_alpha_g1.into_projective()),
                mode,
            )?;
            let h_g2 = P::G2Var::new_variable_omit_prime_order_check(
                r1cs_core::ns!(cs, "h"),
                || Ok(vk.h_g2.into_projective()),
                mode,
            )?;
            let h_beta_g2 = P::G2Var::new_variable_omit_prime_order_check(
                r1cs_core::ns!(cs, "h_beta"),
                || Ok(vk.h_beta_g2.into_projective()),
                mode,
            )?;
            let g_gamma_g1 = P::G1Var::new_variable_omit_prime_order_check(
                r1cs_core::ns!(cs, "g_gamma"),
                || Ok(vk.g_gamma_g1.into_projective()),
                mode,
            )?;
            let h_gamma_g2 = P::G2Var::new_variable_omit_prime_order_check(
                r1cs_core::ns!(cs, "h_gamma"),
                || Ok(vk.h_gamma_g2.into_projective()),
                mode,
            )?;

            let query = vk
                .query
                .iter()
                .map(|g| {
                    P::G1Var::new_variable_omit_prime_order_check(
                        r1cs_core::ns!(cs, "g"),
                        || Ok(g.into_projective()),
                        mode,
                    )
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(VerifyingKeyVar {
                g_alpha_g1,
                h_g2,
                h_beta_g2,
                g_gamma_g1,
                h_gamma_g2,
                query,
            })
        })
    }

    #[tracing::instrument(target = "r1cs", skip(vk, input, proof))]
    fn verify<'a, T: 'a + ToBitsGadget<E::Fq> + ?Sized>(
        vk: &Self::VerificationKeyVar,
        input: impl IntoIterator<Item = &'a T>,
        proof: &Self::ProofVar,
    ) -> Result<Boolean<E::Fq>, SynthesisError> {
        let pvk = vk.prepare()?;
        <Self as NIZKVerifierGadget<Gm17<E, C, V>, E::Fq>>::verify_prepared(&pvk, input, proof)
    }

    #[tracing::instrument(target = "r1cs", skip(pvk, input, proof))]
    fn verify_prepared<'a, T: 'a + ToBitsGadget<E::Fq> + ?Sized>(
        pvk: &Self::PreparedVerificationKeyVar,
        input: impl IntoIterator<Item = &'a T>,
        proof: &Self::ProofVar,
    ) -> Result<Boolean<E::Fq>, SynthesisError> {
        let pvk = pvk.clone();
        // e(A*G^{alpha}, B*H^{beta}) = e(G^{alpha}, H^{beta}) * e(G^{psi}, H^{gamma}) *
        // e(C, H) where psi = \sum_{i=0}^l input_i pvk.query[i]
        let g_psi = {
            let mut g_psi = pvk.query[0].clone();
            let mut input_len = 1;
            let mut input = input.into_iter();
            for (input, b) in input.by_ref().zip(pvk.query.iter().skip(1)) {
                let input_bits = input.to_bits_le()?;
                g_psi += b.scalar_mul_le(input_bits.iter())?;
                input_len += 1;
            }
            // Check that the input and the query in the verification are of the
            // same length.
            assert!(input_len == pvk.query.len() && input.next().is_none());
            g_psi
        };

        let mut test1_a_g_alpha = proof.a.clone();
        test1_a_g_alpha += pvk.g_alpha.clone();
        let mut test1_b_h_beta = proof.b.clone();
        test1_b_h_beta += pvk.h_beta.clone();

        let test1_exp = {
            test1_a_g_alpha = test1_a_g_alpha.negate()?;
            let test1_a_g_alpha_prep = P::prepare_g1(&test1_a_g_alpha)?;
            let test1_b_h_beta_prep = P::prepare_g2(&test1_b_h_beta)?;

            let g_psi_prep = P::prepare_g1(&g_psi)?;

            let c_prep = P::prepare_g1(&proof.c)?;

            P::miller_loop(
                &[
                    test1_a_g_alpha_prep,
                    g_psi_prep,
                    c_prep,
                    pvk.g_alpha_pc.clone(),
                ],
                &[
                    test1_b_h_beta_prep,
                    pvk.h_gamma_pc.clone(),
                    pvk.h_pc.clone(),
                    pvk.h_beta_pc.clone(),
                ],
            )?
        };

        let test1 = P::final_exponentiation(&test1_exp).unwrap();

        // e(A, H^{gamma}) = e(G^{gamma}, B)
        let test2_exp = {
            let a_prep = P::prepare_g1(&proof.a)?;
            // pvk.h_gamma_pc
            //&pvk.g_gamma_pc
            let proof_b = proof.b.negate()?;
            let b_prep = P::prepare_g2(&proof_b)?;
            P::miller_loop(&[a_prep, pvk.g_gamma_pc.clone()], &[pvk.h_gamma_pc, b_prep])?
        };
        let test2 = P::final_exponentiation(&test2_exp)?;

        let one = P::GTVar::one();
        test1.is_eq(&one)?.and(&test2.is_eq(&one)?)
    }
}

impl<E, P> AllocVar<PreparedVerifyingKey<E>, E::Fq> for PreparedVerifyingKeyVar<E, P>
where
    E: PairingEngine,
    P: PairingVar<E>,
    P::G1PreparedVar: AllocVar<E::G1Prepared, E::Fq>,
    P::G2PreparedVar: AllocVar<E::G2Prepared, E::Fq>,
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
            let g_alpha =
                P::G1Var::new_variable(r1cs_core::ns!(cs, "g_alpha"), || Ok(pvk.g_alpha), mode)?;
            let h_beta =
                P::G2Var::new_variable(r1cs_core::ns!(cs, "h_beta"), || Ok(pvk.h_beta), mode)?;
            let g_alpha_pc = P::G1PreparedVar::new_variable(
                r1cs_core::ns!(cs, "g_alpha_pc"),
                || Ok(pvk.g_alpha.into()),
                mode,
            )?;
            let h_beta_pc = P::G2PreparedVar::new_variable(
                r1cs_core::ns!(cs, "h_beta_pc"),
                || Ok(pvk.h_beta.into()),
                mode,
            )?;
            let g_gamma_pc = P::G1PreparedVar::new_variable(
                r1cs_core::ns!(cs, "g_gamma_pc"),
                || Ok(&pvk.g_gamma_pc),
                mode,
            )?;
            let h_gamma_pc = P::G2PreparedVar::new_variable(
                r1cs_core::ns!(cs, "h_gamma_pc"),
                || Ok(&pvk.h_gamma_pc),
                mode,
            )?;
            let h_pc =
                P::G2PreparedVar::new_variable(r1cs_core::ns!(cs, "h_pc"), || Ok(&pvk.h_pc), mode)?;
            let query =
                Vec::new_variable(r1cs_core::ns!(cs, "query"), || Ok(pvk.query.clone()), mode)?;

            Ok(Self {
                g_alpha,
                h_beta,
                g_alpha_pc,
                h_beta_pc,
                g_gamma_pc,
                h_gamma_pc,
                h_pc,
                query,
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
            let vk = vk.borrow();
            let g_alpha_g1 =
                P::G1Var::new_variable(r1cs_core::ns!(cs, "g_alpha"), || Ok(vk.g_alpha_g1), mode)?;
            let h_g2 = P::G2Var::new_variable(r1cs_core::ns!(cs, "h"), || Ok(vk.h_g2), mode)?;
            let h_beta_g2 =
                P::G2Var::new_variable(r1cs_core::ns!(cs, "h_beta"), || Ok(vk.h_beta_g2), mode)?;
            let g_gamma_g1 =
                P::G1Var::new_variable(r1cs_core::ns!(cs, "g_gamma"), || Ok(&vk.g_gamma_g1), mode)?;
            let h_gamma_g2 =
                P::G2Var::new_variable(r1cs_core::ns!(cs, "h_gamma"), || Ok(&vk.h_gamma_g2), mode)?;
            let query =
                Vec::new_variable(r1cs_core::ns!(cs, "query"), || Ok(vk.query.clone()), mode)?;
            Ok(Self {
                h_g2,
                g_alpha_g1,
                h_beta_g2,
                g_gamma_g1,
                h_gamma_g2,
                query,
            })
        })
    }
}

impl<E, P> AllocVar<Proof<E>, E::Fq> for ProofVar<E, P>
where
    E: PairingEngine,

    P: PairingVar<E>,
{
    #[inline]
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
            let a = P::G1Var::new_variable(cs.clone(), || Ok(a), mode)?;
            let b = P::G2Var::new_variable(cs.clone(), || Ok(b), mode)?;
            let c = P::G1Var::new_variable(r1cs_core::ns!(cs, "c"), || Ok(c), mode)?;
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
    fn to_bytes(&self) -> Result<Vec<UInt8<E::Fq>>, SynthesisError> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.h_g2.to_bytes()?);
        bytes.extend_from_slice(&self.g_alpha_g1.to_bytes()?);
        bytes.extend_from_slice(&self.h_beta_g2.to_bytes()?);
        bytes.extend_from_slice(&self.g_gamma_g1.to_bytes()?);
        bytes.extend_from_slice(&self.h_gamma_g2.to_bytes()?);
        for q in &self.query {
            bytes.extend_from_slice(&q.to_bytes()?);
        }
        Ok(bytes)
    }
}

#[cfg(test)]
mod test {
    use gm17::*;
    use r1cs_core::{
        lc, ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError,
    };

    use super::*;
    use algebra::{
        bls12_377::{Bls12_377, Fq, Fr},
        test_rng, BitIteratorLE, Field, PrimeField,
    };
    use r1cs_std::{bls12_377::PairingVar as Bls12_377PairingVar, boolean::Boolean, Assignment};
    use rand::Rng;

    type TestProofSystem = Gm17<Bls12_377, Bench<Fr>, Fr>;
    type TestVerifierGadget = Gm17VerifierGadget<Bls12_377, Bls12_377PairingVar>;
    type TestProofVar = ProofVar<Bls12_377, Bls12_377PairingVar>;
    type TestVkVar = VerifyingKeyVar<Bls12_377, Bls12_377PairingVar>;

    struct Bench<F: Field> {
        inputs: Vec<Option<F>>,
        num_constraints: usize,
    }

    impl<F: Field> ConstraintSynthesizer<F> for Bench<F> {
        fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
            assert!(self.inputs.len() >= 2);
            assert!(self.num_constraints >= self.inputs.len());

            let mut variables: Vec<_> = Vec::with_capacity(self.inputs.len());
            for input in self.inputs {
                let input_var = cs.new_input_variable(|| input.get())?;
                variables.push((input, input_var));
            }

            for i in 0..self.num_constraints {
                let new_entry = {
                    let (input_1_val, input_1_var) = variables[i];
                    let (input_2_val, input_2_var) = variables[i + 1];
                    let result_val = input_1_val
                        .and_then(|input_1| input_2_val.map(|input_2| input_1 * &input_2));
                    let result_var = cs.new_witness_variable(|| {
                        result_val.ok_or(SynthesisError::AssignmentMissing)
                    })?;
                    cs.enforce_constraint(
                        lc!() + input_1_var,
                        lc!() + input_2_var,
                        lc!() + result_var,
                    )
                    .unwrap();
                    (result_val, result_var)
                };
                variables.push(new_entry);
            }
            Ok(())
        }
    }

    #[test]
    fn gm17_verifier_test() {
        let num_inputs = 100;
        let num_constraints = num_inputs;
        let rng = &mut test_rng();
        let mut inputs: Vec<Option<Fr>> = Vec::with_capacity(num_inputs);
        for _ in 0..num_inputs {
            inputs.push(Some(rng.gen()));
        }
        let params = {
            let c = Bench::<Fr> {
                inputs: vec![None; num_inputs],
                num_constraints,
            };

            generate_random_parameters(c, rng).unwrap()
        };

        {
            let proof = {
                // Create an instance of our circuit (with the
                // witness)
                let c = Bench {
                    inputs: inputs.clone(),
                    num_constraints,
                };
                // Create a gm17 proof with our parameters.
                create_random_proof(c, &params, rng).unwrap()
            };

            // assert!(!verify_proof(&pvk, &proof, &[a]).unwrap());
            let cs = ConstraintSystem::<Fq>::new_ref();

            let inputs: Vec<_> = inputs.into_iter().map(|input| input.unwrap()).collect();
            let mut input_gadgets = Vec::new();

            {
                for input in inputs.into_iter() {
                    let input_bits: Vec<_> = BitIteratorLE::new(input.into_repr()).collect();
                    let input_bits =
                        Vec::<Boolean<Fq>>::new_input(r1cs_core::ns!(cs, "Input"), || {
                            Ok(input_bits)
                        })
                        .unwrap();
                    input_gadgets.push(input_bits);
                }
            }

            let vk_gadget =
                TestVkVar::new_input(r1cs_core::ns!(cs, "Vk"), || Ok(&params.vk)).unwrap();
            let proof_gadget =
                TestProofVar::new_witness(r1cs_core::ns!(cs, "Proof"), || Ok(proof.clone()))
                    .unwrap();
            println!("Time to verify!\n\n\n\n");
            <TestVerifierGadget as NIZKVerifierGadget<TestProofSystem, Fq>>::verify(
                &vk_gadget,
                &input_gadgets,
                &proof_gadget,
            )
            .unwrap()
            .enforce_equal(&Boolean::TRUE)
            .unwrap();

            if !cs.is_satisfied().unwrap() {
                println!("=========================================================");
                println!("Unsatisfied constraints:");
                println!("{:?}", cs.which_is_unsatisfied().unwrap());
                println!("=========================================================");
            }

            // cs.print_named_objects();
            assert!(cs.is_satisfied().unwrap());
        }
    }
}

#[cfg(test)]
mod test_recursive {
    use gm17::*;
    use r1cs_core::{
        lc, ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, SynthesisError,
    };

    use super::*;
    use algebra::{
        fields::{FftParameters, FpParameters},
        mnt4_298::{Fq as MNT4Fq, FqParameters as MNT4FqParameters, Fr as MNT4Fr, MNT4_298},
        mnt6_298::{Fq as MNT6Fq, FqParameters as MNT6FqParameters, Fr as MNT6Fr, MNT6_298},
        test_rng, BigInteger, Field, PrimeField,
    };
    use r1cs_std::{
        fields::fp::FpVar, mnt4_298::PairingVar as MNT4_298PairingVar,
        mnt6_298::PairingVar as MNT6_298PairingVar, uint8::UInt8, Assignment,
    };
    use rand::Rng;

    type TestProofSystem1 = Gm17<MNT6_298, Bench<MNT4Fq>, MNT6Fr>;
    type TestVerifierGadget1 = Gm17VerifierGadget<MNT6_298, MNT6_298PairingVar>;
    type TestProofVar1 = ProofVar<MNT6_298, MNT6_298PairingVar>;
    type TestVkVar1 = VerifyingKeyVar<MNT6_298, MNT6_298PairingVar>;

    type TestProofSystem2 = Gm17<MNT4_298, Wrapper, MNT4Fr>;
    type TestVerifierGadget2 = Gm17VerifierGadget<MNT4_298, MNT4_298PairingVar>;
    type TestProofVar2 = ProofVar<MNT4_298, MNT4_298PairingVar>;
    type TestVkVar2 = VerifyingKeyVar<MNT4_298, MNT4_298PairingVar>;

    #[derive(Clone)]
    struct Bench<F: Field> {
        inputs: Vec<Option<F>>,
        num_constraints: usize,
    }

    impl<F: Field> ConstraintSynthesizer<F> for Bench<F> {
        fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError> {
            assert!(self.inputs.len() >= 2);
            assert!(self.num_constraints >= self.inputs.len());

            let mut variables: Vec<_> = Vec::with_capacity(self.inputs.len());
            for input in self.inputs {
                let input_var = cs.new_input_variable(|| input.get())?;
                variables.push((input, input_var));
            }

            for i in 0..self.num_constraints {
                let new_entry = {
                    let (input_1_val, input_1_var) = variables[i];
                    let (input_2_val, input_2_var) = variables[i + 1];
                    let result_val = input_1_val
                        .and_then(|input_1| input_2_val.map(|input_2| input_1 * &input_2));
                    let result_var = cs.new_witness_variable(|| {
                        result_val.ok_or(SynthesisError::AssignmentMissing)
                    })?;
                    cs.enforce_constraint(
                        lc!() + input_1_var,
                        lc!() + input_2_var,
                        lc!() + result_var,
                    )
                    .unwrap();
                    (result_val, result_var)
                };
                variables.push(new_entry);
            }
            Ok(())
        }
    }

    struct Wrapper {
        inputs: Vec<Option<MNT4Fq>>,
        params: Parameters<MNT6_298>,
        proof: Proof<MNT6_298>,
    }

    impl ConstraintSynthesizer<MNT6Fq> for Wrapper {
        fn generate_constraints(
            self,
            cs: ConstraintSystemRef<MNT6Fq>,
        ) -> Result<(), SynthesisError> {
            let params = self.params;
            let proof = self.proof;
            let inputs: Vec<_> = self
                .inputs
                .into_iter()
                .map(|input| input.unwrap())
                .collect();
            let input_gadgets;

            {
                // Chain all input values in one large byte array.
                let input_bytes = inputs
                    .clone()
                    .into_iter()
                    .flat_map(|input| {
                        input
                            .into_repr()
                            .as_ref()
                            .iter()
                            .flat_map(|l| l.to_le_bytes().to_vec())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();

                // Allocate this byte array as input packed into field elements.
                let input_bytes =
                    UInt8::new_input_vec(r1cs_core::ns!(cs, "Input"), &input_bytes[..])?;
                // 40 byte
                let element_size = <MNT4FqParameters as FftParameters>::BigInt::NUM_LIMBS * 8;
                input_gadgets = input_bytes
                    .chunks(element_size)
                    .map(|chunk| {
                        chunk
                            .iter()
                            .flat_map(|byte| byte.to_bits_le().unwrap())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
            }

            let vk_gadget = TestVkVar1::new_witness(r1cs_core::ns!(cs, "Vk"), || Ok(&params.vk))?;
            let proof_gadget =
                TestProofVar1::new_witness(r1cs_core::ns!(cs, "Proof"), || Ok(proof.clone()))
                    .unwrap();
            <TestVerifierGadget1 as NIZKVerifierGadget<TestProofSystem1, MNT6Fq>>::verify(
                &vk_gadget,
                &input_gadgets,
                &proof_gadget,
            )?
            .enforce_equal(&Boolean::TRUE)?;
            Ok(())
        }
    }

    #[test]
    fn gm17_recursive_verifier_test() {
        let num_inputs = 5;
        let num_constraints = num_inputs;
        let rng = &mut test_rng();
        let mut inputs: Vec<Option<MNT4Fq>> = Vec::with_capacity(num_inputs);
        for _ in 0..num_inputs {
            inputs.push(Some(rng.gen()));
        }

        // Generate inner params and proof.
        let inner_params = {
            let c = Bench::<MNT4Fq> {
                inputs: vec![None; num_inputs],
                num_constraints,
            };

            generate_random_parameters(c, rng).unwrap()
        };

        let inner_proof = {
            // Create an instance of our circuit (with the
            // witness)
            let c = Bench {
                inputs: inputs.clone(),
                num_constraints,
            };
            // Create a gm17 proof with our parameters.
            create_random_proof(c, &inner_params, rng).unwrap()
        };

        // Generate outer params and proof.
        let params = {
            let c = Wrapper {
                inputs: inputs.clone(),
                params: inner_params.clone(),
                proof: inner_proof.clone(),
            };

            generate_random_parameters(c, rng).unwrap()
        };

        {
            let proof = {
                // Create an instance of our circuit (with the
                // witness)
                let c = Wrapper {
                    inputs: inputs.clone(),
                    params: inner_params.clone(),
                    proof: inner_proof.clone(),
                };
                // Create a gm17 proof with our parameters.
                create_random_proof(c, &params, rng).unwrap()
            };

            let cs = ConstraintSystem::<MNT4Fq>::new_ref();

            let inputs: Vec<_> = inputs.into_iter().map(|input| input.unwrap()).collect();
            let mut input_gadgets = Vec::new();

            {
                let bigint_size = <MNT4FqParameters as FftParameters>::BigInt::NUM_LIMBS * 64;
                let mut input_bits = Vec::new();
                for input in inputs.into_iter() {
                    let input_gadget =
                        FpVar::new_input(r1cs_core::ns!(cs, "Input"), || Ok(input)).unwrap();
                    let mut fp_bits = input_gadget.to_bits_le().unwrap();

                    // Use 320 bits per element.
                    for _ in fp_bits.len()..bigint_size {
                        fp_bits.push(Boolean::constant(false));
                    }
                    input_bits.extend_from_slice(&fp_bits);
                }

                // Pack input bits into field elements of the underlying circuit.
                let max_size = 8 * (<MNT6FqParameters as FpParameters>::CAPACITY / 8) as usize;
                let max_size = max_size as usize;
                let bigint_size = <MNT6FqParameters as FftParameters>::BigInt::NUM_LIMBS * 64;
                for chunk in input_bits.chunks(max_size) {
                    let mut chunk = chunk.to_vec();
                    let len = chunk.len();
                    for _ in len..bigint_size {
                        chunk.push(Boolean::constant(false));
                    }
                    input_gadgets.push(chunk);
                }
                // assert!(!verify_proof(&pvk, &proof, &[a]).unwrap());
            }

            let vk_gadget =
                TestVkVar2::new_input(r1cs_core::ns!(cs, "Vk"), || Ok(&params.vk)).unwrap();
            let proof_gadget =
                TestProofVar2::new_witness(r1cs_core::ns!(cs, "Proof"), || Ok(proof.clone()))
                    .unwrap();
            println!("Time to verify!\n\n\n\n");
            <TestVerifierGadget2 as NIZKVerifierGadget<TestProofSystem2, MNT4Fq>>::verify(
                &vk_gadget,
                &input_gadgets,
                &proof_gadget,
            )
            .unwrap()
            .enforce_equal(&Boolean::TRUE)
            .unwrap();
            if !cs.is_satisfied().unwrap() {
                println!("=========================================================");
                println!("Unsatisfied constraints:");
                println!("{:?}", cs.which_is_unsatisfied().unwrap());
                println!("=========================================================");
            }

            // cs.print_named_objects();
            assert!(cs.is_satisfied().unwrap());
        }
    }
}
