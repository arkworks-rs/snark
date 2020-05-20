use crate::{
    nizk::{gm17::Gm17, NIZKVerifierGadget},
    Vec,
};
use algebra_core::{AffineCurve, Field, PairingEngine, ToConstraintField};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use core::{borrow::Borrow, marker::PhantomData};
use gm17::{Proof, VerifyingKey};

#[derive(Derivative)]
#[derivative(Clone(bound = "P::G1Gadget: Clone, P::G2Gadget: Clone"))]
pub struct ProofGadget<
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
> {
    pub a: P::G1Gadget,
    pub b: P::G2Gadget,
    pub c: P::G1Gadget,
}

#[derive(Derivative)]
#[derivative(Clone(
    bound = "P::G1Gadget: Clone, P::GTGadget: Clone, P::G1PreparedGadget: Clone, \
             P::G2PreparedGadget: Clone, "
))]
pub struct VerifyingKeyGadget<
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
> {
    pub h_g2: P::G2Gadget,
    pub g_alpha_g1: P::G1Gadget,
    pub h_beta_g2: P::G2Gadget,
    pub g_gamma_g1: P::G1Gadget,
    pub h_gamma_g2: P::G2Gadget,
    pub query: Vec<P::G1Gadget>,
}

impl<PairingE: PairingEngine, ConstraintF: Field, P: PairingGadget<PairingE, ConstraintF>>
    VerifyingKeyGadget<PairingE, ConstraintF, P>
{
    pub fn prepare<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<PreparedVerifyingKeyGadget<PairingE, ConstraintF, P>, SynthesisError> {
        let mut cs = cs.ns(|| "Preparing verifying key");
        let g_alpha_pc = P::prepare_g1(&mut cs.ns(|| "Prepare g_alpha_g1"), &self.g_alpha_g1)?;
        let h_beta_pc = P::prepare_g2(&mut cs.ns(|| "Prepare h_beta_g2"), &self.h_beta_g2)?;
        let g_gamma_pc = P::prepare_g1(&mut cs.ns(|| "Prepare g_gamma_pc"), &self.g_gamma_g1)?;
        let h_gamma_pc = P::prepare_g2(&mut cs.ns(|| "Prepare h_gamma_pc"), &self.h_gamma_g2)?;
        let h_pc = P::prepare_g2(&mut cs.ns(|| "Prepare h_pc"), &self.h_g2)?;
        Ok(PreparedVerifyingKeyGadget {
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
#[derivative(Clone(
    bound = "P::G1Gadget: Clone, P::GTGadget: Clone, P::G1PreparedGadget: Clone, \
             P::G2PreparedGadget: Clone, "
))]
pub struct PreparedVerifyingKeyGadget<
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
> {
    pub g_alpha: P::G1Gadget,
    pub h_beta: P::G2Gadget,
    pub g_alpha_pc: P::G1PreparedGadget,
    pub h_beta_pc: P::G2PreparedGadget,
    pub g_gamma_pc: P::G1PreparedGadget,
    pub h_gamma_pc: P::G2PreparedGadget,
    pub h_pc: P::G2PreparedGadget,
    pub query: Vec<P::G1Gadget>,
}

pub struct Gm17VerifierGadget<PairingE, ConstraintF, P>
where
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
{
    _pairing_engine: PhantomData<PairingE>,
    _engine: PhantomData<ConstraintF>,
    _pairing_gadget: PhantomData<P>,
}

impl<PairingE, ConstraintF, P, C, V> NIZKVerifierGadget<Gm17<PairingE, C, V>, ConstraintF>
    for Gm17VerifierGadget<PairingE, ConstraintF, P>
where
    PairingE: PairingEngine,
    ConstraintF: Field,
    C: ConstraintSynthesizer<PairingE::Fr>,
    V: ToConstraintField<PairingE::Fr>,
    P: PairingGadget<PairingE, ConstraintF>,
{
    type VerificationKeyGadget = VerifyingKeyGadget<PairingE, ConstraintF, P>;
    type ProofGadget = ProofGadget<PairingE, ConstraintF, P>;

    fn check_verify<'a, CS, I, T>(
        cs: CS,
        vk: &Self::VerificationKeyGadget,
        public_inputs: I,
        proof: &Self::ProofGadget,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<ConstraintF>,
        I: Iterator<Item = &'a T>,
        T: 'a + ToBitsGadget<ConstraintF> + ?Sized,
    {
        <Self as NIZKVerifierGadget<Gm17<PairingE, C, V>, ConstraintF>>::conditional_check_verify(
            cs,
            vk,
            public_inputs,
            proof,
            &Boolean::constant(true),
        )
    }

    fn conditional_check_verify<'a, CS, I, T>(
        mut cs: CS,
        vk: &Self::VerificationKeyGadget,
        mut public_inputs: I,
        proof: &Self::ProofGadget,
        condition: &Boolean,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<ConstraintF>,
        I: Iterator<Item = &'a T>,
        T: 'a + ToBitsGadget<ConstraintF> + ?Sized,
    {
        let pvk = vk.prepare(&mut cs.ns(|| "Prepare vk"))?;
        // e(A*G^{alpha}, B*H^{beta}) = e(G^{alpha}, H^{beta}) * e(G^{psi}, H^{gamma}) *
        // e(C, H) where psi = \sum_{i=0}^l input_i pvk.query[i]

        let g_psi = {
            let mut cs = cs.ns(|| "Process input");
            let mut g_psi = pvk.query[0].clone();
            let mut input_len = 1;
            for (i, (input, b)) in public_inputs
                .by_ref()
                .zip(pvk.query.iter().skip(1))
                .enumerate()
            {
                let input_bits = input.to_bits(cs.ns(|| format!("Input {}", i)))?;
                g_psi = b.mul_bits(cs.ns(|| format!("Mul {}", i)), &g_psi, input_bits.iter())?;
                input_len += 1;
            }
            // Check that the input and the query in the verification are of the
            // same length.
            assert!(input_len == pvk.query.len() && public_inputs.next().is_none());
            g_psi
        };

        let mut test1_a_g_alpha = proof.a.add(cs.ns(|| "A * G^{alpha}"), &pvk.g_alpha)?;
        let test1_b_h_beta = proof.b.add(cs.ns(|| "B * H^{beta}"), &pvk.h_beta)?;

        let test1_exp = {
            test1_a_g_alpha = test1_a_g_alpha.negate(cs.ns(|| "neg 1"))?;
            let test1_a_g_alpha_prep = P::prepare_g1(cs.ns(|| "First prep"), &test1_a_g_alpha)?;
            let test1_b_h_beta_prep = P::prepare_g2(cs.ns(|| "Second prep"), &test1_b_h_beta)?;

            let g_psi_prep = P::prepare_g1(cs.ns(|| "Third prep"), &g_psi)?;

            let c_prep = P::prepare_g1(cs.ns(|| "Fourth prep"), &proof.c)?;

            P::miller_loop(
                cs.ns(|| "Miller loop 1"),
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

        let test1 = P::final_exponentiation(cs.ns(|| "Final Exp 1"), &test1_exp).unwrap();

        // e(A, H^{gamma}) = e(G^{gamma}, B)
        let test2_exp = {
            let a_prep = P::prepare_g1(cs.ns(|| "Fifth prep"), &proof.a)?;
            // pvk.h_gamma_pc
            //&pvk.g_gamma_pc
            let proof_b = proof.b.negate(cs.ns(|| "Negate b"))?;
            let b_prep = P::prepare_g2(cs.ns(|| "Sixth prep"), &proof_b)?;
            P::miller_loop(
                cs.ns(|| "Miller loop 4"),
                &[a_prep, pvk.g_gamma_pc.clone()],
                &[pvk.h_gamma_pc, b_prep],
            )?
        };
        let test2 = P::final_exponentiation(cs.ns(|| "Final Exp 2"), &test2_exp)?;

        let one = P::GTGadget::one(cs.ns(|| "GT One"))?;
        test1.conditional_enforce_equal(cs.ns(|| "Test 1"), &one, condition)?;
        test2.conditional_enforce_equal(cs.ns(|| "Test 2"), &one, condition)?;
        Ok(())
    }
}

impl<PairingE, ConstraintF, P> AllocGadget<VerifyingKey<PairingE>, ConstraintF>
    for VerifyingKeyGadget<PairingE, ConstraintF, P>
where
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
{
    #[inline]
    fn alloc_constant<T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        val: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<VerifyingKey<PairingE>>,
    {
        let VerifyingKey {
            h_g2,
            g_alpha_g1,
            h_beta_g2,
            g_gamma_g1,
            h_gamma_g2,
            query,
        } = val.borrow().clone();
        let h_g2 = P::G2Gadget::alloc_constant(cs.ns(|| "h_g2"), h_g2.into_projective())?;
        let g_alpha_g1 =
            P::G1Gadget::alloc_constant(cs.ns(|| "g_alpha"), g_alpha_g1.into_projective())?;
        let h_beta_g2 =
            P::G2Gadget::alloc_constant(cs.ns(|| "h_beta"), h_beta_g2.into_projective())?;
        let g_gamma_g1 =
            P::G1Gadget::alloc_constant(cs.ns(|| "g_gamma_g1"), g_gamma_g1.into_projective())?;
        let h_gamma_g2 =
            P::G2Gadget::alloc_constant(cs.ns(|| "h_gamma_g2"), h_gamma_g2.into_projective())?;

        let query = query
            .into_iter()
            .enumerate()
            .map(|(i, query_i)| {
                P::G1Gadget::alloc_constant(
                    cs.ns(|| format!("query_{}", i)),
                    query_i.into_projective(),
                )
            })
            .collect::<Vec<_>>()
            .into_iter()
            .collect::<Result<_, _>>()?;
        Ok(Self {
            h_g2,
            g_alpha_g1,
            h_beta_g2,
            g_gamma_g1,
            h_gamma_g2,
            query,
        })
    }

    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<VerifyingKey<PairingE>>,
    {
        value_gen().and_then(|vk| {
            let VerifyingKey {
                h_g2,
                g_alpha_g1,
                h_beta_g2,
                g_gamma_g1,
                h_gamma_g2,
                query,
            } = vk.borrow().clone();
            let h_g2 = P::G2Gadget::alloc(cs.ns(|| "h_g2"), || Ok(h_g2.into_projective()))?;
            let g_alpha_g1 =
                P::G1Gadget::alloc(cs.ns(|| "g_alpha"), || Ok(g_alpha_g1.into_projective()))?;
            let h_beta_g2 =
                P::G2Gadget::alloc(cs.ns(|| "h_beta"), || Ok(h_beta_g2.into_projective()))?;
            let g_gamma_g1 =
                P::G1Gadget::alloc(cs.ns(|| "g_gamma_g1"), || Ok(g_gamma_g1.into_projective()))?;
            let h_gamma_g2 =
                P::G2Gadget::alloc(cs.ns(|| "h_gamma_g2"), || Ok(h_gamma_g2.into_projective()))?;

            let query = query
                .into_iter()
                .enumerate()
                .map(|(i, query_i)| {
                    P::G1Gadget::alloc(cs.ns(|| format!("query_{}", i)), || {
                        Ok(query_i.into_projective())
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .collect::<Result<_, _>>()?;
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

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<VerifyingKey<PairingE>>,
    {
        value_gen().and_then(|vk| {
            let VerifyingKey {
                h_g2,
                g_alpha_g1,
                h_beta_g2,
                g_gamma_g1,
                h_gamma_g2,
                query,
            } = vk.borrow().clone();
            let h_g2 = P::G2Gadget::alloc_input(cs.ns(|| "h_g2"), || Ok(h_g2.into_projective()))?;
            let g_alpha_g1 =
                P::G1Gadget::alloc_input(cs.ns(|| "g_alpha"), || Ok(g_alpha_g1.into_projective()))?;
            let h_beta_g2 =
                P::G2Gadget::alloc_input(cs.ns(|| "h_beta"), || Ok(h_beta_g2.into_projective()))?;
            let g_gamma_g1 = P::G1Gadget::alloc_input(cs.ns(|| "g_gamma_g1"), || {
                Ok(g_gamma_g1.into_projective())
            })?;
            let h_gamma_g2 = P::G2Gadget::alloc_input(cs.ns(|| "h_gamma_g2"), || {
                Ok(h_gamma_g2.into_projective())
            })?;

            let query = query
                .into_iter()
                .enumerate()
                .map(|(i, query_i)| {
                    P::G1Gadget::alloc_input(cs.ns(|| format!("query_{}", i)), || {
                        Ok(query_i.into_projective())
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .collect::<Result<_, _>>()?;
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

impl<PairingE, ConstraintF, P> AllocGadget<Proof<PairingE>, ConstraintF>
    for ProofGadget<PairingE, ConstraintF, P>
where
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
{
    #[inline]
    fn alloc_constant<T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        val: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<Proof<PairingE>>,
    {
        let Proof { a, b, c } = val.borrow().clone();
        let a = P::G1Gadget::alloc_constant(cs.ns(|| "a"), a.into_projective())?;
        let b = P::G2Gadget::alloc_constant(cs.ns(|| "b"), b.into_projective())?;
        let c = P::G1Gadget::alloc_constant(cs.ns(|| "c"), c.into_projective())?;
        Ok(Self { a, b, c })
    }

    #[inline]
    fn alloc<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Proof<PairingE>>,
    {
        value_gen().and_then(|proof| {
            let Proof { a, b, c } = proof.borrow().clone();
            let a = P::G1Gadget::alloc_checked(cs.ns(|| "a"), || Ok(a.into_projective()))?;
            let b = P::G2Gadget::alloc_checked(cs.ns(|| "b"), || Ok(b.into_projective()))?;
            let c = P::G1Gadget::alloc_checked(cs.ns(|| "c"), || Ok(c.into_projective()))?;
            Ok(Self { a, b, c })
        })
    }

    #[inline]
    fn alloc_input<FN, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: FN,
    ) -> Result<Self, SynthesisError>
    where
        FN: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<Proof<PairingE>>,
    {
        value_gen().and_then(|proof| {
            let Proof { a, b, c } = proof.borrow().clone();
            // We don't need to check here because the prime order check can be performed
            // in plain.
            let a = P::G1Gadget::alloc_input(cs.ns(|| "a"), || Ok(a.into_projective()))?;
            let b = P::G2Gadget::alloc_input(cs.ns(|| "b"), || Ok(b.into_projective()))?;
            let c = P::G1Gadget::alloc_input(cs.ns(|| "c"), || Ok(c.into_projective()))?;
            Ok(Self { a, b, c })
        })
    }
}

impl<PairingE, ConstraintF, P> ToBytesGadget<ConstraintF>
    for VerifyingKeyGadget<PairingE, ConstraintF, P>
where
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
{
    #[inline]
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&self.h_g2.to_bytes(&mut cs.ns(|| "h_g2 to bytes"))?);
        bytes.extend_from_slice(
            &self
                .g_alpha_g1
                .to_bytes(&mut cs.ns(|| "g_alpha_g1 to bytes"))?,
        );
        bytes.extend_from_slice(
            &self
                .h_beta_g2
                .to_bytes(&mut cs.ns(|| "h_beta_g2 to bytes"))?,
        );
        bytes.extend_from_slice(
            &self
                .g_gamma_g1
                .to_bytes(&mut cs.ns(|| "g_gamma_g1 to bytes"))?,
        );
        bytes.extend_from_slice(
            &self
                .h_gamma_g2
                .to_bytes(&mut cs.ns(|| "h_gamma_g2 to bytes"))?,
        );
        for (i, q) in self.query.iter().enumerate() {
            let mut cs = cs.ns(|| format!("Iteration {}", i));
            bytes.extend_from_slice(&q.to_bytes(&mut cs.ns(|| "q"))?);
        }
        Ok(bytes)
    }
}

#[cfg(test)]
mod test {
    use gm17::*;
    use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};

    use super::*;
    use algebra::{
        bls12_377::{Bls12_377, Fq, Fr},
        test_rng, BitIterator, PrimeField,
    };
    use r1cs_std::{
        bls12_377::PairingGadget as Bls12_377PairingGadget, boolean::Boolean,
        test_constraint_system::TestConstraintSystem,
    };
    use rand::Rng;

    type TestProofSystem = Gm17<Bls12_377, Bench<Fr>, Fr>;
    type TestVerifierGadget = Gm17VerifierGadget<Bls12_377, Fq, Bls12_377PairingGadget>;
    type TestProofGadget = ProofGadget<Bls12_377, Fq, Bls12_377PairingGadget>;
    type TestVkGadget = VerifyingKeyGadget<Bls12_377, Fq, Bls12_377PairingGadget>;

    struct Bench<F: Field> {
        inputs: Vec<Option<F>>,
        num_constraints: usize,
    }

    impl<F: Field> ConstraintSynthesizer<F> for Bench<F> {
        fn generate_constraints<CS: ConstraintSystem<F>>(
            self,
            cs: &mut CS,
        ) -> Result<(), SynthesisError> {
            assert!(self.inputs.len() >= 2);
            assert!(self.num_constraints >= self.inputs.len());

            let mut variables: Vec<_> = Vec::with_capacity(self.inputs.len());
            for (i, input) in self.inputs.into_iter().enumerate() {
                let input_var = cs.alloc_input(
                    || format!("Input {}", i),
                    || input.ok_or(SynthesisError::AssignmentMissing),
                )?;
                variables.push((input, input_var));
            }

            for i in 0..self.num_constraints {
                let new_entry = {
                    let (input_1_val, input_1_var) = variables[i];
                    let (input_2_val, input_2_var) = variables[i + 1];
                    let result_val = input_1_val
                        .and_then(|input_1| input_2_val.map(|input_2| input_1 * &input_2));
                    let result_var = cs.alloc(
                        || format!("Result {}", i),
                        || result_val.ok_or(SynthesisError::AssignmentMissing),
                    )?;
                    cs.enforce(
                        || format!("Enforce constraint {}", i),
                        |lc| lc + input_1_var,
                        |lc| lc + input_2_var,
                        |lc| lc + result_var,
                    );
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
            let mut cs = TestConstraintSystem::<Fq>::new();

            let inputs: Vec<_> = inputs.into_iter().map(|input| input.unwrap()).collect();
            let mut input_gadgets = Vec::new();

            {
                let mut cs = cs.ns(|| "Allocate Input");
                for (i, input) in inputs.into_iter().enumerate() {
                    let mut input_bits = BitIterator::new(input.into_repr()).collect::<Vec<_>>();
                    // Input must be in little-endian, but BitIterator outputs in big-endian.
                    input_bits.reverse();

                    let input_bits =
                        Vec::<Boolean>::alloc_input(cs.ns(|| format!("Input {}", i)), || {
                            Ok(input_bits)
                        })
                        .unwrap();
                    input_gadgets.push(input_bits);
                }
            }

            let vk_gadget = TestVkGadget::alloc_input(cs.ns(|| "Vk"), || Ok(&params.vk)).unwrap();
            let proof_gadget =
                TestProofGadget::alloc(cs.ns(|| "Proof"), || Ok(proof.clone())).unwrap();
            println!("Time to verify!\n\n\n\n");
            <TestVerifierGadget as NIZKVerifierGadget<TestProofSystem, Fq>>::check_verify(
                cs.ns(|| "Verify"),
                &vk_gadget,
                input_gadgets.iter(),
                &proof_gadget,
            )
            .unwrap();
            if !cs.is_satisfied() {
                println!("=========================================================");
                println!("Unsatisfied constraints:");
                println!("{:?}", cs.which_is_unsatisfied().unwrap());
                println!("=========================================================");
            }

            // cs.print_named_objects();
            assert!(cs.is_satisfied());
        }
    }
}

#[cfg(test)]
mod test_recursive {
    use gm17::*;
    use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};

    use super::*;
    use algebra::{
        fields::{FftParameters, FpParameters},
        mnt4_298::{Fq as MNT4Fq, FqParameters as MNT4FqParameters, Fr as MNT4Fr, MNT4_298},
        mnt6_298::{Fq as MNT6Fq, FqParameters as MNT6FqParameters, Fr as MNT6Fr, MNT6_298},
        test_rng, BigInteger, PrimeField,
    };
    use r1cs_std::{
        fields::fp::FpGadget, mnt4_298::PairingGadget as MNT4_298PairingGadget,
        mnt6_298::PairingGadget as MNT6_298PairingGadget,
        test_constraint_system::TestConstraintSystem, uint8::UInt8,
    };
    use rand::Rng;

    type TestProofSystem1 = Gm17<MNT6_298, Bench<MNT4Fq>, MNT6Fr>;
    type TestVerifierGadget1 = Gm17VerifierGadget<MNT6_298, MNT6Fq, MNT6_298PairingGadget>;
    type TestProofGadget1 = ProofGadget<MNT6_298, MNT6Fq, MNT6_298PairingGadget>;
    type TestVkGadget1 = VerifyingKeyGadget<MNT6_298, MNT6Fq, MNT6_298PairingGadget>;

    type TestProofSystem2 = Gm17<MNT4_298, Wrapper, MNT4Fr>;
    type TestVerifierGadget2 = Gm17VerifierGadget<MNT4_298, MNT4Fq, MNT4_298PairingGadget>;
    type TestProofGadget2 = ProofGadget<MNT4_298, MNT4Fq, MNT4_298PairingGadget>;
    type TestVkGadget2 = VerifyingKeyGadget<MNT4_298, MNT4Fq, MNT4_298PairingGadget>;

    #[derive(Clone)]
    struct Bench<F: Field> {
        inputs: Vec<Option<F>>,
        num_constraints: usize,
    }

    impl<F: Field> ConstraintSynthesizer<F> for Bench<F> {
        fn generate_constraints<CS: ConstraintSystem<F>>(
            self,
            cs: &mut CS,
        ) -> Result<(), SynthesisError> {
            assert!(self.inputs.len() >= 2);
            assert!(self.num_constraints >= self.inputs.len());

            let mut variables: Vec<_> = Vec::with_capacity(self.inputs.len());
            for (i, input) in self.inputs.into_iter().enumerate() {
                let input_var = cs.alloc_input(
                    || format!("Input {}", i),
                    || input.ok_or(SynthesisError::AssignmentMissing),
                )?;
                variables.push((input, input_var));
            }

            for i in 0..self.num_constraints {
                let new_entry = {
                    let (input_1_val, input_1_var) = variables[i];
                    let (input_2_val, input_2_var) = variables[i + 1];
                    let result_val = input_1_val
                        .and_then(|input_1| input_2_val.map(|input_2| input_1 * &input_2));
                    let result_var = cs.alloc(
                        || format!("Result {}", i),
                        || result_val.ok_or(SynthesisError::AssignmentMissing),
                    )?;
                    cs.enforce(
                        || format!("Enforce constraint {}", i),
                        |lc| lc + input_1_var,
                        |lc| lc + input_2_var,
                        |lc| lc + result_var,
                    );
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
        fn generate_constraints<CS: ConstraintSystem<MNT6Fq>>(
            self,
            cs: &mut CS,
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
                let mut cs = cs.ns(|| "Allocate Input");
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
                let input_bytes = UInt8::alloc_input_vec(cs.ns(|| "Input"), &input_bytes[..])?;
                // 40 byte
                let element_size = <MNT4FqParameters as FftParameters>::BigInt::NUM_LIMBS * 8;
                input_gadgets = input_bytes
                    .chunks(element_size)
                    .map(|chunk| {
                        chunk
                            .iter()
                            .flat_map(|byte| byte.into_bits_le())
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>();
            }

            let vk_gadget = TestVkGadget1::alloc(cs.ns(|| "Vk"), || Ok(&params.vk))?;
            let proof_gadget =
                TestProofGadget1::alloc(cs.ns(|| "Proof"), || Ok(proof.clone())).unwrap();
            <TestVerifierGadget1 as NIZKVerifierGadget<TestProofSystem1, MNT6Fq>>::check_verify(
                cs.ns(|| "Verify"),
                &vk_gadget,
                input_gadgets.iter(),
                &proof_gadget,
            )?;
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
            // Create a groth16 proof with our parameters.
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
                // Create a groth16 proof with our parameters.
                create_random_proof(c, &params, rng).unwrap()
            };

            let mut cs = TestConstraintSystem::<MNT4Fq>::new();

            let inputs: Vec<_> = inputs.into_iter().map(|input| input.unwrap()).collect();
            let mut input_gadgets = Vec::new();

            {
                let bigint_size = <MNT4FqParameters as FftParameters>::BigInt::NUM_LIMBS * 64;
                let mut input_bits = Vec::new();
                let mut cs = cs.ns(|| "Allocate Input");
                for (i, input) in inputs.into_iter().enumerate() {
                    let input_gadget =
                        FpGadget::alloc_input(cs.ns(|| format!("Input {}", i)), || Ok(input))
                            .unwrap();
                    let mut fp_bits = input_gadget
                        .to_bits(cs.ns(|| format!("To bits {}", i)))
                        .unwrap();

                    // FpGadget::to_bits outputs a big-endian binary representation of
                    // fe_gadget's value, so we have to reverse it to get the little-endian
                    // form.
                    fp_bits.reverse();

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

            let vk_gadget = TestVkGadget2::alloc_input(cs.ns(|| "Vk"), || Ok(&params.vk)).unwrap();
            let proof_gadget =
                TestProofGadget2::alloc(cs.ns(|| "Proof"), || Ok(proof.clone())).unwrap();
            println!("Time to verify!\n\n\n\n");
            <TestVerifierGadget2 as NIZKVerifierGadget<TestProofSystem2, MNT4Fq>>::check_verify(
                cs.ns(|| "Verify"),
                &vk_gadget,
                input_gadgets.iter(),
                &proof_gadget,
            )
            .unwrap();
            if !cs.is_satisfied() {
                println!("=========================================================");
                println!("Unsatisfied constraints:");
                println!("{:?}", cs.which_is_unsatisfied().unwrap());
                println!("=========================================================");
            }

            // cs.print_named_objects();
            assert!(cs.is_satisfied());
        }
    }
}
