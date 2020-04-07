use crate::nizk::{groth16::Groth16, NIZKVerifierGadget};
use algebra_core::{AffineCurve, Field, PairingEngine, ToConstraintField};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

use core::{borrow::Borrow, marker::PhantomData};
use groth16::{Proof, VerifyingKey};

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
    pub alpha_g1: P::G1Gadget,
    pub beta_g2: P::G2Gadget,
    pub gamma_g2: P::G2Gadget,
    pub delta_g2: P::G2Gadget,
    pub gamma_abc_g1: Vec<P::G1Gadget>,
}

impl<PairingE: PairingEngine, ConstraintF: Field, P: PairingGadget<PairingE, ConstraintF>>
    VerifyingKeyGadget<PairingE, ConstraintF, P>
{
    pub fn prepare<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<PreparedVerifyingKeyGadget<PairingE, ConstraintF, P>, SynthesisError> {
        let mut cs = cs.ns(|| "Preparing verifying key");
        let alpha_g1_pc = P::prepare_g1(&mut cs.ns(|| "Prepare alpha_g1"), &self.alpha_g1)?;
        let beta_g2_pc = P::prepare_g2(&mut cs.ns(|| "Prepare beta_g2"), &self.beta_g2)?;

        let alpha_g1_beta_g2 = P::pairing(
            &mut cs.ns(|| "Precompute e(alpha_g1, beta_g2)"),
            alpha_g1_pc,
            beta_g2_pc,
        )?;

        let gamma_g2_neg = self.gamma_g2.negate(&mut cs.ns(|| "Negate gamma_g2"))?;
        let gamma_g2_neg_pc = P::prepare_g2(&mut cs.ns(|| "Prepare gamma_g2_neg"), &gamma_g2_neg)?;

        let delta_g2_neg = self.delta_g2.negate(&mut cs.ns(|| "Negate delta_g2"))?;
        let delta_g2_neg_pc = P::prepare_g2(&mut cs.ns(|| "Prepare delta_g2_neg"), &delta_g2_neg)?;

        Ok(PreparedVerifyingKeyGadget {
            alpha_g1_beta_g2,
            gamma_g2_neg_pc,
            delta_g2_neg_pc,
            gamma_abc_g1: self.gamma_abc_g1.clone(),
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
    pub alpha_g1_beta_g2: P::GTGadget,
    pub gamma_g2_neg_pc: P::G2PreparedGadget,
    pub delta_g2_neg_pc: P::G2PreparedGadget,
    pub gamma_abc_g1: Vec<P::G1Gadget>,
}

pub struct Groth16VerifierGadget<PairingE, ConstraintF, P>
where
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
{
    _pairing_engine: PhantomData<PairingE>,
    _engine: PhantomData<ConstraintF>,
    _pairing_gadget: PhantomData<P>,
}

impl<PairingE, ConstraintF, P, C, V> NIZKVerifierGadget<Groth16<PairingE, C, V>, ConstraintF>
    for Groth16VerifierGadget<PairingE, ConstraintF, P>
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
        mut cs: CS,
        vk: &Self::VerificationKeyGadget,
        mut public_inputs: I,
        proof: &Self::ProofGadget,
    ) -> Result<(), SynthesisError>
    where
        CS: ConstraintSystem<ConstraintF>,
        I: Iterator<Item = &'a T>,
        T: 'a + ToBitsGadget<ConstraintF> + ?Sized,
    {
        let pvk = vk.prepare(&mut cs.ns(|| "Prepare vk"))?;

        let g_ic = {
            let mut cs = cs.ns(|| "Process input");
            let mut g_ic = pvk.gamma_abc_g1[0].clone();
            let mut input_len = 1;
            for (i, (input, b)) in public_inputs
                .by_ref()
                .zip(pvk.gamma_abc_g1.iter().skip(1))
                .enumerate()
            {
                let input_bits = input.to_bits(cs.ns(|| format!("Input {}", i)))?;
                g_ic = b.mul_bits(cs.ns(|| format!("Mul {}", i)), &g_ic, input_bits.iter())?;
                input_len += 1;
            }
            // Check that the input and the query in the verification are of the
            // same length.
            assert!(input_len == pvk.gamma_abc_g1.len() && public_inputs.next().is_none());
            g_ic
        };

        let test_exp = {
            let proof_a_prep = P::prepare_g1(cs.ns(|| "Prepare proof a"), &proof.a)?;
            let proof_b_prep = P::prepare_g2(cs.ns(|| "Prepare proof b"), &proof.b)?;
            let proof_c_prep = P::prepare_g1(cs.ns(|| "Prepare proof c"), &proof.c)?;

            let g_ic_prep = P::prepare_g1(cs.ns(|| "Prepare g_ic"), &g_ic)?;

            P::miller_loop(
                cs.ns(|| "Miller loop 1"),
                &[proof_a_prep, g_ic_prep, proof_c_prep],
                &[
                    proof_b_prep,
                    pvk.gamma_g2_neg_pc.clone(),
                    pvk.delta_g2_neg_pc.clone(),
                ],
            )?
        };

        let test = P::final_exponentiation(cs.ns(|| "Final Exp"), &test_exp).unwrap();

        test.enforce_equal(cs.ns(|| "Test 1"), &pvk.alpha_g1_beta_g2)?;
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
                alpha_g1,
                beta_g2,
                gamma_g2,
                delta_g2,
                gamma_abc_g1,
            } = vk.borrow().clone();
            let alpha_g1 =
                P::G1Gadget::alloc(cs.ns(|| "alpha_g1"), || Ok(alpha_g1.into_projective()))?;
            let beta_g2 =
                P::G2Gadget::alloc(cs.ns(|| "beta_g2"), || Ok(beta_g2.into_projective()))?;
            let gamma_g2 =
                P::G2Gadget::alloc(cs.ns(|| "gamma_g2"), || Ok(gamma_g2.into_projective()))?;
            let delta_g2 =
                P::G2Gadget::alloc(cs.ns(|| "delta_g2"), || Ok(delta_g2.into_projective()))?;

            let gamma_abc_g1 = gamma_abc_g1
                .into_iter()
                .enumerate()
                .map(|(i, gamma_abc_i)| {
                    P::G1Gadget::alloc(cs.ns(|| format!("gamma_abc_{}", i)), || {
                        Ok(gamma_abc_i.into_projective())
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .collect::<Result<_, _>>()?;
            Ok(Self {
                alpha_g1,
                beta_g2,
                gamma_g2,
                delta_g2,
                gamma_abc_g1,
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
                alpha_g1,
                beta_g2,
                gamma_g2,
                delta_g2,
                gamma_abc_g1,
            } = vk.borrow().clone();
            let alpha_g1 =
                P::G1Gadget::alloc_input(cs.ns(|| "alpha_g1"), || Ok(alpha_g1.into_projective()))?;
            let beta_g2 =
                P::G2Gadget::alloc_input(cs.ns(|| "beta_g2"), || Ok(beta_g2.into_projective()))?;
            let gamma_g2 =
                P::G2Gadget::alloc_input(cs.ns(|| "gamma_g2"), || Ok(gamma_g2.into_projective()))?;
            let delta_g2 =
                P::G2Gadget::alloc_input(cs.ns(|| "delta_g2"), || Ok(delta_g2.into_projective()))?;

            let gamma_abc_g1 = gamma_abc_g1
                .into_iter()
                .enumerate()
                .map(|(i, gamma_abc_i)| {
                    P::G1Gadget::alloc_input(cs.ns(|| format!("gamma_abc_{}", i)), || {
                        Ok(gamma_abc_i.into_projective())
                    })
                })
                .collect::<Vec<_>>()
                .into_iter()
                .collect::<Result<_, _>>()?;

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

impl<PairingE, ConstraintF, P> AllocGadget<Proof<PairingE>, ConstraintF>
    for ProofGadget<PairingE, ConstraintF, P>
where
    PairingE: PairingEngine,
    ConstraintF: Field,
    P: PairingGadget<PairingE, ConstraintF>,
{
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
        bytes.extend_from_slice(&self.alpha_g1.to_bytes(&mut cs.ns(|| "alpha_g1 to bytes"))?);
        bytes.extend_from_slice(&self.beta_g2.to_bytes(&mut cs.ns(|| "beta_g2 to bytes"))?);
        bytes.extend_from_slice(&self.gamma_g2.to_bytes(&mut cs.ns(|| "gamma_g2 to bytes"))?);
        bytes.extend_from_slice(&self.delta_g2.to_bytes(&mut cs.ns(|| "delta_g2 to bytes"))?);
        for (i, g) in self.gamma_abc_g1.iter().enumerate() {
            let mut cs = cs.ns(|| format!("Iteration {}", i));
            bytes.extend_from_slice(&g.to_bytes(&mut cs.ns(|| "g"))?);
        }
        Ok(bytes)
    }
}

#[cfg(test)]
mod test {
    use groth16::*;
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

    type TestProofSystem = Groth16<Bls12_377, Bench<Fr>, Fr>;
    type TestVerifierGadget = Groth16VerifierGadget<Bls12_377, Fq, Bls12_377PairingGadget>;
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
    fn groth16_verifier_test() {
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
                // Create a groth16 proof with our parameters.
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
    use groth16::*;
    use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};

    use super::*;
    use algebra::{
        fields::FpParameters,
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

    type TestProofSystem1 = Groth16<MNT6_298, Bench<MNT4Fq>, MNT6Fr>;
    type TestVerifierGadget1 = Groth16VerifierGadget<MNT6_298, MNT6Fq, MNT6_298PairingGadget>;
    type TestProofGadget1 = ProofGadget<MNT6_298, MNT6Fq, MNT6_298PairingGadget>;
    type TestVkGadget1 = VerifyingKeyGadget<MNT6_298, MNT6Fq, MNT6_298PairingGadget>;

    type TestProofSystem2 = Groth16<MNT4_298, Wrapper, MNT4Fr>;
    type TestVerifierGadget2 = Groth16VerifierGadget<MNT4_298, MNT4Fq, MNT4_298PairingGadget>;
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
                let element_size = <MNT4FqParameters as FpParameters>::BigInt::NUM_LIMBS * 8;
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
    fn groth16_recursive_verifier_test() {
        let num_inputs = 100;
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
                let bigint_size = <MNT4FqParameters as FpParameters>::BigInt::NUM_LIMBS * 64;
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
                let bigint_size = <MNT6FqParameters as FpParameters>::BigInt::NUM_LIMBS * 64;
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
