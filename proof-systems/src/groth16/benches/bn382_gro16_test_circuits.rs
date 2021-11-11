use algebra::{curves::bn_382::Bn382, fields::bn_382::Fr, PrimeField, UniformRand};
use proof_systems::groth16::{create_random_proof, generate_random_parameters};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, SynthesisError};

use criterion::Criterion;
use criterion::{BatchSize, BenchmarkId};
use r1cs_std::alloc::AllocGadget;
use r1cs_std::eq::EqGadget;
use r1cs_std::fields::fp::FpGadget;
use r1cs_std::fields::FieldGadget;
use r1cs_std::Assignment;

use rand::{rngs::OsRng, thread_rng};

use std::time::{SystemTime, UNIX_EPOCH};

#[macro_use]
extern crate criterion;

#[macro_use]
extern crate bench_utils;

/// Has R1CS density d = 1 and keeps the synthesizer costs low (no field inversion needed).
pub struct TestCircuit1<F: PrimeField> {
    num_constraints: usize,
    a: Option<F>,
    b: Option<F>,
}

/// Hash R1CS density d = 1 and that every second constraint costs the synthesizer a field
/// inversion (this should simulate EC arithmetics over prime fields).
pub struct TestCircuit2<F: PrimeField> {
    num_constraints: usize,
    a: Option<F>,
    b: Option<F>,
}

impl<F: PrimeField> ConstraintSynthesizer<F> for TestCircuit1<F> {
    fn generate_constraints<CS: ConstraintSystem<F>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
        let mut a_k_minus_1 = FpGadget::<F>::alloc_input(cs.ns(|| "alloc a"), || {
            self.a.ok_or(SynthesisError::AssignmentMissing)
        })?;

        let mut b_k_minus_1 = FpGadget::<F>::alloc_input(cs.ns(|| "alloc b"), || {
            self.b.ok_or(SynthesisError::AssignmentMissing)
        })?;

        let zero = FpGadget::<F>::zero(cs.ns(|| "alloc zero"))?;

        a_k_minus_1.enforce_not_equal(cs.ns(|| "a_0 != 0"), &zero)?;
        b_k_minus_1.enforce_not_equal(cs.ns(|| "b_0 != 0"), &zero)?;

        for k in 0..(self.num_constraints - 5) / 2 {
            let a_k = FpGadget::<F>::alloc(cs.ns(|| format!("alloc a_{}", k)), || {
                Ok(a_k_minus_1.value.get()? * &b_k_minus_1.value.get()?)
            })?;

            let b_k = FpGadget::<F>::alloc(cs.ns(|| format!("alloc b_{}", k)), || {
                Ok(b_k_minus_1.value.get()? * &a_k_minus_1.value.get()?)
            })?;

            a_k_minus_1.mul_equals(
                cs.ns(|| format!("a_{} * b_{} = a_{}", k - 1, k - 1, k)),
                &b_k_minus_1,
                &a_k,
            )?;

            b_k_minus_1.mul_equals(
                cs.ns(|| format!("b_{} * a_{} = b_{}", k - 1, k - 1, k)),
                &a_k_minus_1,
                &b_k,
            )?;

            a_k_minus_1 = a_k;
            b_k_minus_1 = b_k;
        }

        Ok(())
    }
}

impl<F: PrimeField> ConstraintSynthesizer<F> for TestCircuit2<F> {
    fn generate_constraints<CS: ConstraintSystem<F>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
        let mut a_k_minus_1 = FpGadget::<F>::alloc_input(cs.ns(|| "alloc a"), || {
            self.a.ok_or(SynthesisError::AssignmentMissing)
        })?;

        let mut b_k_minus_1 = FpGadget::<F>::alloc_input(cs.ns(|| "alloc b"), || {
            self.b.ok_or(SynthesisError::AssignmentMissing)
        })?;

        let zero = FpGadget::<F>::zero(cs.ns(|| "alloc zero"))?;

        a_k_minus_1.enforce_not_equal(cs.ns(|| "a_0 != 0"), &zero)?;
        b_k_minus_1.enforce_not_equal(cs.ns(|| "b_0 != 0"), &zero)?;

        for k in 0..(self.num_constraints - 5) / 2 {
            let a_k = FpGadget::<F>::alloc(cs.ns(|| format!("alloc a_{}", k)), || {
                Ok(a_k_minus_1.value.get()? * &b_k_minus_1.value.get()?.inverse().get()?)
            })?;

            let b_k = FpGadget::<F>::alloc(cs.ns(|| format!("alloc b_{}", k)), || {
                Ok(b_k_minus_1.value.get()? * &a_k_minus_1.value.get()?)
            })?;

            a_k.mul_equals(
                cs.ns(|| format!("a_{} * b_{} = a_{}", k, k - 1, k - 1)),
                &b_k_minus_1,
                &a_k_minus_1,
            )?;

            b_k_minus_1.mul_equals(
                cs.ns(|| format!("b_{} * a_{} = b_{}", k - 1, k - 1, k)),
                &a_k_minus_1,
                &b_k,
            )?;

            a_k_minus_1 = a_k;
            b_k_minus_1 = b_k;
        }

        Ok(())
    }
}

fn bench_prover_circuit1(c: &mut Criterion) {
    let mut rng = thread_rng();

    let mut group = c.benchmark_group("gro16-bn382-test circuit 1-variable constraints");

    let num_constraints = (14..=22).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    for &num_constraints in num_constraints.iter() {
        let params = {
            let c = TestCircuit1::<Fr> {
                num_constraints,
                a: None,
                b: None,
            };
            generate_random_parameters::<Bn382, _, _>(c, &mut rng).unwrap()
        };

        add_to_trace!(
            || format!("****************{}*******************", num_constraints),
            || format!(
                "--->START TIMESTAMP: {:?}",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_secs()
            )
        );

        group.bench_with_input(
            BenchmarkId::from_parameter(num_constraints),
            &num_constraints,
            |bn, _constraints| {
                bn.iter_batched(
                    || {
                        let mut rng = OsRng::default();
                        let a = Fr::rand(&mut rng);
                        let b = Fr::rand(&mut rng);
                        (a, b)
                    },
                    |(a, b)| {
                        let c = TestCircuit1 {
                            num_constraints,
                            a: Some(a),
                            b: Some(b),
                        };

                        create_random_proof(c, &params, &mut rng).unwrap();
                    },
                    BatchSize::PerIteration,
                );
            },
        );
        add_to_trace!(
            || format!("****************{}*******************", num_constraints),
            || format!(
                "--->END TIMESTAMP: {:?}",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_secs()
            )
        );
    }
    group.finish();
}

fn bench_prover_circuit2(c: &mut Criterion) {
    let mut rng = thread_rng();

    let mut group = c.benchmark_group("gro16-bn382-test circuit 2-variable constraints");

    let num_constraints = (14..=22).map(|i| 2usize.pow(i)).collect::<Vec<_>>();

    for &num_constraints in num_constraints.iter() {
        let params = {
            let c = TestCircuit2::<Fr> {
                num_constraints,
                a: None,
                b: None,
            };
            generate_random_parameters::<Bn382, _, _>(c, &mut rng).unwrap()
        };

        add_to_trace!(
            || format!("****************{}*******************", num_constraints),
            || format!(
                "--->START TIMESTAMP: {:?}",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_secs()
            )
        );
        group.bench_with_input(
            BenchmarkId::from_parameter(num_constraints),
            &num_constraints,
            |bn, _constraints| {
                bn.iter_batched(
                    || {
                        let mut rng = OsRng::default();
                        let a = Fr::rand(&mut rng);
                        let b = Fr::rand(&mut rng);
                        (a, b)
                    },
                    |(a, b)| {
                        let c = TestCircuit2 {
                            num_constraints,
                            a: Some(a),
                            b: Some(b),
                        };

                        create_random_proof(c, &params, &mut rng).unwrap();
                    },
                    BatchSize::PerIteration,
                );
            },
        );
        add_to_trace!(
            || format!("****************{}*******************", num_constraints),
            || format!(
                "--->END TIMESTAMP: {:?}",
                SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_secs()
            )
        );
    }
    group.finish();
}

criterion_group!(
name = bn382_gro16_test_circuits;
config = Criterion::default().sample_size(10);
targets = bench_prover_circuit1, bench_prover_circuit2
);

criterion_main!(bn382_gro16_test_circuits);
