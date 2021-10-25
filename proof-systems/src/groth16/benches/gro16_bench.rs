#[macro_use]
extern crate criterion;

use algebra::{curves::bn_382::Bn382, fields::bn_382::Fr, Field, PrimeField, UniformRand};
use proof_systems::groth16::{create_random_proof, generate_random_parameters};
use r1cs_core::{ConstraintSynthesizer, ConstraintSystem, LinearCombination, SynthesisError};

use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use criterion::Criterion;
use criterion::{BatchSize, BenchmarkId};

use std::marker::PhantomData;

/// Circuit designed to have low R1CS density (d = 1) and synthesization cost (no field inversions).
pub struct Benchmark<F: PrimeField> {
    inputs: Vec<Option<F>>,
    num_constraints: usize,
}

impl<F: PrimeField> ConstraintSynthesizer<F> for Benchmark<F> {
    fn generate_constraints<CS: ConstraintSystem<F>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
        assert!(self.inputs.len() >= 2);
        assert!(self.num_constraints >= self.inputs.len());

        let mut variables: Vec<_> = Vec::with_capacity(self.inputs.len() + self.num_constraints);
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
                let result_val =
                    input_1_val.and_then(|input_1| input_2_val.map(|input_2| input_1 * &input_2));
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

pub struct BenchmarkHighDensities<F: Field> {
    num_constraints: usize,
    _engine: PhantomData<F>,
}

impl<F: Field> ConstraintSynthesizer<F> for BenchmarkHighDensities<F> {
    fn generate_constraints<CS: ConstraintSystem<F>>(
        self,
        cs: &mut CS,
    ) -> Result<(), SynthesisError> {
        let mut assignments = Vec::new();

        let mut a_val = F::one();
        let mut a_var = cs.alloc_input(|| "a", || Ok(a_val))?;
        assignments.push((a_val, a_var));

        let mut b_val = F::one();
        let mut b_var = cs.alloc_input(|| "b", || Ok(b_val))?;
        assignments.push((a_val, a_var));

        for i in 0..self.num_constraints - 1 {
            if i % 2 != 0 {
                let c_val = a_val * &b_val;
                let c_var = cs.alloc(|| format!("{}", i), || Ok(c_val))?;

                cs.enforce(
                    || format!("{}: a * b = c", i),
                    |lc| lc + a_var,
                    |lc| lc + b_var,
                    |lc| lc + c_var,
                );

                assignments.push((c_val, c_var));
                a_val = b_val;
                a_var = b_var;
                b_val = c_val;
                b_var = c_var;
            } else {
                let c_val = a_val + &b_val;
                let c_var = cs.alloc(|| format!("{}", i), || Ok(c_val))?;

                cs.enforce(
                    || format!("{}: a + b = c", i),
                    |lc| lc + a_var + b_var,
                    |lc| lc + CS::one(),
                    |lc| lc + c_var,
                );

                assignments.push((c_val, c_var));
                a_val = b_val;
                a_var = b_var;
                b_val = c_val;
                b_var = c_var;
            }
        }

        let mut a_lc = LinearCombination::zero();
        let mut b_lc = LinearCombination::zero();
        let mut c_val = F::zero();

        for (val, var) in assignments {
            a_lc = a_lc + var;
            b_lc = b_lc + var;
            c_val += &val;
        }
        c_val = c_val.square();

        let c_var = cs.alloc(|| "c_val", || Ok(c_val))?;

        cs.enforce(
            || "assignments.sum().square()",
            |_| a_lc,
            |_| b_lc,
            |lc| lc + c_var,
        );

        Ok(())
    }
}

fn bench_prover_circuit_high_densities(c: &mut Criterion) {
    let mut rng = XorShiftRng::seed_from_u64(1234567890u64);
    let mut group = c.benchmark_group("bench gro16 prover varying the number of constraints");

    let num_constraints = (15..=23).map(|i| 2usize.pow(i) - 3).collect::<Vec<_>>();

    for &num_constraints in num_constraints.iter() {
        println!(
            "************************{}************************",
            num_constraints
        );
        let params = {
            let c = BenchmarkHighDensities::<Fr> {
                num_constraints,
                _engine: PhantomData,
            };
            generate_random_parameters::<Bn382, _, _>(c, &mut rng).unwrap()
        };

        group.bench_with_input(
            BenchmarkId::from_parameter(num_constraints),
            &num_constraints,
            |bn, _constraints| {
                bn.iter(|| {
                    let c = BenchmarkHighDensities::<Fr> {
                        num_constraints,
                        _engine: PhantomData,
                    };
                    create_random_proof(c, &params, &mut rng).unwrap();
                });
            },
        );
    }
    group.finish();
}

fn bench_prover_circuit(c: &mut Criterion) {
    let mut rng = XorShiftRng::seed_from_u64(1234567890u64);
    let mut group =
        c.benchmark_group("bench gro16 prover varying the number of constraints high densities");

    let num_inputs = 2;
    let num_constraints = (15..=23)
        .map(|i| 2usize.pow(i) - (num_inputs + 1))
        .collect::<Vec<_>>();

    for &num_constraints in num_constraints.iter() {
        println!(
            "************************{}************************",
            num_constraints
        );
        let params = {
            let c = Benchmark::<Fr> {
                num_constraints,
                inputs: vec![None; num_inputs],
            };
            generate_random_parameters::<Bn382, _, _>(c, &mut rng).unwrap()
        };

        group.bench_with_input(
            BenchmarkId::from_parameter(num_constraints),
            &num_constraints,
            |bn, _constraints| {
                bn.iter_batched(
                    || {
                        let mut rng = XorShiftRng::seed_from_u64(num_constraints as u64);
                        let mut v = Vec::with_capacity(num_inputs);
                        for _ in 0..num_inputs {
                            v.push(Some(Fr::rand(&mut rng)))
                        }
                        v
                    },
                    |v| {
                        let c = Benchmark::<Fr> {
                            num_constraints,
                            inputs: v,
                        };
                        create_random_proof(c, &params, &mut rng).unwrap();
                    },
                    BatchSize::PerIteration,
                );
            },
        );
    }
    group.finish();
}

criterion_group!(
name = gro16_bench;
config = Criterion::default().sample_size(10);
targets = bench_prover_circuit, bench_prover_circuit_high_densities
);

criterion_main!(gro16_bench);
