#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use ark_ff::{Field, UniformRand};
use ark_relations::gr1cs::{
    ConstraintSynthesizer, ConstraintSystem, ConstraintSystemRef, LinearCombination,
    OptimizationGoal, SynthesisMode,
};
use ark_std::rand::{Rng, SeedableRng};
use ark_test_curves::bls12_381::Fr;
// use jemallocator::Jemalloc;
use mimalloc::MiMalloc;

pub const NUM_COEFFS_IN_LC: usize = 10;

struct BenchCircuit {
    a: Fr,
    b: Fr,
    c: Fr,
    num_constraints: usize,
}

impl ConstraintSynthesizer<Fr> for BenchCircuit {
    fn generate_constraints(self, cs: ConstraintSystemRef<Fr>) -> ark_relations::gr1cs::Result<()> {
        let a = cs.new_witness_variable(|| Ok(self.a)).unwrap();
        let b = cs.new_witness_variable(|| Ok(self.b)).unwrap();
        let c = cs.new_witness_variable(|| Ok(self.c)).unwrap();
        let mut variables = vec![a, b, c];
        variables.reserve(3 * self.num_constraints);
        let mut lcs = Vec::with_capacity(self.num_constraints);

        let mut rng_a = ark_std::rand::rngs::StdRng::seed_from_u64(0_u64);
        let mut rng_b = ark_std::rand::rngs::StdRng::seed_from_u64(rng_a.gen::<u64>());
        let mut rng_c = ark_std::rand::rngs::StdRng::seed_from_u64(rng_a.gen::<u64>());

        for i in 0..self.num_constraints {
            let cur_num_vars = ark_std::cmp::min(variables.len(), 10);
            let lower = variables.len().saturating_sub(cur_num_vars);
            let upper = variables.len();

            let a_i_lc_size = rng_a.gen_range(1..=NUM_COEFFS_IN_LC);
            let mut a_i = || {
                LinearCombination(
                    (0..a_i_lc_size)
                        .map(|_| (Fr::ONE, variables[rng_a.gen_range(lower..upper)]))
                        .collect::<Vec<_>>(),
                )
            };

            let b_i_lc_size = rng_b.gen_range(1..=NUM_COEFFS_IN_LC);
            let b_i = || {
                LinearCombination(
                    (0..b_i_lc_size)
                        .map(|_| (Fr::ONE, variables[rng_b.gen_range(lower..upper)]))
                        .collect::<Vec<_>>(),
                )
            };

            let c_i = variables[rng_c.gen_range(lower..upper)];

            if i % 2 == 0 {
                let extra_lc = || {
                    LinearCombination(
                        (0..a_i_lc_size)
                            .map(|_| (Fr::ONE, variables[rng_c.gen_range(lower..upper)]))
                            .collect::<Vec<_>>(),
                    )
                };
                let extra_lc = cs.new_lc(extra_lc)?;
                lcs.push(extra_lc);
                cs.enforce_r1cs_constraint(|| a_i() + extra_lc, b_i, || c_i.into())?;
            } else {
                cs.enforce_r1cs_constraint(a_i, b_i, || c_i.into())?;
            };
            let v1 = cs.new_witness_variable(|| Ok(self.a))?;
            let v2 = cs.new_witness_variable(|| Ok(self.a))?;
            let v3 = cs.new_witness_variable(|| Ok(self.a))?;
            variables.push(v1);
            variables.push(v2);
            variables.push(v3);
        }

        Ok(())
    }
}

fn main() {
    for log_num_constraints in 23..24 {
        let circuit = BenchCircuit {
            a: Fr::rand(&mut ark_std::test_rng()),
            b: Fr::rand(&mut ark_std::test_rng()),
            c: Fr::rand(&mut ark_std::test_rng()),
            num_constraints: 2usize.pow(log_num_constraints),
        };
        let cs = ConstraintSystem::<Fr>::new_ref();
        cs.set_optimization_goal(OptimizationGoal::Constraints);
        cs.set_mode(SynthesisMode::Prove {
            construct_matrices: false,
            generate_lc_assignments: false,
        });
        let start = std::time::Instant::now();
        circuit.generate_constraints(cs.clone()).unwrap();
        // cs.finalize();
        let elapsed = start.elapsed();
        println!(
            "Synthesizing 2^{} constraints took {:?}",
            log_num_constraints,
            elapsed.as_secs_f32()
        );
    }
}
