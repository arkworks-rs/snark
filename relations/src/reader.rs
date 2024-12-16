use crate::r1cs::{
    ConstraintSystem, ConstraintSystemRef, LinearCombination, SynthesisError, Variable,
};
use ark_circom::{CircomBuilder, CircomCircuit, CircomConfig};
use ark_ff::PrimeField;
use ark_std::path::Path;

pub fn read_constraint_system<F: PrimeField>(
    r1cs_file: impl AsRef<Path>,
    wasm_file: impl AsRef<Path>,
) -> ConstraintSystem<F> {
    // Load the WASM and R1CS for witness and proof generation
    let cfg = CircomConfig::<F>::new(wasm_file, r1cs_file).unwrap();

    let builder = CircomBuilder::new(cfg);
    let circom = builder.setup();

    let cs = ConstraintSystem::<F>::new_ref();

    // TODO: replace with `circom.generate_constraints(cs.clone())` once the
    // dependency issue is fixed.
    generate_constraints(circom, cs.clone()).unwrap();

    cs.into_inner().unwrap()
}

// TODO: currently, CircomCircuit::generate_constraints() cannot be used due to
// inconsistent imports. This branch of `ark-relations` has
// `circom-compat/release-0.5` as a dev-dependency, and
// `circom-compat/release-0.5` has `ark-relations v0.5.0` as a dependency, which
// currently does not exist, therefore
// `ark-circom::ark_relations::r1cs::ConstraintSystem` is not the same structure
// as `crate::r1cs::ConstraintSystem`.
fn generate_constraints<F: PrimeField>(
    circom: CircomCircuit<F>,
    cs: ConstraintSystemRef<F>,
) -> Result<(), SynthesisError> {
    let witness = &circom.witness;
    let wire_mapping = &circom.r1cs.wire_mapping;

    // Start from 1 because Arkworks implicitly allocates One for the first input
    for i in 1..circom.r1cs.num_inputs {
        cs.new_input_variable(|| {
            Ok(match witness {
                None => F::from(1u32),
                Some(w) => match wire_mapping {
                    Some(m) => w[m[i]],
                    None => w[i],
                },
            })
        })?;
    }

    for i in 0..circom.r1cs.num_aux {
        cs.new_witness_variable(|| {
            Ok(match witness {
                None => F::from(1u32),
                Some(w) => match wire_mapping {
                    Some(m) => w[m[i + circom.r1cs.num_inputs]],
                    None => w[i + circom.r1cs.num_inputs],
                },
            })
        })?;
    }

    let make_index = |index| {
        if index < circom.r1cs.num_inputs {
            Variable::Instance(index)
        } else {
            Variable::Witness(index - circom.r1cs.num_inputs)
        }
    };
    let make_lc = |lc_data: &[(usize, F)]| {
        lc_data.iter().fold(
            LinearCombination::<F>::zero(),
            |lc: LinearCombination<F>, (index, coeff)| lc + (*coeff, make_index(*index)),
        )
    };

    for constraint in &circom.r1cs.constraints {
        cs.enforce_constraint(
            make_lc(&constraint.0),
            make_lc(&constraint.1),
            make_lc(&constraint.2),
        )?;
    }

    Ok(())
}
