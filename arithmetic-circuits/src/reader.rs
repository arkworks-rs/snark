use ark_circom::{CircomBuilder, CircomConfig};
use ark_ff::PrimeField;
use ark_relations::r1cs::{ConstraintSynthesizer, ConstraintSystem};
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
    circom.generate_constraints(cs.clone()).unwrap();
    cs.into_inner().unwrap()
}
