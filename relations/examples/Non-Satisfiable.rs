use ark_ff::Field;
use ark_relations::{
    gr1cs::{
        trace::{ConstraintLayer, TracingMode},
        ConstraintSystem, ConstraintSystemRef, LinearCombination, SynthesisError, Variable,
    },
    ns,
};
use ark_test_curves::bls12_381::Fr;
use tracing::{info_span, Instrument};
use tracing_subscriber::{fmt, layer::SubscriberExt, Registry};

fn main() {
    // Set up tracing with a ConstraintLayer
    let constraint_layer = ConstraintLayer::new(TracingMode::All);
    let subscriber = Registry::default()
            .with(fmt::layer()) // Optional: Log formatted output to stdout
            .with(constraint_layer);

    tracing::subscriber::set_global_default(subscriber)
        .expect("Failed to set global default subscriber");

    // Initialize a constraint system
    let cs = ConstraintSystem::<Fr>::new_ref();

    // Create input variables
    let p1 = cs.new_input_variable(|| Ok(Fr::from(3u32))).unwrap();
    let p2 = cs.new_input_variable(|| Ok(Fr::from(4u32))).unwrap();
    let p3 = cs.new_input_variable(|| Ok(Fr::from(6u32))).unwrap();
    let p4 = cs.new_input_variable(|| Ok(Fr::from(7u32))).unwrap();

    // Create witness variables
    // Note that w3 has a wrong value
    let w1 = cs.new_witness_variable(|| Ok(Fr::from(2u32))).unwrap();
    let w2 = cs.new_witness_variable(|| Ok(Fr::from(5u32))).unwrap();
    let w3 = cs.new_witness_variable(|| Ok(Fr::from(8u32))).unwrap();
    let w4 = cs.new_witness_variable(|| Ok(Fr::from(9u32))).unwrap();
    // Create expected result as a witness
    let expected_result = cs.new_input_variable(|| Ok(Fr::from(198))).unwrap();

    // Build the circuit
    let _result_var =
        generate_constraints(cs.clone(), p1, p2, p3, p4, w1, w2, w3, w4, expected_result).unwrap();
    let trace = cs.which_predicate_is_unsatisfied().unwrap().unwrap();
    println!("This is the trace of non-satisfied scenario, Check out the trace:\n{}", trace)
}
// Function to enforce the overall constraints by combining subcircuits
#[tracing::instrument(target = "gr1cs")]
fn generate_constraints<F: Field>(
    cs: ConstraintSystemRef<F>,
    p1: Variable,
    p2: Variable,
    p3: Variable,
    p4: Variable,
    w1: Variable,
    w2: Variable,
    w3: Variable,
    w4: Variable,
    expected_result: Variable,
) -> Result<Variable, SynthesisError> {
    // Subcircuit 1
    let subcircuit1_result = enforce_subcircuit1(cs.clone(), p1, p2, w1, w2)?;

    // Subcircuit 2
    let subcircuit2_result = enforce_subcircuit2(cs.clone(), p3, p4, w3, w4)?;

    // Final operation: sum the results of the two subcircuits
    let final_result = enforce_addition(cs.clone(), subcircuit1_result, subcircuit2_result)?;
    // Verify that the final result matches the expected result
    cs.enforce_r1cs_constraint(
        LinearCombination::from(final_result),
        LinearCombination::from(Variable::One),
        LinearCombination::from(expected_result),
    )?;

    Ok(final_result)
}

// Enforces the constraints for Subcircuit 1
#[tracing::instrument(target = "gr1cs")]
fn enforce_subcircuit1<F: Field>(
    cs: ConstraintSystemRef<F>,
    p1: Variable,
    p2: Variable,
    w1: Variable,
    w2: Variable,
) -> Result<Variable, SynthesisError> {
    // let cs = ns!(cs, "subcircuit1").cs();

    // Multiplication gate: p1 * p2
    let product = enforce_multiplication(cs.clone(), p1, p2)?;

    // Addition gate: w1 + w2
    let sum = enforce_addition(cs.clone(), w1, w2)?;

    // Final multiplication: sum * product
    let result = enforce_multiplication(cs, sum, product)?;

    Ok(result)
}

// Enforces the constraints for Subcircuit 2
#[tracing::instrument(target = "gr1cs")]
fn enforce_subcircuit2<F: Field>(
    cs: ConstraintSystemRef<F>,
    p3: Variable,
    p4: Variable,
    w3: Variable,
    w4: Variable,
) -> Result<Variable, SynthesisError> {
    // let cs = ns!(cs, "subcircuit2").cs();

    // Multiplication gate 1: p3 * w3
    let product1 = enforce_multiplication(cs.clone(), p3, p4)?;

    // Multiplication gate 2: p4 * w4
    let product2 = enforce_multiplication(cs.clone(), w3, w4)?;

    // Final addition: product1 + product2
    let result = enforce_addition(cs, product1, product2)?;

    Ok(result)
}

// Function to enforce a multiplication constraint
#[tracing::instrument(target = "gr1cs")]
fn enforce_multiplication<F: Field>(
    cs: ConstraintSystemRef<F>,
    left: Variable,
    right: Variable,
) -> Result<Variable, SynthesisError> {
    // let cs = ns!(cs, "multiplication").cs();

    let product = cs.new_witness_variable(|| {
        Ok(cs.assigned_value(left).unwrap() * cs.assigned_value(right).unwrap())
    })?; // Placeholder for product

    cs.enforce_r1cs_constraint(
        LinearCombination::from(left),
        LinearCombination::from(right),
        LinearCombination::from(product),
    )?;

    Ok(product)
}

// Function to enforce an addition constraint
#[tracing::instrument(target = "gr1cs")]
fn enforce_addition<F: Field>(
    cs: ConstraintSystemRef<F>,
    left: Variable,
    right: Variable,
) -> Result<Variable, SynthesisError> {
    // let cs = ns!(cs, "addition").cs();

    // Instead of + we use *, This is intentioanlly made to fail
    let sum = cs.new_witness_variable(|| {
        Ok(cs.assigned_value(left).unwrap() * cs.assigned_value(right).unwrap())
    })?; // Placeholder for sum

    cs.enforce_r1cs_constraint(
        LinearCombination::from(left) + LinearCombination::from(right),
        LinearCombination::from(Variable::One),
        LinearCombination::from(sum),
    )?;

    Ok(sum)
}



