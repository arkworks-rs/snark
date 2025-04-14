mod circuit1;
mod circuit2;

use crate::{
    ark_std::string::ToString,
    gr1cs::instance_outliner::{outline_r1cs, InstanceOutliner},
};
use ark_std::rc::Rc;

use super::*;
use ark_ff::{One, Zero};
use ark_test_curves::bls12_381::Fr;
pub use circuit1::Circuit1;
pub use circuit2::Circuit2;
use constraint_system::ConstraintSystem;

#[test]
fn test_circuit1_sat() {
    let c = Circuit1 {
        x1: Fr::from(1u8),
        x2: Fr::from(2u8),
        x3: Fr::from(3u8),
        x4: Fr::from(0u8),
        x5: Fr::from(1255254u32),
        w1: Fr::from(4u8),
        w2: Fr::from(2u8),
        w3: Fr::from(5u8),
        w4: Fr::from(29u8),
        w5: Fr::from(28u8),
        w6: Fr::from(10u8),
        w7: Fr::from(57u8),
        w8: Fr::from(22022u32),
    };
    let cs = ConstraintSystem::<Fr>::new_ref();
    c.clone().generate_constraints(cs.clone()).unwrap();
    cs.finalize();
    assert!(cs.is_satisfied().unwrap());

    let cs = ConstraintSystem::<Fr>::new_ref();
    cs.set_optimization_goal(OptimizationGoal::Constraints);
    c.clone().generate_constraints(cs.clone()).unwrap();
    cs.finalize();
    assert!(cs.is_satisfied().unwrap());

    let cs = ConstraintSystem::<Fr>::new_ref();
    cs.set_optimization_goal(OptimizationGoal::Weight);
    c.clone().generate_constraints(cs.clone()).unwrap();
    cs.finalize();
    assert!(cs.is_satisfied().unwrap());
}

/// Test the circuit with non-satisfying inputs
/// The first input is changed comparing to sat test
#[test]
fn test_circuit1_non_sat() {
    // use crate::gr1cs::TracingMode;
    // use tracing_subscriber::prelude::*;
    // let mut layer = ConstraintLayer::default();
    // layer.mode = TracingMode::OnlyConstraints;
    // let subscriber = tracing_subscriber::Registry::default().with(layer);
    // tracing::subscriber::set_global_default(subscriber).unwrap();

    let c = Circuit1 {
        x1: Fr::from(4u8),
        x2: Fr::from(2u8),
        x3: Fr::from(3u8),
        x4: Fr::from(0u8),
        x5: Fr::from(1255254u32),
        w1: Fr::from(4u8),
        w2: Fr::from(2u8),
        w3: Fr::from(5u8),
        w4: Fr::from(29u8),
        w5: Fr::from(28u8),
        w6: Fr::from(10u8),
        w7: Fr::from(57u8),
        w8: Fr::from(22022u32),
    };

    let cs = ConstraintSystem::<Fr>::new_ref();
    c.generate_constraints(cs.clone()).unwrap();
    assert!(!cs.is_satisfied().unwrap()); // This fails as expected.
}

#[test]
fn test_circuit1_matrices() {
    let c = Circuit1 {
        x1: Fr::zero(),
        x2: Fr::zero(),
        x3: Fr::zero(),
        x4: Fr::zero(),
        x5: Fr::zero(),
        w1: Fr::zero(),
        w2: Fr::zero(),
        w3: Fr::zero(),
        w4: Fr::zero(),
        w5: Fr::zero(),
        w6: Fr::zero(),
        w7: Fr::zero(),
        w8: Fr::zero(),
    };
    let cs = ConstraintSystem::<Fr>::new_ref();
    c.clone().generate_constraints(cs.clone()).unwrap();
    assert_eq!(Circuit1::get_matrices(), cs.clone().to_matrices().unwrap());
    cs.set_instance_outliner(InstanceOutliner {
        pred_label: R1CS_PREDICATE_LABEL.to_string(),
        func: Rc::new(outline_r1cs),
    });
    cs.finalize();
}

#[test]
fn test_circuit1_instance_outlined() {
    let c = Circuit1 {
        x1: Fr::zero(),
        x2: Fr::zero(),
        x3: Fr::zero(),
        x4: Fr::zero(),
        x5: Fr::zero(),
        w1: Fr::zero(),
        w2: Fr::zero(),
        w3: Fr::zero(),
        w4: Fr::zero(),
        w5: Fr::zero(),
        w6: Fr::zero(),
        w7: Fr::zero(),
        w8: Fr::zero(),
    };
    let cs = ConstraintSystem::<Fr>::new_ref();
    c.clone().generate_constraints(cs.clone()).unwrap();
    let num_instance = cs.num_instance_variables();
    let prev_num_witness = cs.num_witness_variables();
    cs.set_instance_outliner(InstanceOutliner {
        pred_label: R1CS_PREDICATE_LABEL.to_string(),
        func: Rc::new(outline_r1cs),
    });
    cs.finalize();
    let new_num_witness = cs.num_witness_variables();
    assert_eq!(num_instance, new_num_witness - prev_num_witness);
}

/// This is the legacy test for R1CS from the previous version of the library
#[test]
fn test_circuit2_matrices() {
    let c = Circuit2 {
        a: Fr::one(),
        b: Fr::one(),
        c: Fr::one() + Fr::one(),
    };
    let cs = ConstraintSystem::<Fr>::new_ref();
    c.clone().generate_constraints(cs.clone()).unwrap();
    cs.finalize();
    assert_eq!(Circuit2::get_matrices(), cs.clone().to_matrices().unwrap());
}
