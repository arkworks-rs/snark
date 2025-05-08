use ark_ff::Field;

use super::{
    lc, predicate::polynomial_constraint::SR1CS_PREDICATE_LABEL, ConstraintSystem, Label,
    SynthesisError, Variable, R1CS_PREDICATE_LABEL,
};
use ark_std::rc::Rc;
use core::fmt::Debug;

/// A type alias for the instance outlining function
pub type InstanceOutliningFunction<F> = dyn Fn(&mut ConstraintSystem<F>, &[Variable]) -> Result<(), SynthesisError>;

/// An instance outliner is a strategy for reducing the number of constraints
/// that public input/instance variables are involved in.
/// It does this as follows:
/// 1. Create new public input variables that correspond to the original input variables
/// 2. Replace the original input variables with corresponding new witness variables
/// 3. Enforce equality between the new input variables and the new witness variables
/// 4. In every constraint that involves the original input variables, replace them with the new witness variables.
#[derive(Clone)]
pub struct InstanceOutliner<F: Field> {
    /// The label for the predicate that is used to enforce equality between
    /// the new input/instance variables and the new witness variables.
    pub pred_label: Label,
    /// The strategy for outlining the instance variables
    /// It takes as input the constraint system, and a map from the new
    /// instance variables to the new witness variables.
    pub func: Rc<InstanceOutliningFunction<F>>,
}

impl<F: Field> Debug for InstanceOutliner<F> {
    fn fmt(&self, f: &mut ark_std::fmt::Formatter<'_>) -> ark_std::fmt::Result {
        write!(
            f,
            "InstanceOutliner {{ pred_label: {:?} }}",
            self.pred_label
        )
    }
}

/// The outlining strategy for R1CS constraints.
pub fn outline_r1cs<F: Field>(
    cs: &mut ConstraintSystem<F>,
    instance_witness_map: &[Variable],
) -> crate::gr1cs::Result<()> {
    // Now, enforce the equality between the instance and the corresponding witness
    // variable This is done by iterating over the instance-witness map
    // which contains the unique instance-witness pairs The equality
    // constraints are enforced with r1cs constraints, it is assumed that a
    // constraint system has a default r1cs predicate registered
    let one = instance_witness_map[0];
    cs.enforce_constraint(
        R1CS_PREDICATE_LABEL,
        [lc!() + one, lc!() + one, lc!() + Variable::One],
    )?;
    for (instance, witness) in instance_witness_map.iter().enumerate().skip(1) {
        cs.enforce_constraint(
            R1CS_PREDICATE_LABEL,
            [
                lc!() + one,
                lc!() + *witness,
                lc!() + Variable::Instance(instance),
            ],
        )?;
    }

    Ok(())
}

/// The outlining strategy for Square R1CS constraints.
pub fn outline_sr1cs<F: Field>(
    cs: &mut ConstraintSystem<F>,
    instance_witness_map: &[Variable],
) -> crate::gr1cs::Result<()> {
    // Now, enforce the equality between the instance and the corresponding witness
    // variable This is done by iterating over the instance-witness map
    // which contains the unique instance-witness pairs The equality
    // constraints are enforced with r1cs constraints, it is assumed that a
    // constraint system has a default r1cs predicate registered
    for (instance, witness) in instance_witness_map.iter().enumerate() {
        cs.enforce_constraint(
            SR1CS_PREDICATE_LABEL,
            [lc!() + Variable::Instance(instance) - witness, lc!()],
        )?;
    }

    Ok(())
} 
