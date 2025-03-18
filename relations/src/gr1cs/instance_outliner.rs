use ark_ff::Field;

use crate::gr1cs::LinearCombination;

use super::{
    lc, predicate::polynomial_constraint::SR1CS_PREDICATE_LABEL, ConstraintSystem,
    Label, SynthesisError, Variable, R1CS_PREDICATE_LABEL,
};
use ark_std::{collections::BTreeMap, rc::Rc, vec::Vec};
use core::fmt::Debug;


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
    pub func: Rc<
        dyn Fn(
            &mut ConstraintSystem<F>,
            BTreeMap<Variable, Variable>,
        ) -> Result<(), SynthesisError>,
    >,
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
    instance_witness_map: BTreeMap<Variable, Variable>,
) -> crate::gr1cs::Result<()> {
    let one_witt = instance_witness_map.get(&Variable::one()).unwrap();
    // Now, enforce the equality between the instance and the corresponding witness
    // variable This is done by iterating over the instance-witness map
    // which contains the unique instance-witness pairs The equality
    // constraints are enforced with r1cs constraints, it is assumed that a
    // constraint system has a default r1cs predicate registered
    for (instance, witness) in instance_witness_map.iter() {
        let r1cs_constraint = if instance.is_one() {
            vec![
                LinearCombination::from(*witness),
                LinearCombination::from(*witness),
                LinearCombination::from(*instance),
            ]
        } else {
            vec![
                LinearCombination::from(*one_witt),
                LinearCombination::from(*witness),
                LinearCombination::from(*instance),
            ]
        };

        cs.enforce_constraint(R1CS_PREDICATE_LABEL, r1cs_constraint)?;
    }

    Ok(())
}

/// The outlining strategy for Square R1CS constraints.
pub fn outline_sr1cs<F: Field>(
    cs: &mut ConstraintSystem<F>,
    instance_witness_map: BTreeMap<Variable, Variable>,
) -> crate::gr1cs::Result<()> {
    // Now, enforce the equality between the instance and the corresponding witness
    // variable This is done by iterating over the instance-witness map
    // which contains the unique instance-witness pairs The equality
    // constraints are enforced with r1cs constraints, it is assumed that a
    // constraint system has a default r1cs predicate registered
    for (instance, witness) in instance_witness_map.iter() {
        let cnstr: Vec<LinearCombination<F>> = vec![lc!() + instance - witness, lc!()];
        cs.enforce_constraint(SR1CS_PREDICATE_LABEL, cnstr)?;
    }

    Ok(())
}
