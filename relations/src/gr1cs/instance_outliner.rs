use ark_ff::{Field, Zero};

use crate::gr1cs::LinearCombination;

use super::{
    lc, predicate::polynomial_constraint::SR1CS_PREDICATE_LABEL, ConstraintSystem,
    ConstraintSystemRef, Label, SynthesisError, Variable, R1CS_PREDICATE_LABEL,
};
use ark_std::{collections::BTreeMap, rc::Rc, vec::Vec};
use core::fmt::Debug;
#[derive(Clone)]
pub struct InstanceOutliner<F: Field> {
    pub pred_label: Label,
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
