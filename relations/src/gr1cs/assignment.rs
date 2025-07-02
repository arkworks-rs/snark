use ark_ff::Field;
use ark_std::vec::Vec;

use crate::{
    gr1cs::{field_interner::FieldInterner, lc_map::LcMap, Variable},
    utils::variable::VarKind,
};

/// Assignments for a GR1CS constraint system.
#[derive(Debug, Clone)]
pub struct Assignments<F> {
    /// Assignments to the public input variables. This is empty if `self.mode
    /// == SynthesisMode::Setup`.
    pub instance_assignment: Vec<F>,
    /// Assignments to the private input variables. This is empty if `self.mode
    /// == SynthesisMode::Setup`.
    pub witness_assignment: Vec<F>,
    /// A cache for the linear combination assignments. It shows evaluation
    /// result of each linear combination
    pub lc_assignment: Vec<F>,
}

impl<F: Field> Assignments<F> {
    /// Obtain the assignment corresponding to the `Variable` `v`.
    #[inline]
    pub fn assigned_value(&self, v: Variable) -> Option<F> {
        let idx = v.index();
        match v.kind() {
            VarKind::Zero => Some(F::ZERO),
            VarKind::One => Some(F::ONE),
            VarKind::Instance => self.instance_assignment.get(idx?).copied(),
            VarKind::Witness => self.witness_assignment.get(idx?).copied(),
            VarKind::SymbolicLc => self.lc_assignment.get(idx?).copied(),
        }
    }

    /// Evaluate the linear combination `lc` with the assigned values and return
    /// the result.
    #[inline]
    pub(super) fn eval_lc(
        &self,
        lc: usize,
        lc_map: &LcMap<F>,
        f_interner: &FieldInterner<F>,
    ) -> Option<F> {
        lc_map.get(lc).map(|lc| {
            lc.map(|(&coeff, &var)| {
                f_interner.value(coeff).unwrap() * self.assigned_value(var).unwrap()
            })
            .sum()
        })
    }
}
