//! This module contains the implementation of the `ConstraintSystem` struct.
//! a constraint system contains multiple predicate constraint systems,
//! each of which enforce have separate predicates and constraints. For more information about the terminology and the structure of the constraint system, refer to section 3.3 of https://eprint.iacr.org/2024/1245

use super::{
    instance_outliner::InstanceOutliner,
    predicate::{
        polynomial_constraint::{R1CS_PREDICATE_LABEL, SR1CS_PREDICATE_LABEL},
        Predicate, PredicateConstraintSystem,
    },
    ConstraintSystemRef, Label, OptimizationGoal, SynthesisMode,
};
#[cfg(feature = "std")]
use crate::gr1cs::ConstraintTrace;
use crate::{
    gr1cs::{LinearCombination, Matrix, SynthesisError, Variable},
    utils::variable::VarKind,
};
use ark_ff::Field;
use ark_std::{
    any::{Any, TypeId},
    boxed::Box,
    cell::RefCell,
    collections::BTreeMap,
    format,
    rc::Rc,
    string::{String, ToString},
    vec::Vec,
};

use crate::utils::IndexMap;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
///////////////////////////////////////////////////////////////////////////////////////

/// A GR1CS `ConstraintSystem`. Enforces constraints of the form  
/// `L_i(⟨M_{i,1}, z⟩ , ⟨M_{i,2}, z⟩ , ..., ⟨M_{i,t_i}, z⟩)=0`
/// More Information: https://eprint.iacr.org/2024/1245
#[derive(Debug, Clone)]
pub struct ConstraintSystem<F: Field> {
    /// The mode in which the constraint system is operating. `self` can either
    /// be in setup mode (i.e., `self.mode == SynthesisMode::Setup`) or in
    /// proving mode (i.e., `self.mode == SynthesisMode::Prove`). If we are
    /// in proving mode, then we have the additional option of whether or
    /// not to construct the A, B, and C matrices of the constraint system
    mode: SynthesisMode,

    /// The number of variables that are "public inputs" to the constraint
    /// system.
    #[doc(hidden)]
    pub num_instance_variables: usize,

    /// The number of variables that are "private inputs" to the constraint
    /// system.
    #[doc(hidden)]
    pub num_witness_variables: usize,

    /// The number of linear combinations
    #[doc(hidden)]
    pub num_linear_combinations: usize,

    /// The parameter we aim to minimize in this constraint system (either the
    /// number of constraints or their total weight).
    optimization_goal: OptimizationGoal,

    /// If true, the constraint system will outline the instances. Outlining the instances is a technique used for verification succinctness in some SNARKs like Polymath, Garuda, and Pari. For more information, refer to https://eprint.iacr.org/2024/1245
    /// It assigns a witness variable to each instance variable and enforces the
    /// equality of the instance and witness variables. Then only uses the
    /// witness variables in the constraints.
    instance_outliner: Option<InstanceOutliner<F>>,

    /// Assignments to the input, witness, and lc variables. This is empty if `self.mode
    /// == SynthesisMode::Setup`.
    pub assignments: Assignments<F>,

    /// Map for gadgets to cache computation results.
    pub cache_map: Rc<RefCell<BTreeMap<TypeId, Box<dyn Any>>>>,

    /// A data structure to store the linear combinations. We use map because
    /// it's easier to inline and outline the linear combinations.
    #[doc(hidden)]
    pub lc_map: LcMap<F>,

    /// A map from the the predicate labels to the predicates
    #[doc(hidden)]
    pub predicate_constraint_systems: BTreeMap<Label, PredicateConstraintSystem<F>>,

    /// data structure to store the traces for each predicate
    #[cfg(feature = "std")]
    pub predicate_traces: BTreeMap<Label, Vec<Option<ConstraintTrace>>>,
}

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
    fn eval_lc(&self, lc: usize, lc_map: &LcMap<F>) -> Option<F> {
        let acc = lc_map[lc]
            .iter()
            .map(|(coeff, var)| *coeff * self.assigned_value(*var).unwrap())
            .sum();
        Some(acc)
    }
}

impl<F: Field> Default for ConstraintSystem<F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<F: Field> ConstraintSystem<F> {
    /// Create a new and empty `ConstraintSystem<F>`.
    /// Note that by default, the constraint system is
    /// registered with the R1CS predicate.
    pub fn new() -> Self {
        let mut lc_map = LcMap::new();
        lc_map.push(LinearCombination::zero());

        let mut cs = Self {
            num_instance_variables: 1,
            num_witness_variables: 0,
            num_linear_combinations: 1,
            instance_outliner: None,
            predicate_constraint_systems: BTreeMap::new(),
            assignments: Assignments {
                instance_assignment: vec![F::one()],
                witness_assignment: Vec::new(),
                lc_assignment: vec![F::zero()],
            },
            cache_map: Rc::new(RefCell::new(BTreeMap::new())),
            lc_map,
            mode: SynthesisMode::Prove {
                construct_matrices: true,
                generate_lc_assignments: true,
            },
            optimization_goal: OptimizationGoal::None,
            #[cfg(feature = "std")]
            predicate_traces: BTreeMap::new(),
        };
        let r1cs_constraint_system = PredicateConstraintSystem::new_r1cs().unwrap();
        let _ = cs.register_predicate(R1CS_PREDICATE_LABEL, r1cs_constraint_system);
        cs
    }

    /// Create a new `ConstraintSystemRef<F>`.
    pub fn new_ref() -> ConstraintSystemRef<F> {
        ConstraintSystemRef::new(Self::new())
    }

    /// Returns a mapping from predicate label to number of constraints for that
    /// predicate
    pub fn get_all_predicates_num_constraints(&self) -> IndexMap<Label, usize> {
        self.predicate_constraint_systems
            .iter()
            .map(|(label, predicate)| (label.clone(), predicate.num_constraints()))
            .collect()
    }

    /// Returns the number of constraints for a given predicate
    pub fn get_predicate_num_constraints(&self, predicate_label: &str) -> Option<usize> {
        self.predicate_constraint_systems
            .get(predicate_label)
            .map(|predicate| predicate.num_constraints())
    }

    /// Returns a mapping from predicate label to arity for that predicate
    pub fn get_all_predicate_arities(&self) -> IndexMap<Label, usize> {
        self.predicate_constraint_systems
            .iter()
            .map(|(label, predicate)| (label.clone(), predicate.get_arity()))
            .collect()
    }

    /// Returns the type of the predicate with the given label
    pub fn get_predicate_arity(&self, predicate_label: &str) -> Option<usize> {
        self.predicate_constraint_systems
            .get(predicate_label)
            .map(|predicate| predicate.get_arity())
    }

    /// Returns a mapping from predicate labels to their types
    pub fn get_all_predicate_types(&self) -> BTreeMap<Label, Predicate<F>> {
        self.predicate_constraint_systems
            .iter()
            .map(|(label, predicate)| (label.clone(), predicate.get_predicate().clone()))
            .collect()
    }

    /// Returns the type of the predicate with the given label
    pub fn get_predicate_type(&self, predicate_label: &str) -> Option<Predicate<F>> {
        self.predicate_constraint_systems
            .get(predicate_label)
            .map(|predicate| predicate.get_predicate().clone())
    }

    /// Returns the assignment to the public input variables of the constraint
    pub fn instance_assignment(&self) -> crate::gr1cs::Result<&[F]> {
        if self.is_in_setup_mode() {
            return Err(SynthesisError::AssignmentMissing);
        }
        Ok(&self.assignments.instance_assignment)
    }

    /// Returns the assignment to the private input variables of the constraint
    pub fn witness_assignment(&self) -> crate::gr1cs::Result<&[F]> {
        if self.is_in_setup_mode() {
            return Err(SynthesisError::AssignmentMissing);
        }
        Ok(&self.assignments.witness_assignment)
    }

    /// Returns the number of constraints which is the sum of the number of
    /// constraints in each predicate.
    pub fn num_constraints(&self) -> usize {
        self.predicate_constraint_systems
            .values()
            .map(|p| p.num_constraints())
            .sum()
    }

    /// Returns the number of instance variables.
    pub fn num_instance_variables(&self) -> usize {
        self.num_instance_variables
    }

    /// Returns the number of witness variables.
    pub fn num_witness_variables(&self) -> usize {
        self.num_witness_variables
    }

    /// Returns the number of witness variables.
    pub fn num_variables(&self) -> usize {
        self.num_witness_variables() + self.num_instance_variables()
    }

    /// Returns the number of predicates in the constraint system
    pub fn num_predicates(&self) -> usize {
        self.predicate_constraint_systems.len()
    }

    /// Enforce a constraint in the constraint system. It takes a
    /// predicate name and enforces a vector of linear combinations of the
    /// length of the arity of the predicate enforces the constraint.
    #[inline]
    pub fn enforce_constraint(
        &mut self,
        predicate_label: &str,
        lcs: impl IntoIterator<
            Item = Box<dyn FnOnce() -> LinearCombination<F>>,
            IntoIter: ExactSizeIterator,
        >,
    ) -> crate::gr1cs::Result<()> {
        if !self.has_predicate(predicate_label) {
            return Err(SynthesisError::PredicateNotFound);
        }

        if self.should_construct_matrices() {
            let should_generate_lc_assignments = self.should_generate_lc_assignments();
            let lc_map = &mut self.lc_map;
            let num_lcs = &mut self.num_linear_combinations;
            let assignments = &mut self.assignments;

            let lc_indices = lcs.into_iter().map(|lc| {
                let lc = lc();
                Self::new_lc_add_helper(
                    num_lcs,
                    lc_map,
                    should_generate_lc_assignments,
                    assignments,
                    lc,
                )
                .unwrap()
            });

            let predicate = self
                .predicate_constraint_systems
                .get_mut(predicate_label)
                .unwrap();

            predicate.enforce_constraint(lc_indices)?;
        }

        #[cfg(feature = "std")]
        match self.predicate_traces.get_mut(predicate_label) {
            Some(traces) => traces.push(ConstraintTrace::capture()),
            None => {
                eprintln!("Constraint trace requires adding the predicate constraint trace")
            },
        }

        Ok(())
    }

    /// Enforce a constraint for a predicate with arity 2.
    pub fn enforce_constraint_arity_2(
        &mut self,
        predicate_label: &str,
        a: impl FnOnce() -> LinearCombination<F>,
        b: impl FnOnce() -> LinearCombination<F>,
    ) -> crate::gr1cs::Result<()> {
        if !self.has_predicate(predicate_label) {
            return Err(SynthesisError::PredicateNotFound);
        }

        if self.should_construct_matrices() {
            let a = self.new_constraint_lc(a)?;
            let b = self.new_constraint_lc(b)?;

            let predicate = self
                .predicate_constraint_systems
                .get_mut(predicate_label)
                .unwrap();

            predicate.enforce_constraint([a, b])?;
        }

        #[cfg(feature = "std")]
        if let Some(traces) = self.predicate_traces.get_mut(predicate_label) {
            traces.push(ConstraintTrace::capture())
        }

        Ok(())
    }

    /// Enforce a constraint for a predicate with arity 3.
    pub fn enforce_constraint_arity_3(
        &mut self,
        predicate_label: &str,
        a: impl FnOnce() -> LinearCombination<F>,
        b: impl FnOnce() -> LinearCombination<F>,
        c: impl FnOnce() -> LinearCombination<F>,
    ) -> crate::gr1cs::Result<()> {
        if !self.has_predicate(predicate_label) {
            return Err(SynthesisError::PredicateNotFound);
        }

        if self.should_construct_matrices() {
            let a = self.new_constraint_lc(a)?;
            let b = self.new_constraint_lc(b)?;
            let c = self.new_constraint_lc(c)?;

            let predicate = self
                .predicate_constraint_systems
                .get_mut(predicate_label)
                .unwrap();

            predicate.enforce_constraint([a, b, c])?;
        }

        #[cfg(feature = "std")]
        if let Some(traces) = self.predicate_traces.get_mut(predicate_label) {
            traces.push(ConstraintTrace::capture())
        }

        Ok(())
    }

    /// Enforce a constraint for a predicate with arity 4.
    pub fn enforce_constraint_arity_4(
        &mut self,
        predicate_label: &str,
        a: impl FnOnce() -> LinearCombination<F>,
        b: impl FnOnce() -> LinearCombination<F>,
        c: impl FnOnce() -> LinearCombination<F>,
        d: impl FnOnce() -> LinearCombination<F>,
    ) -> crate::gr1cs::Result<()> {
        if !self.has_predicate(predicate_label) {
            return Err(SynthesisError::PredicateNotFound);
        }

        if self.should_construct_matrices() {
            let a = self.new_constraint_lc(a)?;
            let b = self.new_constraint_lc(b)?;
            let c = self.new_constraint_lc(c)?;
            let d = self.new_constraint_lc(d)?;

            let predicate = self
                .predicate_constraint_systems
                .get_mut(predicate_label)
                .unwrap();

            predicate.enforce_constraint([a, b, c, d])?;
        }

        #[cfg(feature = "std")]
        if let Some(traces) = self.predicate_traces.get_mut(predicate_label) {
            traces.push(ConstraintTrace::capture())
        }

        Ok(())
    }

    /// Enforce a constraint for a predicate with arity 5.
    pub fn enforce_constraint_arity_5(
        &mut self,
        predicate_label: &str,
        a: impl FnOnce() -> LinearCombination<F>,
        b: impl FnOnce() -> LinearCombination<F>,
        c: impl FnOnce() -> LinearCombination<F>,
        d: impl FnOnce() -> LinearCombination<F>,
        e: impl FnOnce() -> LinearCombination<F>,
    ) -> crate::gr1cs::Result<()> {
        if !self.has_predicate(predicate_label) {
            return Err(SynthesisError::PredicateNotFound);
        }

        if self.should_construct_matrices() {
            let a = self.new_constraint_lc(a)?;
            let b = self.new_constraint_lc(b)?;
            let c = self.new_constraint_lc(c)?;
            let d = self.new_constraint_lc(d)?;
            let e = self.new_constraint_lc(e)?;

            let predicate = self
                .predicate_constraint_systems
                .get_mut(predicate_label)
                .unwrap();

            predicate.enforce_constraint([a, b, c, d, e])?;
        }

        #[cfg(feature = "std")]
        if let Some(traces) = self.predicate_traces.get_mut(predicate_label) {
            traces.push(ConstraintTrace::capture())
        }

        Ok(())
    }

    /// Enforce a constraint in the constraint system. It takes a
    /// predicate name and enforces a vector of linear combinations of the
    /// length of the arity of the predicate enforces the constraint.
    #[inline]
    pub fn enforce_r1cs_constraint(
        &mut self,
        a: impl FnOnce() -> LinearCombination<F>,
        b: impl FnOnce() -> LinearCombination<F>,
        c: impl FnOnce() -> LinearCombination<F>,
    ) -> crate::gr1cs::Result<()> {
        self.enforce_constraint_arity_3(R1CS_PREDICATE_LABEL, a, b, c)
    }

    /// Enforce a constraint in the constraint system. It takes a
    /// predicate name and enforces a vector of linear combinations of the
    /// length of the arity of the predicate enforces the constraint.
    #[inline]
    pub fn enforce_sr1cs_constraint(
        &mut self,
        a: impl FnOnce() -> LinearCombination<F>,
        b: impl FnOnce() -> LinearCombination<F>,
    ) -> crate::gr1cs::Result<()> {
        self.enforce_constraint_arity_2(SR1CS_PREDICATE_LABEL, a, b)
    }

    /// Add a new linear combination to the constraint system.
    /// This linear combination is to be used only for constraints, not for variables.
    #[inline]
    fn new_constraint_lc(
        &mut self,
        lc: impl FnOnce() -> LinearCombination<F>,
    ) -> crate::gr1cs::Result<Variable> {
        if self.should_construct_matrices() {
            self.new_lc_helper(lc)
        } else {
            self.new_lc_without_adding()
        }
    }

    /// Creates a new index for the linear combination without adding the concrete LC expression
    /// to the map.
    #[inline]
    fn new_lc_without_adding(&mut self) -> crate::gr1cs::Result<Variable> {
        let index = self.num_linear_combinations;
        self.num_linear_combinations += 1;
        Ok(Variable::SymbolicLc(index))
    }

    fn new_lc_add_helper(
        cur_num_lcs: &mut usize,
        lc_map: &mut LcMap<F>,
        should_generate_lc_assignments: bool,
        assignments: &mut Assignments<F>,
        lc: LinearCombination<F>,
    ) -> crate::gr1cs::Result<Variable> {
        match lc.0.as_slice() {
            // If the linear combination is empty, we return a symbolic LC with index 0.
            [] | [(_, Variable::Zero)] => Ok(Variable::SymbolicLc(0)),
            // If the linear combination is just another variable
            // with a coefficient of 1, we return the variable directly.
            [(c, var)] if c.is_one() => Ok(*var),
            // In all other cases, we create a new linear combination
            _ => {
                let index = *cur_num_lcs;
                lc_map.push(lc);
                *cur_num_lcs += 1;
                if should_generate_lc_assignments {
                    let value = assignments.eval_lc(index, lc_map).unwrap();
                    assignments.lc_assignment.push(value)
                }
                Ok(Variable::SymbolicLc(index))
            },
        }
    }

    /// Helper function to add a new linear combination to the constraint system.
    #[inline]
    fn new_lc_helper(
        &mut self,
        lc: impl FnOnce() -> LinearCombination<F>,
    ) -> crate::gr1cs::Result<Variable> {
        let should_push = self.should_construct_matrices() || self.should_generate_lc_assignments();
        let should_generate_lc_assignments = self.should_generate_lc_assignments();
        if should_push {
            let lc = lc();
            Self::new_lc_add_helper(
                &mut self.num_linear_combinations,
                &mut self.lc_map,
                should_generate_lc_assignments,
                &mut self.assignments,
                lc,
            )
        } else {
            self.new_lc_without_adding()
        }
    }

    /// Adds a new linear combination to the constraint system.
    #[inline]
    pub fn new_lc(
        &mut self,
        lc: impl FnOnce() -> LinearCombination<F>,
    ) -> crate::gr1cs::Result<Variable> {
        // Because this LC might be used to construct constraints,
        // we need to ensure that it is added to the lc_map whenever
        // `self.should_construct_matrices()` is true.
        // `self.new_lc_helper` will handle this.
        self.new_lc_helper(lc)
    }

    /// Set `self.mode` to `mode`.
    pub fn set_mode(&mut self, mode: SynthesisMode) {
        self.mode = mode;
    }

    /// Check whether `self.mode == SynthesisMode::Setup`.
    /// Returns true if
    /// 1. There is an underlying `ConstraintSystem`, and
    /// 2. It is in setup mode.
    pub fn is_in_setup_mode(&self) -> bool {
        self.mode == SynthesisMode::Setup
    }

    /// Check whether this constraint system aims to optimize weight,
    /// number of constraints, or neither.
    pub fn optimization_goal(&self) -> OptimizationGoal {
        self.optimization_goal
    }

    /// Check whether this constraint system is new, i.e., it is just created
    fn is_new(&self) -> bool {
        self.num_instance_variables == 1
            && self.num_witness_variables == 0
            && self.num_constraints() == 0
            && self.num_linear_combinations == 1
    }

    /// Specify whether this constraint system should aim to optimize weight,
    /// number of constraints, or neither.
    pub fn set_optimization_goal(&mut self, goal: OptimizationGoal) {
        assert!(self.is_new());
        self.optimization_goal = goal;
    }

    /// Check whether or not `self` will construct matrices.
    pub fn should_construct_matrices(&self) -> bool {
        match self.mode {
            SynthesisMode::Setup => true,
            SynthesisMode::Prove {
                construct_matrices, ..
            } => construct_matrices,
        }
    }

    /// Check whether or not `self` will construct matrices.
    pub fn should_generate_lc_assignments(&self) -> bool {
        match self.mode {
            SynthesisMode::Setup => false,
            SynthesisMode::Prove {
                generate_lc_assignments,
                ..
            } => generate_lc_assignments,
        }
    }

    /// Obtain a variable representing a new public instance input
    /// This function takes a closure, this closure returns `Result<F>`
    /// Internally, this function calls new_input_variable of the constraint
    /// system to which it's pointing
    #[inline]
    pub fn new_input_variable<Func>(&mut self, f: Func) -> crate::utils::Result<Variable>
    where
        Func: FnOnce() -> crate::utils::Result<F>,
    {
        let index = self.num_instance_variables;
        self.num_instance_variables += 1;

        if !self.is_in_setup_mode() {
            self.assignments.instance_assignment.push(f()?);
        }
        Ok(Variable::Instance(index))
    }

    /// Obtain a variable representing a new private witness input.
    #[inline]
    pub fn new_witness_variable<Func>(&mut self, f: Func) -> crate::utils::Result<Variable>
    where
        Func: FnOnce() -> crate::utils::Result<F>,
    {
        let index = self.num_witness_variables;
        self.num_witness_variables += 1;

        if !self.is_in_setup_mode() {
            self.assignments.witness_assignment.push(f()?);
        }
        Ok(Variable::Witness(index))
    }

    /// Register a predicate in the constraint system with a given label.
    pub fn register_predicate(
        &mut self,
        predicate_label: &str,
        predicate: PredicateConstraintSystem<F>,
    ) -> crate::utils::Result<()> {
        self.predicate_constraint_systems
            .insert(predicate_label.to_string(), predicate);
        #[cfg(feature = "std")]
        self.predicate_traces
            .insert(predicate_label.to_string(), Vec::new());
        Ok(())
    }

    /// Remove the predicate with the given label from the constraint system.
    pub fn remove_predicate(&mut self, predicate_label: &str) {
        self.predicate_constraint_systems.remove(predicate_label);
    }

    /// check if there is a predicate with the given label
    pub fn has_predicate(&self, predicate_label: &str) -> bool {
        self.predicate_constraint_systems
            .contains_key(predicate_label)
    }

    /// Obtain the assignment corresponding to the `Variable` `v`.
    pub fn assigned_value(&self, v: Variable) -> Option<F> {
        self.assignments.assigned_value(v)
    }

    /// If `self` is satisfied, outputs `Ok(true)`.
    /// If `self` is unsatisfied, outputs `Ok(false)`.
    /// If `self.is_in_setup_mode()` or if `self == None`, outputs `Err(())`.
    pub fn is_satisfied(&self) -> crate::utils::Result<bool> {
        self.which_is_unsatisfied().map(|s| s.is_none())
    }

    /// If `self` is satisfied, outputs `Ok(None)`.
    /// If `self` is unsatisfied, outputs `Some(s,i)`, where `s` is the label of
    /// the unsatisfied prediacate and  `i` is the index of
    /// the first unsatisfied constraint in that predicate.
    /// If `self.is_in_setup_mode()` or `self == None`, outputs `Err(())`.
    pub fn which_is_unsatisfied(&self) -> crate::utils::Result<Option<String>> {
        if self.is_in_setup_mode() {
            Err(SynthesisError::AssignmentMissing)
        } else {
            for (label, predicate) in self.predicate_constraint_systems.iter() {
                if let Some(unsatisfied_constraint) =
                    predicate.which_constraint_is_unsatisfied(self)
                {
                    #[cfg(feature = "std")]
                    let trace = self.predicate_traces[label][unsatisfied_constraint]
                        .as_ref()
                        .map_or_else(
                            || {
                                eprintln!("Constraint trace requires `ConstraintLayer`");
                                format!("{label} - {unsatisfied_constraint}")
                            },
                            |t| format!("{t}"),
                        )
                        .to_string();
                    #[cfg(not(feature = "std"))]
                    let trace = format!("{label} - {unsatisfied_constraint}");
                    return Ok(Some(trace));
                }
            }
            Ok(None)
        }
    }

    /// Finalize the constraint system (either by outlining or inlining,
    /// if an optimization goal is set).
    pub fn finalize(&mut self) {
        let timer_finalize = start_timer!(|| "Finalize GR1CS");
        let timer_inline_ouline_lcs = start_timer!(|| "Inline/Outline LCs");
        self.inline_all_lcs();
        end_timer!(timer_inline_ouline_lcs);
        // check if should outline instance or not
        let timer_instance_outlining = start_timer!(|| "Instance Outlining");
        if let Some(instance_outliner) = self.instance_outliner.take() {
            // Check if the predicate to be outlined is in the constraint system
            if self.has_predicate(&instance_outliner.pred_label) {
                // Outline the instances
                let _ = self.perform_instance_outlining(instance_outliner);
            }
        }
        end_timer!(timer_instance_outlining);
        end_timer!(timer_finalize);
    }

    /// Naively inlines symbolic linear combinations into the linear
    /// combinations that use them.
    ///
    /// Useful for standard pairing-based SNARKs where addition gates are
    /// cheap. For example, in the SNARKs such as [\[Groth16\]](https://eprint.iacr.org/2016/260) and
    /// [\[Groth-Maller17\]](https://eprint.iacr.org/2017/540), addition gates
    /// do not contribute to the size of the multi-scalar multiplication,
    /// which is the dominating cost.
    pub fn inline_all_lcs(&mut self) {
        if !self.should_construct_matrices() {
            return;
        }

        let any_used = self.any_lcs_used();
        if !any_used {
            return;
        }
        let old_lc_map = core::mem::take(&mut self.lc_map);
        let mut inlined_lcs = LcMap::with_capacity(old_lc_map.len());

        let mut out = LinearCombination(Vec::with_capacity(10));
        for lc in old_lc_map.iter() {
            for (coeff, var) in lc {
                if let Some(lc_index) = var.get_lc_index() {
                    // Must already be transformed — guaranteed by ordering.
                    let inlined = &inlined_lcs[lc_index];

                    if coeff.is_one() {
                        out.extend_from_slice(inlined);
                    } else {
                        out.extend(
                            inlined
                                .iter()
                                .filter(|(_, v)| !v.is_zero())
                                .map(|(c, v)| (*coeff * c, *v)),
                        );
                    }
                } else {
                    out.push((*coeff, *var));
                }
            }
            out.compactify();
            inlined_lcs.push(out.clone());
            out.0.clear();
        }
        self.lc_map = inlined_lcs;
    }

    /// Returns whether any linear combinations are used
    /// by other LCs.
    fn any_lcs_used(&self) -> bool {
        cfg_iter!(self.lc_map).any(|lc| lc.iter().any(|(_, v)| v.is_lc()))
    }

    /// Get the matrices corresponding to the predicates.and the
    /// corresponding set of matrices
    pub fn to_matrices(&self) -> crate::gr1cs::Result<BTreeMap<Label, Vec<Matrix<F>>>> {
        let mut matrices = BTreeMap::new();
        for (label, predicate) in self.predicate_constraint_systems.iter() {
            matrices.insert(label.clone(), predicate.to_matrices(self));
        }
        Ok(matrices)
    }

    /// Get the linear combination corresponding to the given `lc_index`.
    pub fn get_lc(&self, var: Variable) -> LinearCombination<F> {
        if var.is_zero() {
            LinearCombination::zero()
        } else if var.is_lc() {
            let lc_index = var.index().unwrap();
            LinearCombination(self.lc_map[lc_index].to_vec())
        } else {
            LinearCombination::from(var)
        }
    }

    /// Given a linear combination, create a row in the matrix
    #[inline]
    pub(crate) fn make_row(&self, l: LinearCombination<F>) -> Vec<(F, usize)> {
        let num_input = self.num_instance_variables();
        l.0.into_iter()
            .filter_map(|(coeff, var)| {
                if coeff.is_zero() || var.is_zero() {
                    None
                } else {
                    let index = var.get_variable_index(num_input);
                    Some((coeff, index.unwrap()))
                }
            })
            .collect()
    }

    /// Sets the flag for outlining the instances
    pub(crate) fn set_instance_outliner(&mut self, instance_outliner: InstanceOutliner<F>) {
        self.instance_outliner = Some(instance_outliner);
    }

    /// Returns the flag for outlining the instances, This is by default set
    /// to false
    pub(crate) fn should_outline_instances(&self) -> bool {
        self.instance_outliner.is_some()
    }

    /// Outlines the instances in the constraint system
    /// This function creates a new witness variable for each instance
    /// variable and uses these witness variables in the constraints
    /// instead of instance variables. This technique is useful for
    /// verifier succinctness in some SNARKs like Garuda, Pari and
    /// PolyMath After the function call, The instances are only
    /// used in the `c` matrix of r1cs
    pub fn perform_instance_outlining(
        &mut self,
        outliner: InstanceOutliner<F>,
    ) -> crate::gr1cs::Result<()> {
        // First build a map from instance variables to witness variables
        let mut instance_to_witness_map = Vec::<Variable>::new();
        // Initialize the map with the one variable, this is done manually because we
        // certainely need this variable and it might not show up in the lc_map
        let one_witness_var = self.new_witness_variable(|| Ok(F::ONE))?;
        instance_to_witness_map.push(one_witness_var);
        let instance_assignment = &self.assignments.instance_assignment.clone();
        // Skip the first one because that is the ONE variable.
        for i in 1..self.num_instance_variables {
            let value = instance_assignment.get(i).copied();
            let witness_var =
                self.new_witness_variable(|| value.ok_or(SynthesisError::AssignmentMissing))?;
            instance_to_witness_map.push(witness_var);
        }

        // Now, Go over all the linear combinations and create a new witness for each
        // instance variable you see
        cfg_iter_mut!(self.lc_map).for_each(|lc| {
            for (_, var) in lc.iter_mut() {
                if var.is_instance() {
                    *var = instance_to_witness_map[var.index().unwrap()];
                } else if var.is_one() {
                    *var = one_witness_var;
                }
            }
        });
        (outliner.func)(self, &instance_to_witness_map)?;
        Ok(())
    }
}

#[derive(Clone, Debug, Default)]
pub struct LcMap<F: Field>(Vec<LinearCombination<F>>);

impl<F: Field> LcMap<F> {
    #[inline(always)]
    pub fn new() -> Self {
        Self(Vec::with_capacity(10))
    }

    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        Self(Vec::with_capacity(capacity))
    }

    #[inline(always)]
    pub fn push(&mut self, v: impl IntoIterator<Item = (F, Variable)>) {
        self.0.push(LinearCombination(v.into_iter().collect()));
    }

    #[inline(always)]
    pub fn push_by_ref(&mut self, v: &[(F, Variable)]) {
        self.0.push(LinearCombination(v.to_vec()));
    }

    #[inline(always)]
    pub fn iter(&self) -> impl Iterator<Item = &[(F, Variable)]> {
        self.0.iter().map(|lc| lc.0.as_slice())
    }

    #[cfg(feature = "parallel")]
    #[inline(always)]
    pub fn par_iter(&self) -> impl ParallelIterator<Item = &[(F, Variable)]> {
        self.0.par_iter().map(|lc| lc.0.as_slice())
    }

    #[inline(always)]
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut [(F, Variable)]> {
        self.0.iter_mut().map(|lc| lc.0.as_mut_slice())
    }

    #[cfg(feature = "parallel")]
    #[inline(always)]
    pub fn par_iter_mut(&mut self) -> impl ParallelIterator<Item = &mut [(F, Variable)]> {
        self.0.par_iter_mut().map(|lc| lc.0.as_mut_slice())
    }

    #[inline(always)]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline(always)]
    #[allow(unsafe_code)]
    pub fn get(&self, idx: usize) -> Option<&[(F, Variable)]> {
        self.0.get(idx).map(|lc| lc.0.as_slice())
    }
}

impl<F: Field> core::ops::Index<usize> for LcMap<F> {
    type Output = [(F, Variable)];

    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}
