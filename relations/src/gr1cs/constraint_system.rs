//! This module contains the implementation of the `ConstraintSystem` struct.
//! a constraint system contains multiple predicate constraint systems,
//! each of which enforce have separate predicates and constraints. For more information about the terminology and the structure of the constraint system, refer to section 3.3 of https://eprint.iacr.org/2024/1245

use super::{
    instance_outliner::InstanceOutliner,
    predicate::{
        polynomial_constraint::R1CS_PREDICATE_LABEL, Predicate, PredicateConstraintSystem,
    },
    ConstraintSystemRef, Label, OptimizationGoal, SynthesisMode,
};
#[cfg(feature = "std")]
use crate::gr1cs::ConstraintTrace;
use crate::gr1cs::{LcIndex, LinearCombination, Matrix, SynthesisError, Variable};
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
    num_instance_variables: usize,

    /// The number of variables that are "private inputs" to the constraint
    /// system.
    num_witness_variables: usize,

    /// The number of linear combinations
    num_linear_combinations: usize,

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
    lc_map: Vec<Option<LinearCombination<F>>>,

    /// A map from the the predicate labels to the predicates
    predicate_constraint_systems: BTreeMap<Label, PredicateConstraintSystem<F>>,

    /// data structure to store the traces for each predicate
    #[cfg(feature = "std")]
    pub predicate_traces: BTreeMap<Label, Vec<Option<ConstraintTrace>>>,
    
    /// A mapping between unoptimized and optimized witnesses.
    /// This is used to avoid regenerating constraints during proving.
    /// The key is the index of the unoptimized witness, and the value is the index of the optimized witness.
    witness_mapping: BTreeMap<usize, usize>,
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
        match v {
            Variable::One => Some(F::one()),
            Variable::Zero => Some(F::zero()),
            Variable::Witness(idx) => self.witness_assignment.get(idx).copied(),
            Variable::Instance(idx) => self.instance_assignment.get(idx).copied(),
            Variable::SymbolicLc(idx) => self.lc_assignment.get(idx.0).copied(),
        }
    }

    /// Evaluate the linear combination `lc` with the assigned values and return
    /// the result.
    #[inline]
    fn eval_lc(&self, lc: LcIndex, lc_map: &[Option<LinearCombination<F>>]) -> Option<F> {
        let acc = lc_map[lc.0]
            .as_ref()
            .unwrap()
            .0
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
        let mut cs = Self {
            num_instance_variables: 1,
            num_witness_variables: 0,
            num_linear_combinations: 0,
            instance_outliner: None,
            predicate_constraint_systems: BTreeMap::new(),
            assignments: Assignments {
                instance_assignment: vec![F::one()],
                witness_assignment: Vec::new(),
                lc_assignment: Vec::new(),
            },
            cache_map: Rc::new(RefCell::new(BTreeMap::new())),
            lc_map: Vec::new(),
            mode: SynthesisMode::Prove {
                construct_matrices: true,
                generate_lc_assignments: true,
            },
            optimization_goal: OptimizationGoal::None,
            #[cfg(feature = "std")]
            predicate_traces: BTreeMap::new(),
            witness_mapping: BTreeMap::new(),
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
        lcs: impl IntoIterator<Item = LinearCombination<F>, IntoIter: ExactSizeIterator>,
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
                let index = LcIndex(*num_lcs);
                lc_map.push(Some(lc));
                *num_lcs += 1;
                if should_generate_lc_assignments {
                    let value = assignments.eval_lc(index, lc_map).unwrap();
                    assignments.lc_assignment.push(value);
                }
                index
            });

            let predicate = self
                .predicate_constraint_systems
                .get_mut(predicate_label)
                .unwrap();

            predicate.enforce_constraint(lc_indices)?;
        } else if self.should_generate_lc_assignments() {
            // We're in proving mode but not constructing matrices
            // We need to evaluate the linear combinations to add them to the assignment
            let assignments = &mut self.assignments;
            let lc_map = &mut self.lc_map;
            let num_lcs = &mut self.num_linear_combinations;
            let witness_mapping = &self.witness_mapping;
            
            for lc in lcs {
                // Replace any instance variables with their mapped witness variables
                let mut mapped_lc = LinearCombination(Vec::new());
                for (coeff, var) in lc.0 {
                    let mapped_var = match var {
                        Variable::Instance(i) => {
                            if let Some(&witness_idx) = witness_mapping.get(&i) {
                                Variable::Witness(witness_idx)
                            } else {
                                var
                            }
                        },
                        _ => var,
                    };
                    mapped_lc += (coeff, mapped_var);
                }
                
                // Create a temporary index for evaluation
                let temp_index = LcIndex(*num_lcs);
                lc_map.push(Some(mapped_lc));
                *num_lcs += 1;
                
                // Evaluate and store the result
                let value = assignments.eval_lc(temp_index, lc_map).unwrap();
                assignments.lc_assignment.push(value);
            }
        }

        #[cfg(feature = "std")]
        {
            let trace = ConstraintTrace::capture();
            match self.predicate_traces.get_mut(predicate_label) {
                Some(traces) => traces.push(trace),
                None => {
                    eprintln!("Constraint trace requires adding the predicate constraint trace")
                },
            }
        }

        Ok(())
    }

    /// Adds a new linear combination to the constraint system.
    /// If the linear combination is already in the map, return the
    /// corresponding index (using bimap)
    #[inline]
    pub fn new_lc(&mut self, lc: LinearCombination<F>) -> crate::gr1cs::Result<Variable> {
        // Note: update also enforce_constraint if you change this logic.
        let index = LcIndex(self.num_linear_combinations);
        self.lc_map.push(Some(lc));
        self.num_linear_combinations += 1;
        if self.should_generate_lc_assignments() {
            let value = self
                .eval_lc(index)
                .ok_or(SynthesisError::AssignmentMissing)?;
            self.assignments.lc_assignment.push(value)
        }
        Ok(Variable::SymbolicLc(index))
    }

    /// Set `self.mode` to `mode`.
    pub fn set_mode(&mut self, mode: SynthesisMode) {
        self.mode = mode;
    }

    /// Check whether `self.mode == SynthesisMode::Setup`.
    /// Returns true if 1- There is an underlying ConstraintSystem and
    /// 2- It's in setup mode
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
            && self.num_linear_combinations == 0
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
        {
            self.predicate_traces
                .insert(predicate_label.to_string(), Vec::new());
        }
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

    /// Evaluate the linear combination `lc` with the assigned values and return
    /// the result.
    fn eval_lc(&self, lc: LcIndex) -> Option<F> {
        self.assignments.eval_lc(lc, &self.lc_map)
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
                    let trace;
                    #[cfg(feature = "std")]
                    {
                        trace = self.predicate_traces[label][unsatisfied_constraint]
                            .as_ref()
                            .map_or_else(
                                || {
                                    eprintln!(
                                        "Constraint trace requires enabling `ConstraintLayer`"
                                    );
                                    format!("{label} - {unsatisfied_constraint}")
                                },
                                |t| format!("{t}"),
                            )
                            .to_string();
                    }
                    #[cfg(not(feature = "std"))]
                    {
                        trace = format!("{label} - {unsatisfied_constraint}");
                    }
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

        let mut num_times_used = self.lc_num_times_used();
        let old_lc_map = core::mem::take(&mut self.lc_map);
        let mut inlined_lcs: Vec<Option<LinearCombination<_>>> =
            Vec::with_capacity(old_lc_map.len());

        for lc_opt in old_lc_map.into_iter() {
            let lc = lc_opt.expect("LC should never be None");
            let mut out = LinearCombination(Vec::with_capacity(lc.len()));

            for (coeff, var) in lc.0 {
                if let Some(lc_index) = var.get_lc_index() {
                    // Must already be transformed — guaranteed by ordering.
                    let inlined = inlined_lcs[lc_index.0]
                        .as_ref()
                        .expect("inlined LC must exist");

                    if coeff.is_one() {
                        out.extend_from_slice(&inlined.0);
                    } else {
                        out.extend(inlined.iter().map(|(c, v)| (coeff * c, *v)));
                    }
                    // Decrement usage and prune if no longer needed
                    num_times_used[lc_index.0] -= 1;
                    if num_times_used[lc_index.0] == 0 {
                        inlined_lcs[lc_index.0] = None;
                    }
                } else {
                    out.push((coeff, var));
                }
            }
            out.compactify();
            inlined_lcs.push(Some(out));
        }
        self.lc_map = inlined_lcs;
    }

    /// Count the number of times each linear combination is used.
    fn lc_num_times_used(&self) -> Vec<usize> {
        let mut num_times_used = vec![0; self.lc_map.len()];

        // Iterate over every lc in constraint system
        for lc in &self.lc_map {
            // Increment the counter for each lc that this lc has a direct dependency on.
            for &(_, var) in lc.as_ref().unwrap().iter() {
                if let Some(lc_index) = var.get_lc_index() {
                    num_times_used[lc_index.0] += 1;
                }
            }
        }
        num_times_used
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
    pub fn get_lc(&self, lc_index: LcIndex) -> crate::gr1cs::Result<&LinearCombination<F>> {
        self.lc_map
            .get(lc_index.0)
            .map(|e| e.as_ref())
            .flatten()
            .ok_or(SynthesisError::LcNotFound(lc_index))
    }

    /// Given a linear combination, create a row in the matrix
    #[inline]
    pub(crate) fn make_row(&self, l: &LinearCombination<F>) -> Vec<(F, usize)> {
        let num_input = self.num_instance_variables();
        l.0.iter()
            .filter_map(|(coeff, var)| {
                if coeff.is_zero() {
                    None
                } else {
                    Some((
                        *coeff,
                        var.get_index_unchecked(num_input).expect("no symbolic LCs"),
                    ))
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
            
            // Store the mapping from instance to witness variable
            if let Variable::Witness(witness_idx) = witness_var {
                self.witness_mapping.insert(i, witness_idx);
            }
        }

        // Now, Go over all the linear combinations and create a new witness for each
        // instance variable you see
        cfg_iter_mut!(self.lc_map)
            .filter(|lc| lc.is_some())
            .for_each(|lc| {
                for (_, var) in lc.as_mut().unwrap().iter_mut() {
                    if let Variable::Instance(i) = var {
                        *var = instance_to_witness_map[*i];
                    } else if let Variable::One = var {
                        *var = one_witness_var;
                    }
                }
            });
        (outliner.func)(self, &instance_to_witness_map)?;
        Ok(())
    }

    /// Set the witness mapping from an external source.
    /// This is used to avoid regenerating constraints during proving.
    /// The mapping should be between unoptimized instance variables and optimized witness variables.
    pub fn set_witness_mapping(&mut self, mapping: BTreeMap<usize, usize>) {
        self.witness_mapping = mapping;
    }

    /// Get the current witness mapping.
    /// This can be used to store the mapping in an index file.
    pub fn get_witness_mapping(&self) -> &BTreeMap<usize, usize> {
        &self.witness_mapping
    }
}
