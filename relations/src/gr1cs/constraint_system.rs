//! This module contains the implementation of the `ConstraintSystem` struct.
//! a constraint system contains multiple predicate constraint systems,
//! each of which enforce have seperate predicates and constraints. For more infomation about the terminology and the structure of the constraint system, refer to section 3.3 of https://eprint.iacr.org/2024/1245

use super::{
    predicate::{
        polynomial_constraint::R1CS_PREDICATE_LABEL, PredicateConstraintSystem, Predicate,
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
use hashbrown::HashMap;
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
    outline_instances: bool,

    /// Assignments to the public input variables. This is empty if `self.mode
    /// == SynthesisMode::Setup`.
    instance_assignment: Vec<F>,

    /// Assignments to the private input variables. This is empty if `self.mode
    /// == SynthesisMode::Setup`.
    witness_assignment: Vec<F>,

    /// Map for gadgets to cache computation results.
    pub cache_map: Rc<RefCell<BTreeMap<TypeId, Box<dyn Any>>>>,

    /// A data structure to store the linear combinations. We use map because
    /// it's easier to inline and outline the linear combinations.
    lc_map: Vec<LinearCombination<F>>,

    /// A map from the the predicate labels to the predicates
    predicate_constraint_systems: BTreeMap<Label, PredicateConstraintSystem<F>>,

    /// A cache for the linear combination assignments. It shows evaluation
    /// result of each linear combination
    lc_assignment_cache: Rc<RefCell<HashMap<LcIndex, F>>>,

    /// data structure to store the traces for each predicate
    #[cfg(feature = "std")]
    pub predicate_traces: BTreeMap<Label, Vec<Option<ConstraintTrace>>>,
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
            outline_instances: false,
            predicate_constraint_systems: BTreeMap::new(),
            instance_assignment: vec![F::one()],
            witness_assignment: Vec::new(),
            cache_map: Rc::new(RefCell::new(BTreeMap::new())),
            lc_map: Vec::new(),
            lc_assignment_cache: Rc::new(RefCell::new(HashMap::new())),
            mode: SynthesisMode::Prove {
                construct_matrices: true,
            },
            optimization_goal: OptimizationGoal::Constraints,
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
    pub fn get_all_predicates_num_constraints(&self) -> HashMap<Label, usize> {
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
    pub fn get_all_predicate_arities(&self) -> HashMap<Label, usize> {
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
        Ok(&self.instance_assignment)
    }

    /// Returns the assignment to the private input variables of the constraint
    pub fn witness_assignment(&self) -> crate::gr1cs::Result<&[F]> {
        if self.is_in_setup_mode() {
            return Err(SynthesisError::AssignmentMissing);
        }
        Ok(&self.witness_assignment)
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
        lc_vec: impl IntoIterator<Item = LinearCombination<F>>,
    ) -> crate::gr1cs::Result<()> {
        if !self.has_predicate(predicate_label) {
            return Err(SynthesisError::PredicateNotFound);
        }
        if self.should_construct_matrices() {
            // TODO: Find a way to get rid of the following collect, barrier1: we need the
            // size, possible to use ExactSizeIterator, barrier2: We will need lifetimes
            // which leads to having two &mut refs to self

            let lc_indices = lc_vec.into_iter().map(|lc| {
                let var = {
                    let index = LcIndex(self.num_linear_combinations);
                    self.lc_map.push(lc);
                    self.num_linear_combinations += 1;
                    Variable::SymbolicLc(index)
                };
                match var {
                    Variable::SymbolicLc(index) => index,
                    _ => panic!("Unexpected variable type"),
                }
            });
            let predicate = self
                .predicate_constraint_systems
                .get_mut(predicate_label)
                .unwrap();

            predicate.enforce_constraint(lc_indices)?;
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
        self.lc_map.push(lc);
        self.num_linear_combinations += 1;
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
            SynthesisMode::Prove { construct_matrices } => construct_matrices,
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
            self.instance_assignment.push(f()?);
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
            self.witness_assignment.push(f()?);
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

    /// check if there is a predicate with the given label
    pub fn has_predicate(&self, predicate_label: &str) -> bool {
        self.predicate_constraint_systems
            .contains_key(predicate_label)
    }

    /// Obtain the assignment corresponding to the `Variable` `v`.
    pub fn assigned_value(&self, v: Variable) -> Option<F> {
        match v {
            Variable::One => Some(F::one()),
            Variable::Zero => Some(F::zero()),
            Variable::Witness(idx) => self.witness_assignment.get(idx).copied(),
            Variable::Instance(idx) => self.instance_assignment.get(idx).copied(),
            Variable::SymbolicLc(idx) => {
                let value = self.lc_assignment_cache.borrow().get(&idx).copied();
                if value.is_some() {
                    value
                } else {
                    let value = self.eval_lc(idx)?;
                    self.lc_assignment_cache.borrow_mut().insert(idx, value);
                    Some(value)
                }
            },
        }
    }

    /// Evaluate the linear combination `lc` with the assigned values and return
    /// the result.
    fn eval_lc(&self, lc: LcIndex) -> Option<F> {
        let lc = &self.lc_map[lc.0];
        let mut acc = F::zero();
        for (coeff, var) in lc.iter() {
            acc += *coeff * self.assigned_value(*var)?;
        }
        Some(acc)
    }

    /// If `self` is satisfied, outputs `Ok(true)`.
    /// If `self` is unsatisfied, outputs `Ok(false)`.
    /// If `self.is_in_setup_mode()` or if `self == None`, outputs `Err(())`.
    pub fn is_satisfied(&self) -> crate::utils::Result<bool> {
        self.which_predicate_is_unsatisfied().map(|s| s.is_none())
    }

    /// If `self` is satisfied, outputs `Ok(None)`.
    /// If `self` is unsatisfied, outputs `Some(s,i)`, where `s` is the label of
    /// the unsatisfied prediacate and  `i` is the index of
    /// the first unsatisfied constraint in that predicate.
    /// If `self.is_in_setup_mode()` or `self == None`, outputs `Err(())`.
    pub fn which_predicate_is_unsatisfied(&self) -> crate::utils::Result<Option<String>> {
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
        match self.optimization_goal {
            OptimizationGoal::None => self.inline_all_lcs(),
            OptimizationGoal::Constraints => self.inline_all_lcs(),
            OptimizationGoal::Weight => self.outline_lcs(),
        };
        if self.outline_instances {
            let _ = self.do_outline_instances();
        }
    }

    /// Naively inlines symbolic linear combinations into the linear
    /// combinations that use them.
    ///
    /// Useful for standard pairing-based SNARKs where addition gates are cheap.
    /// For example, in the SNARKs such as [\[Groth16\]](https://eprint.iacr.org/2016/260) and
    /// [\[Groth-Maller17\]](https://eprint.iacr.org/2017/540), addition gates
    /// do not contribute to the size of the multi-scalar multiplication, which
    /// is the dominating cost.
    pub fn inline_all_lcs(&mut self) {
        // Only inline when a matrix representing R1CS is needed.
        if !self.should_construct_matrices() {
            return;
        }

        // A dummy closure is used, which means that
        // - it does not modify the inlined LC.
        // - it does not add new witness variables.
        self.transform_lc_map(&mut |_, _, _| (0, None));
    }

    /// If a `SymbolicLc` is used in more than one location and has sufficient
    /// length, this method makes a new variable for that `SymbolicLc`, adds
    /// a constraint ensuring the equality of the variable and the linear
    /// combination, and then uses that variable in every location the
    /// `SymbolicLc` is used.
    ///
    /// Useful for SNARKs like [\[Marlin\]](https://eprint.iacr.org/2019/1047) or
    /// [\[Fractal\]](https://eprint.iacr.org/2019/1076), where addition gates
    /// are not cheap.
    fn outline_lcs(&mut self) {
        // Only inline when a matrix representing R1CS is needed.
        if !self.should_construct_matrices() {
            return;
        }

        // Store information about new witness variables created
        // for outlining. New constraints will be added after the
        // transformation of the LC map.
        let mut new_witness_linear_combinations = Vec::new();
        let mut new_witness_indices = Vec::new();

        // It goes through all the LCs in the map, starting from
        // the early ones, and decides whether or not to dedicate a witness
        // variable for this LC.
        //
        // If true, the LC is replaced with 1 * this witness variable.
        // Otherwise, the LC is inlined.
        //
        // Each iteration first updates the LC according to outlinings in prior
        // iterations, and then sees if it should be outlined, and if so adds
        // the outlining to the map.
        //
        self.transform_lc_map(&mut |cs, num_times_used, inlined_lc| {
            let mut should_dedicate_a_witness_variable = false;
            let mut new_witness_index = None;
            let mut new_witness_assignment = Vec::new();

            // Check if it is worthwhile to dedicate a witness variable.
            let this_used_times = num_times_used + 1;
            let this_len = inlined_lc.len();

            // Cost with no outlining = `lc_len * number of usages`
            // Cost with outlining is one constraint for `(lc_len) * 1 = {new variable}` and
            // using that single new variable in each of the prior usages.
            // This has total cost `number_of_usages + lc_len + 2`
            if this_used_times * this_len > this_used_times + 2 + this_len {
                should_dedicate_a_witness_variable = true;
            }

            // If it is worthwhile to dedicate a witness variable,
            if should_dedicate_a_witness_variable {
                // Add a new witness (the value of the linear combination).
                // This part follows the same logic of `new_witness_variable`.
                let witness_index = cs.num_witness_variables;
                new_witness_index = Some(witness_index);

                // Compute the witness assignment.
                if !cs.is_in_setup_mode() {
                    let mut acc = F::zero();
                    for (coeff, var) in inlined_lc.iter() {
                        acc += *coeff * cs.assigned_value(*var).unwrap();
                    }
                    new_witness_assignment.push(acc);
                }

                // Add a new constraint for this new witness.
                new_witness_linear_combinations.push(inlined_lc.clone());
                new_witness_indices.push(witness_index);

                // Replace the linear combination with (1 * this new witness).
                *inlined_lc = LinearCombination::from(Variable::Witness(witness_index));
            }
            // Otherwise, the LC remains unchanged.

            // Return information about new witness variables.
            if new_witness_index.is_some() {
                (1, Some(new_witness_assignment))
            } else {
                (0, None)
            }
        });

        // Add the constraints for the newly added witness variables.
        for (new_witness_linear_combination, new_witness_variable) in
            new_witness_linear_combinations
                .iter()
                .zip(new_witness_indices.iter())
        {
            let r1cs_constraint: Vec<LinearCombination<F>> = vec![
                new_witness_linear_combination.clone(),
                LinearCombination::from(Variable::one()),
                LinearCombination::from(Variable::Witness(*new_witness_variable)),
            ];
            // Add a new constraint
            self.enforce_constraint(R1CS_PREDICATE_LABEL, r1cs_constraint)
                .unwrap();
        }
    }

    /// Transform the map of linear combinations.
    /// Specifically, allow the creation of additional witness assignments.
    ///
    /// This method is used as a subroutine of `inline_all_lcs` and
    /// `outline_lcs`.
    ///
    /// The transformer function is given a references of this constraint system
    /// (&self), number of times used, and a mutable reference of the linear
    /// combination to be transformed.
    ///
    /// The transformer function returns the number of new witness variables
    /// needed and a vector of new witness assignments (if not in the setup
    /// mode).
    pub fn transform_lc_map(
        &mut self,
        transformer: &mut dyn FnMut(
            &ConstraintSystem<F>,
            usize,
            &mut LinearCombination<F>,
        ) -> (usize, Option<Vec<F>>),
    ) {
        // `transformed_lc_map` stores the transformed linear combinations.
        let mut transformed_lc_map = BTreeMap::<LcIndex, LinearCombination<F>>::new();
        let mut num_times_used = self.lc_num_times_used(false);

        // This loop goes through all the LCs in the map, starting from
        // the early ones. The transformer function is applied to the
        // inlined LC, where new witness variables can be created.
        for (index, lc) in self.lc_map.iter().enumerate() {
            let mut transformed_lc = LinearCombination::new();

            // Inline the LC, unwrapping symbolic LCs that may constitute it,
            // and updating them according to transformations in prior iterations.
            for &(coeff, var) in lc.iter() {
                if var.is_lc() {
                    let lc_index = var.get_lc_index().expect("should be lc");

                    // If `var` is a `SymbolicLc`, fetch the corresponding
                    // inlined LC, and substitute it in.
                    //
                    // We have the guarantee that `lc_index` must exist in
                    // `new_lc_map` since a LC can only depend on other
                    // LCs with lower indices, which we have transformed.
                    //
                    let lc = transformed_lc_map
                        .get(&lc_index)
                        .expect("should be inlined");
                    transformed_lc.extend((lc * coeff).0.into_iter());

                    // Delete linear combinations that are no longer used.
                    //
                    // Deletion is safe for both outlining and inlining:
                    // * Inlining: the LC is substituted directly into all use sites, and so once it
                    //   is fully inlined, it is redundant.
                    //
                    // * Outlining: the LC is associated with a new variable `w`, and a new
                    //   constraint of the form `lc_data * 1 = w`, where `lc_data` is the actual
                    //   data in the linear combination. Furthermore, we replace its entry in
                    //   `new_lc_map` with `(1, w)`. Once `w` is fully inlined, then we can delete
                    //   the entry from `new_lc_map`
                    //
                    num_times_used[lc_index.0] -= 1;
                    if num_times_used[lc_index.0] == 0 {
                        // This lc is not used any more, so remove it.
                        transformed_lc_map.remove(&lc_index);
                    }
                } else {
                    // Otherwise, it's a concrete variable and so we
                    // substitute it in directly.
                    transformed_lc.push((coeff, var));
                }
            }
            transformed_lc.compactify();

            // Call the transformer function.
            let (num_new_witness_variables, new_witness_assignments) =
                transformer(self, num_times_used[index], &mut transformed_lc);

            // Insert the transformed LC.
            transformed_lc_map.insert(LcIndex(index), transformed_lc);

            // Update the witness counter.
            self.num_witness_variables += num_new_witness_variables;

            // Supply additional witness assignments if not in the
            // setup mode and if new witness variables are created.
            if !self.is_in_setup_mode() && num_new_witness_variables > 0 {
                assert!(new_witness_assignments.is_some());
                if let Some(new_witness_assignments) = new_witness_assignments {
                    assert_eq!(new_witness_assignments.len(), num_new_witness_variables);
                    self.witness_assignment
                        .extend_from_slice(&new_witness_assignments);
                }
            }
        }
        // Replace the LC map.
        self.lc_map = transformed_lc_map.into_values().collect();
    }

    /// Count the number of times each linear combination is used.
    fn lc_num_times_used(&self, count_sinks: bool) -> Vec<usize> {
        let mut num_times_used = vec![0; self.lc_map.len()];

        // Iterate over every lc in constraint system
        for (index, lc) in self.lc_map.iter().enumerate() {
            num_times_used[index] += count_sinks as usize;

            // Increment the counter for each lc that this lc has a direct dependency on.
            for &(_, var) in lc.iter() {
                if var.is_lc() {
                    let lc_index = var.get_lc_index().expect("should be lc");
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
    /// TODO: This function should return a reference to the linear combination
    /// and not clone it.
    pub fn get_lc(&self, lc_index: LcIndex) -> crate::gr1cs::Result<LinearCombination<F>> {
        dbg!(self.lc_map.len());
        self.lc_map
            .get(lc_index.0)
            .cloned()
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
    pub(crate) fn outline_instances(&mut self) {
        self.outline_instances = true;
    }

    /// Returns the flag for outlining the instances, This is by default set to
    /// false
    pub(crate) fn should_outline_instances(&self) -> bool {
        self.outline_instances
    }

    /// Outlines the instances in the constraint system
    /// This function creates a new witness variable for each instance variable
    /// and uses these witness variables in the constraints instead of instance
    /// variables. This technique is useful for verifier succinctness in some
    /// SNARKs like Garuda, Pari and PolyMath
    pub(crate) fn do_outline_instances(&mut self) -> crate::gr1cs::Result<()> {
        let mut instance_to_witness_map = BTreeMap::<Variable, Variable>::new();
        instance_to_witness_map
            .insert(Variable::One, Variable::Witness(self.num_witness_variables));
        self.num_witness_variables += 1;

        for lc in self.lc_map.iter() {
            for (_, var) in lc.iter() {
                if var.is_instance() {
                    let _witness = instance_to_witness_map
                        .entry(*var)
                        .or_insert(Variable::Witness(self.num_witness_variables));
                    self.num_witness_variables += 1;
                }
            }
        }

        for lc in self.lc_map.iter_mut() {
            for (_, var) in lc.iter_mut() {
                if var.is_instance() {
                    *var = instance_to_witness_map[var];
                }
            }
        }

        if !self.is_in_setup_mode() {
            self.witness_assignment.resize(
                self.witness_assignment.len() + instance_to_witness_map.len(),
                F::zero(),
            );
            for (instance, witness) in instance_to_witness_map.iter() {
                let instance_value = self.assigned_value(*instance).unwrap();

                let witness_index = match witness {
                    Variable::Witness(index) => *index,
                    _ => unreachable!(),
                };

                self.witness_assignment[witness_index] = instance_value;
            }
        }

        let one_witt = instance_to_witness_map.get(&Variable::One).unwrap();
        for (instance, witness) in instance_to_witness_map.iter() {
            let r1cs_constraint = vec![
                LinearCombination::from(*instance),
                LinearCombination::from(*one_witt),
                LinearCombination::from(*witness),
            ];
            self.enforce_constraint(R1CS_PREDICATE_LABEL, r1cs_constraint)?;
        }

        Ok(())
    }
}
