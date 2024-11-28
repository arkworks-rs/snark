//! This module contains the implementation of a sharef reference to
//! `ConstraintSystem` struct. Refer to the `ConstraintSystem` struct for the
//! inner struct. Most of the functions of `ConstraintSystemRef` are just
//! wrappers around the functions of `ConstraintSystem`.

use ark_std::{collections::BTreeMap, rc::Weak};
use core::cell::RefCell;

use super::{
    constraint_system::ConstraintSystem,
    local_predicate::{polynomial_constraint::R1CS_PREDICATE_LABEL, PredicateConstraintSystem},
    Label, LcIndex, LinearCombination, Matrix, OptimizationGoal, SynthesisError, SynthesisMode,
    Variable,
};
use ark_ff::Field;
use ark_std::{rc::Rc, string::String, vec::Vec};

/// A shared reference to a constraint system that can be stored in high level
/// variables.
#[derive(Debug, Clone)]
pub enum ConstraintSystemRef<F: Field> {
    /// Represents the case where we *don't* need to allocate variables or
    /// enforce constraints. Encountered when operating over constant
    /// values.
    None,
    /// Represents the case where we *do* allocate variables or enforce
    /// constraints.
    CS(Rc<RefCell<ConstraintSystem<F>>>),
}

impl<F: Field> PartialEq for ConstraintSystemRef<F> {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::None, Self::None) => true,
            (..) => false,
        }
    }
}

impl<F: Field> Eq for ConstraintSystemRef<F> {}

impl<F: Field> ConstraintSystemRef<F> {
    /// Construct a `ConstraintSystemRef` from a `ConstraintSystem`.
    #[inline]
    pub fn new(inner: ConstraintSystem<F>) -> Self {
        Self::CS(Rc::new(RefCell::new(inner)))
    }

    // /// Returns the instance assignment of the constraint system
    // /// TODO:Fix the panic
    // pub fn instance_assignment(&self) -> Ref<[F]> {
    //     match self {
    //         ConstraintSystemRef::None => panic!(),
    //         ConstraintSystemRef::CS(cs) => {
    //             // Borrow the RefCell immutably and map to instance_assignment
    // slice             Ref::map(cs.borrow(), |cs_inner|
    // cs_inner.instance_assignment())         },
    //     }
    // }

    /// Returns the maximum arity of the local predicates.
    /// Maximum arity is used when stacking the local predicates as Garuda does
    pub fn max_arity(&self) -> usize {
        self.inner().map_or(0, |cs| cs.borrow().max_arity())
    }

    /// Returns the number of constraints which is the sum of the number of
    /// constraints in each local predicate.    #[inline]
    pub fn num_constraints(&self) -> usize {
        self.inner().map_or(0, |cs| cs.borrow().num_constraints())
    }

    /// Returns the number of instance variables.
    #[inline]
    pub fn num_instance_variables(&self) -> usize {
        self.inner()
            .map_or(0, |cs| cs.borrow().num_instance_variables())
    }

    /// Returns the number of witness variables.
    #[inline]
    pub fn num_witness_variables(&self) -> usize {
        self.inner()
            .map_or(0, |cs| cs.borrow().num_witness_variables())
    }

    /// Enforce a constraint in the constraint system. It takes a local
    /// predicate name and enforces a vector of linear combinations of the
    /// length of the arity of the local predicate enforces the constraint.
    #[inline]
    pub fn enforce_constraint(
        &self,
        local_predicate_label: &str,
        lc_vec: impl IntoIterator<Item = LinearCombination<F>>,
    ) -> crate::gr1cs::Result<()> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| {
                cs.borrow_mut()
                    .enforce_constraint(local_predicate_label, lc_vec)
            })
    }

    /// Enforce an r1cs constraint in the constraint system. It takes a, b, and
    /// c and enforces `a * b = c`. If R1CS predicate does not exist in the
    /// constraint system, It will create one. This function is a special case
    /// of `enforce_constraint` and is used as the legacy R1CS API to be bacward
    /// compatible with R1CS gadgets.
    /// TODO: Delete the inner enforce_r1cs_constraint function
    #[inline]
    pub fn enforce_r1cs_constraint(
        &self,
        a: LinearCombination<F>,
        b: LinearCombination<F>,
        c: LinearCombination<F>,
    ) -> crate::gr1cs::Result<()> {
        if !self.has_predicate(R1CS_PREDICATE_LABEL) {
            let r1cs_constraint_system =
                PredicateConstraintSystem::new_r1cs_predicate(self.clone())?;
            self.register_predicate(R1CS_PREDICATE_LABEL, r1cs_constraint_system)?;
        }
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| {
                cs.borrow_mut()
                    .enforce_constraint(R1CS_PREDICATE_LABEL, [a, b, c])
            })
    }

    /// Obtain a variable representing a linear combination.
    #[inline]
    pub fn new_lc(&self, lc: LinearCombination<F>) -> crate::gr1cs::Result<Variable> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| cs.borrow_mut().new_lc(lc))
    }
    /// Set `self.mode` to `mode`.
    /// Sets the mode if there exists an underlying ConstrainSystem
    pub fn set_mode(&self, mode: SynthesisMode) {
        self.inner().map_or((), |cs| cs.borrow_mut().set_mode(mode))
    }

    /// Check whether `self.mode == SynthesisMode::Setup`.
    /// Returns true if 1- There is an underlying ConstraintSystem and
    /// 2- It's in setup mode
    #[inline]
    pub fn is_in_setup_mode(&self) -> bool {
        self.inner()
            .is_some_and(|cs| cs.borrow().is_in_setup_mode())
    }

    /// Check whether this constraint system aims to optimize weight,
    /// number of constraints, or neither.
    #[inline]
    pub fn optimization_goal(&self) -> OptimizationGoal {
        self.inner().map_or(OptimizationGoal::Constraints, |cs| {
            cs.borrow().optimization_goal()
        })
    }

    /// Specify whether this constraint system should aim to optimize weight,
    /// number of constraints, or neither.
    #[inline]
    pub fn set_optimization_goal(&self, goal: OptimizationGoal) {
        self.inner()
            .map_or((), |cs| cs.borrow_mut().set_optimization_goal(goal))
    }

    /// Check whether or not `self` will construct matrices.
    #[inline]
    pub fn should_construct_matrices(&self) -> bool {
        self.inner()
            .is_some_and(|cs| cs.borrow().should_construct_matrices())
    }

    /// Obtain a variable representing a new public instance input
    /// This function takes a closure, this closure returns `Result<F>`
    /// Internally, this function calls new_input_variable of the constraint
    /// system to which it's pointing
    #[inline]
    pub fn new_input_variable<Func>(&self, f: Func) -> crate::utils::Result<Variable>
    where
        Func: FnOnce() -> crate::utils::Result<F>,
    {
        // mem_dbg("Before new_input_variable");
        let a = self
            .inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| {
                if !self.is_in_setup_mode() {
                    // This is needed to avoid double-borrows, because `f`
                    // might itself mutably borrow `cs` (eg: `f = || g.value()`).
                    let value = f();
                    cs.borrow_mut().new_input_variable(|| value)
                } else {
                    cs.borrow_mut().new_input_variable(f)
                }
            });
        // mem_dbg("After new_input_variable");
        a
    }

    /// Obtain a variable representing a new private witness input.
    #[inline]
    pub fn new_witness_variable<Func>(&self, f: Func) -> crate::utils::Result<Variable>
    where
        Func: FnOnce() -> crate::utils::Result<F>,
    {
        // mem_dbg("Before new_witness_variable");
        let a = self
            .inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| {
                if !self.is_in_setup_mode() {
                    // This is needed to avoid double-borrows, because `f`
                    // might itself mutably borrow `cs` (eg: `f = || g.value()`).
                    let value = f();
                    cs.borrow_mut().new_witness_variable(|| value)
                } else {
                    cs.borrow_mut().new_witness_variable(f)
                }
            });
        // mem_dbg("After new_witness_variable");
        a
    }

    /// Register a local predicate in the constraint system with a given label.
    pub fn register_predicate(
        &self,
        predicate_label: &str,
        predicate: PredicateConstraintSystem<F>,
    ) -> crate::utils::Result<()> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| {
                cs.borrow_mut()
                    .register_predicate(predicate_label, predicate)
            })
    }

    /// Checks if there is a predicate with the given label in the constraint
    pub fn has_predicate(&self, predicate_label: &str) -> bool {
        self.inner()
            .is_some_and(|cs| cs.borrow().has_predicate(predicate_label))
    }

    /// Obtain the assignment corresponding to the `Variable` `v`.
    pub fn assigned_value(&self, v: Variable) -> Option<F> {
        self.inner().and_then(|cs| cs.borrow().assigned_value(v))
    }

    /// If `self` is satisfied, outputs `Ok(true)`.
    /// If `self` is unsatisfied, outputs `Ok(false)`.
    /// If `self.is_in_setup_mode()` or if `self == None`, outputs `Err(())`.
    pub fn is_satisfied(&self) -> crate::utils::Result<bool> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| cs.borrow().is_satisfied())
    }

    /// If `self` is satisfied, outputs `Ok(None)`.
    /// If `self` is unsatisfied, outputs `Some(s,i)`, where `s` is the label of
    /// the unsatisfied local prediacate and  `i` is the index of
    /// the first unsatisfied constraint in that local predicate.
    /// If `self.is_in_setup_mode()` or `self == None`, outputs `Err(())`.
    pub fn which_predicate_is_unsatisfied(&self) -> crate::utils::Result<Option<String>> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| cs.borrow().which_predicate_is_unsatisfied())
    }

    /// Finalize the constraint system (either by outlining or inlining,
    /// if an optimization goal is set).
    pub fn finalize(&self) {
        if let Some(cs) = self.inner() {
            cs.borrow_mut().finalize()
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
    pub fn inline_all_lcs(&self) {
        if let Some(cs) = self.inner() {
            cs.borrow_mut().inline_all_lcs()
        }
    }

    /// Returns `self` if `!self.is_none()`, otherwise returns `other`.
    pub fn or(self, other: Self) -> Self {
        match self {
            ConstraintSystemRef::None => other,
            _ => self,
        }
    }

    /// Returns `true` is `self == ConstraintSystemRef::None`.
    pub fn is_none(&self) -> bool {
        matches!(self, ConstraintSystemRef::None)
    }

    pub(crate) fn inner(&self) -> Option<&Rc<RefCell<ConstraintSystem<F>>>> {
        match self {
            Self::CS(a) => Some(a),
            Self::None => None,
        }
    }

    /// Consumes self to return the inner `ConstraintSystem<F>`. Returns
    /// `None` if `Self::CS` is `None` or if any other references to
    /// `Self::CS` exist.  
    pub fn into_inner(self) -> Option<ConstraintSystem<F>> {
        match self {
            Self::CS(a) => Rc::try_unwrap(a).ok().map(|s| s.into_inner()),
            Self::None => None,
        }
    }

    /// Get the matrices corresponding to the local predicates.and the
    /// corresponding set of matrices
    #[inline]
    pub fn to_matrices(&self) -> crate::gr1cs::Result<BTreeMap<Label, Vec<Matrix<F>>>> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| cs.borrow().to_matrices())
    }

    /// Get the linear combination corresponding to the given `lc_index`.
    /// TODO: This function should return a reference to the linear combination
    /// and not clone it.
    pub fn get_lc(&self, lc_index: LcIndex) -> crate::utils::Result<LinearCombination<F>> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| cs.borrow().get_lc(lc_index))
    }


    pub fn make_row(&self, lc: &LinearCombination<F>) -> crate::utils::Result<Vec<(F, usize)>> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| Ok(cs.borrow().make_row(lc)))
    }

    // TODO: Implement this function
    // /// Get trace information about all constraints in the system
    // pub fn constraint_names(&self) -> Option<Vec<String>> {
    //     #[cfg(feature = "std")]
    //     {
    //         self.inner().and_then(|cs| {
    //             cs.borrow()
    //                 .constraint_traces
    //                 .iter()
    //                 .map(|trace| {
    //                     let mut constraint_path = String::new();
    //                     let mut prev_module_path = "";
    //                     let mut prefixes = ark_std::collections::BTreeSet::new();
    //                     for step in trace.as_ref()?.path() {
    //                         let module_path = if prev_module_path ==
    // step.module_path {
    // prefixes.insert(step.module_path.to_string());
    // String::new()                         } else {
    //                             let mut parts = step
    //                                 .module_path
    //                                 .split("::")
    //                                 .filter(|&part| part != "r1cs_std" && part !=
    // "constraints");                             let mut path_so_far =
    // String::new();                             for part in parts.by_ref() {
    //                                 if path_so_far.is_empty() {
    //                                     path_so_far += part;
    //                                 } else {
    //                                     path_so_far += &["::", part].join("");
    //                                 }
    //                                 if prefixes.contains(&path_so_far) {
    //                                     continue;
    //                                 } else {
    //                                     prefixes.insert(path_so_far.clone());
    //                                     break;
    //                                 }
    //                             }
    //                             parts.collect::<Vec<_>>().join("::") + "::"
    //                         };
    //                         prev_module_path = step.module_path;
    //                         constraint_path += &["/", &module_path,
    // step.name].join("");                     }
    //                     Some(constraint_path)
    //                 })
    //                 .collect::<Option<Vec<_>>>()
    //         })
    //     }
    //     #[cfg(not(feature = "std"))]
    //     {
    //         None
    //     }
    // }
}

/// A shared reference to a constraint system that can be stored in high level
/// variables.
#[derive(Debug, Clone)]
pub enum WeakConstraintSystemRef<F: Field> {
    /// Represents the case where we *don't* need to allocate variables or
    /// enforce constraints. Encountered when operating over constant
    /// values.
    None,
    /// Represents the case where we *do* allocate variables or enforce
    /// constraints.
    CS(Weak<RefCell<ConstraintSystem<F>>>),
}

impl<F: Field> WeakConstraintSystemRef<F> {
    /// Create a new `WeakConstraintSystemRef` from a
    /// `Weak<RefCell<ConstraintSystem<F>>>`
    pub fn from(cs: ConstraintSystemRef<F>) -> Self {
        match cs.inner() {
            None => WeakConstraintSystemRef::None,
            Some(cs) => WeakConstraintSystemRef::CS(Rc::downgrade(cs)),
        }
    }

    /// Convert a `WeakConstraintSystemRef` to a `ConstraintSystemRef`
    /// If the `WeakConstraintSystemRef` is `None`, it will return `None`
    /// Otherwise, it will upgrade the weak reference
    pub fn to_constraint_system_ref(&self) -> ConstraintSystemRef<F> {
        match self {
            WeakConstraintSystemRef::CS(ref cs) => ConstraintSystemRef::CS(cs.upgrade().unwrap()),
            WeakConstraintSystemRef::None => ConstraintSystemRef::None,
        }
    }
}
