#[cfg(feature = "std")]
use crate::ConstraintTrace;
use crate::{
    format, vec, BTreeMap, LcIndex, LinearCombination, Matrix, Rc, String, SynthesisError,
    Variable, Vec,
};
use algebra_core::Field;
use core::cell::{Ref, RefCell, RefMut};

/// Computations are expressed in terms of rank-1 constraint systems (R1CS).
/// The `generate_constraints` method is called to generate constraints for
/// both CRS generation and for proving.
// TODO: Think: should we replace this with just a closure?
pub trait ConstraintSynthesizer<F: Field> {
    /// Drives generation of new constraints inside `cs`.
    fn generate_constraints(self, cs: ConstraintSystemRef<F>) -> Result<(), SynthesisError>;
}

/// An Rank-One `ConstraintSystem`. Enforces constraints of the form
/// `⟨a_i, z⟩ ⋅ ⟨b_i, z⟩ = ⟨c_i, z⟩`, where `a_i`, `b_i`, and `c_i` are linear
/// combinations over variables, and `z` is the concrete assignment to these
/// variables.
#[derive(Debug, Clone)]
pub struct ConstraintSystem<F: Field> {
    /// The mode in which the constraint system is operating. `self` can either
    /// be in setup mode (i.e., `self.mode == SynthesisMode::Setup`) or in
    /// proving mode (i.e., `self.mode == SynthesisMode::Prove`). If we are
    /// in proving mode, then we have the additional option of whether or
    /// not to construct the A, B, and C matrices of the constraint system
    /// (see below).
    pub mode: SynthesisMode,
    /// The number of variables that are "public inputs" to the constraint
    /// system.
    pub num_instance_variables: usize,
    /// The number of variables that are "private inputs" to the constraint
    /// system.
    pub num_witness_variables: usize,
    /// The number of constraints in the constraint system.
    pub num_constraints: usize,
    /// The number of linear combinations
    pub num_linear_combinations: usize,

    /// Assignments to the public input variables. This is empty if `self.mode
    /// == SynthesisMode::Setup`.
    pub instance_assignment: Vec<F>,
    /// Assignments to the private input variables. This is empty if `self.mode
    /// == SynthesisMode::Setup`.
    pub witness_assignment: Vec<F>,

    lc_map: BTreeMap<LcIndex, LinearCombination<F>>,

    #[cfg(feature = "std")]
    constraint_traces: Vec<Option<ConstraintTrace>>,

    a_constraints: Vec<LcIndex>,
    b_constraints: Vec<LcIndex>,
    c_constraints: Vec<LcIndex>,

    lc_assignment_cache: Rc<RefCell<BTreeMap<LcIndex, F>>>,
}

/// Defines the mode of operation of a `ConstraintSystem`.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum SynthesisMode {
    /// Indicate to the `ConstraintSystem` that it should only generate
    /// constraint matrices and not populate the variable assignments.
    Setup,
    /// Indicate to the `ConstraintSystem` that it populate the variable
    /// assignments. If additionally `construct_matrices == true`, then generate
    /// the matrices as in the `Setup` case.
    Prove {
        /// If `construct_matrices == true`, then generate
        /// the matrices as in the `Setup` case.
        construct_matrices: bool,
    },
}

impl<F: Field> ConstraintSystem<F> {
    #[inline]
    fn make_row(&self, l: &LinearCombination<F>) -> Vec<(F, usize)> {
        let num_input = self.num_instance_variables;
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

    /// Construct an ampty `ConstraintSystem`.
    pub fn new() -> Self {
        Self {
            num_instance_variables: 1,
            num_witness_variables: 0,
            num_constraints: 0,
            num_linear_combinations: 0,
            a_constraints: Vec::new(),
            b_constraints: Vec::new(),
            c_constraints: Vec::new(),
            instance_assignment: vec![F::one()],
            witness_assignment: Vec::new(),
            #[cfg(feature = "std")]
            constraint_traces: Vec::new(),

            lc_map: BTreeMap::new(),
            lc_assignment_cache: Rc::new(RefCell::new(BTreeMap::new())),

            mode: SynthesisMode::Prove {
                construct_matrices: true,
            },
        }
    }

    /// Create a new `ConstraintSystemRef<F>`.
    pub fn new_ref() -> ConstraintSystemRef<F> {
        ConstraintSystemRef::new(Self::new())
    }

    /// Set `self.mode` to `mode`.
    pub fn set_mode(&mut self, mode: SynthesisMode) {
        self.mode = mode;
    }

    /// Check whether `self.mode == SynthesisMode::Setup`.
    pub fn is_in_setup_mode(&self) -> bool {
        self.mode == SynthesisMode::Setup
    }

    /// Check whether or not `self` will construct matrices.
    pub fn should_construct_matrices(&self) -> bool {
        match self.mode {
            SynthesisMode::Setup => true,
            SynthesisMode::Prove { construct_matrices } => construct_matrices,
        }
    }

    /// Return a variable representing the constant "zero" inside the constraint
    /// system.
    #[inline]
    pub fn zero() -> Variable {
        Variable::Zero
    }

    /// Return a variable representing the constant "one" inside the constraint
    /// system.
    #[inline]
    pub fn one() -> Variable {
        Variable::One
    }

    /// Obtain a variable representing a new public instance input.
    #[inline]
    pub fn new_input_variable<Func>(&mut self, f: Func) -> Result<Variable, SynthesisError>
    where
        Func: FnOnce() -> Result<F, SynthesisError>,
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
    pub fn new_witness_variable<Func>(&mut self, f: Func) -> Result<Variable, SynthesisError>
    where
        Func: FnOnce() -> Result<F, SynthesisError>,
    {
        let index = self.num_witness_variables;
        self.num_witness_variables += 1;

        if !self.is_in_setup_mode() {
            self.witness_assignment.push(f()?);
        }
        Ok(Variable::Witness(index))
    }

    /// Obtain a variable representing a linear combination.
    #[inline]
    pub fn new_lc(&mut self, lc: LinearCombination<F>) -> Result<Variable, SynthesisError> {
        let index = LcIndex(self.num_linear_combinations);
        let var = Variable::SymbolicLc(index);

        self.lc_map.insert(index, lc);

        self.num_linear_combinations += 1;
        Ok(var)
    }

    /// Enforce a R1CS constraint with the name `name`.
    #[inline]
    pub fn enforce_constraint(
        &mut self,
        a: LinearCombination<F>,
        b: LinearCombination<F>,
        c: LinearCombination<F>,
    ) -> Result<(), SynthesisError> {
        if self.should_construct_matrices() {
            let a_index = self.new_lc(a)?.get_lc_index().unwrap();
            let b_index = self.new_lc(b)?.get_lc_index().unwrap();
            let c_index = self.new_lc(c)?.get_lc_index().unwrap();
            self.a_constraints.push(a_index);
            self.b_constraints.push(b_index);
            self.c_constraints.push(c_index);
        }
        self.num_constraints += 1;
        #[cfg(feature = "std")]
        {
            let trace = ConstraintTrace::capture();
            self.constraint_traces.push(trace);
        }
        Ok(())
    }

    /// Naively inlines symbolic linear combinations into the linear
    /// combinations that use them.
    ///
    /// Useful for standard pairing-based SNARKs where addition gates are cheap.
    /// For example, in the SNARKs such as [[Groth16]](https://eprint.iacr.org/2016/260) and
    /// [[Groth-Maller17]](https://eprint.iacr.org/2017/540), addition gates
    /// do not contribute to the size of the multi-scalar multiplication, which
    /// is the dominating cost.
    pub fn inline_all_lcs(&mut self) {
        let mut inlined_lcs = BTreeMap::new();
        for (&index, lc) in &self.lc_map {
            let mut inlined_lc = LinearCombination::new();
            for &(coeff, var) in lc.iter() {
                if var.is_lc() {
                    let lc_index = var.get_lc_index().expect("should be lc");
                    // If `var` is a `SymbolicLc`, fetch the corresponding
                    // inlined LC, and substitute it in.
                    let lc = inlined_lcs.get(&lc_index).expect("should be inlined");
                    inlined_lc.extend((lc * coeff).0.into_iter());
                } else {
                    // Otherwise, it's a concrete variable and so we
                    // substitute it in directly.
                    inlined_lc.push((coeff, var));
                }
            }
            inlined_lc.compactify();
            inlined_lcs.insert(index, inlined_lc);
        }
        self.lc_map = inlined_lcs;
    }

    /// If a `SymbolicLc` is used in more than one location, this method makes a
    /// new variable for that `SymbolicLc`, adds a constraint ensuring the
    /// equality of the variable and the linear combination, and then uses
    /// that variable in every location the `SymbolicLc` is used.
    ///
    /// Useful for SNARKs like `Marlin` or `Fractal`, where addition gates
    /// are not cheap.
    pub fn outline_lcs(&mut self) {
        unimplemented!()
    }

    /// This step must be called after constraint generation has completed, and
    /// after all symbolic LCs have been inlined into the places that they
    /// are used.
    pub fn to_matrices(&self) -> Option<ConstraintMatrices<F>> {
        if let SynthesisMode::Prove {
            construct_matrices: false,
        } = self.mode
        {
            None
        } else {
            let a: Vec<_> = self
                .a_constraints
                .iter()
                .map(|index| self.make_row(self.lc_map.get(index).unwrap()))
                .collect();
            let b: Vec<_> = self
                .b_constraints
                .iter()
                .map(|index| self.make_row(self.lc_map.get(index).unwrap()))
                .collect();
            let c: Vec<_> = self
                .c_constraints
                .iter()
                .map(|index| self.make_row(self.lc_map.get(index).unwrap()))
                .collect();

            let a_num_non_zero: usize = a.iter().map(|lc| lc.len()).sum();
            let b_num_non_zero: usize = b.iter().map(|lc| lc.len()).sum();
            let c_num_non_zero: usize = c.iter().map(|lc| lc.len()).sum();
            let matrices = ConstraintMatrices {
                num_instance_variables: self.num_instance_variables,
                num_witness_variables: self.num_witness_variables,
                num_constraints: self.num_constraints,

                a_num_non_zero,
                b_num_non_zero,
                c_num_non_zero,

                a,
                b,
                c,
            };
            Some(matrices)
        }
    }

    fn eval_lc(&self, lc: LcIndex) -> Option<F> {
        let lc = self.lc_map.get(&lc)?;
        let mut acc = F::zero();
        for (coeff, var) in lc.iter() {
            acc += *coeff * &self.assigned_value(*var)?;
        }
        Some(acc)
    }

    /// If `self` is satisfied, outputs `Ok(true)`.
    /// If `self` is unsatisfied, outputs `Ok(false)`.
    /// If `self.is_in_setup_mode()`, outputs `Err(())`.
    pub fn is_satisfied(&self) -> Result<bool, SynthesisError> {
        self.which_is_unsatisfied().map(|s| s.is_none())
    }

    /// If `self` is satisfied, outputs `Ok(None)`.
    /// If `self` is unsatisfied, outputs `Some(i)`, where `i` is the index of
    /// the first unsatisfied constraint. If `self.is_in_setup_mode()`, outputs
    /// `Err(())`.
    pub fn which_is_unsatisfied(&self) -> Result<Option<String>, SynthesisError> {
        if self.is_in_setup_mode() {
            Err(SynthesisError::AssignmentMissing)
        } else {
            for i in 0..self.num_constraints {
                let a = self
                    .eval_lc(self.a_constraints[i])
                    .ok_or(SynthesisError::AssignmentMissing)?;
                let b = self
                    .eval_lc(self.b_constraints[i])
                    .ok_or(SynthesisError::AssignmentMissing)?;
                let c = self
                    .eval_lc(self.c_constraints[i])
                    .ok_or(SynthesisError::AssignmentMissing)?;
                if a * b != c {
                    let trace;
                    #[cfg(feature = "std")]
                    {
                        trace = self.constraint_traces[i].as_ref().map_or_else(
                            || {
                                eprintln!("Constraint trace requires enabling `ConstraintLayer`");
                                format!("{}", i)
                            },
                            |t| format!("{}", t),
                        );
                    }
                    #[cfg(not(feature = "std"))]
                    {
                        trace = format!("{}", i);
                    }
                    return Ok(Some(trace));
                }
            }
            Ok(None)
        }
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
            }
        }
    }
}
/// The A, B and C matrices of a Rank-One `ConstraintSystem`.
/// Also contains metadata on the structure of the constraint system
/// and the matrices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstraintMatrices<F: Field> {
    /// The number of variables that are "public instances" to the constraint
    /// system.
    pub num_instance_variables: usize,
    /// The number of variables that are "private witnesses" to the constraint
    /// system.
    pub num_witness_variables: usize,
    /// The number of constraints in the constraint system.
    pub num_constraints: usize,
    /// The number of non_zero entries in the A matrix.
    pub a_num_non_zero: usize,
    /// The number of non_zero entries in the B matrix.
    pub b_num_non_zero: usize,
    /// The number of non_zero entries in the C matrix.
    pub c_num_non_zero: usize,

    /// The A constraint matrix. This is empty when
    /// `self.mode == SynthesisMode::Prove { construct_matrices = false }`.
    pub a: Matrix<F>,
    /// The B constraint matrix. This is empty when
    /// `self.mode == SynthesisMode::Prove { construct_matrices = false }`.
    pub b: Matrix<F>,
    /// The C constraint matrix. This is empty when
    /// `self.mode == SynthesisMode::Prove { construct_matrices = false }`.
    pub c: Matrix<F>,
}

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

/// A namespaced `ConstraintSystemRef`.
#[derive(Debug, Clone)]
pub struct Namespace<F: Field> {
    inner: ConstraintSystemRef<F>,
    #[cfg(feature = "std")]
    id: Option<tracing::Id>,
}

impl<F: Field> From<ConstraintSystemRef<F>> for Namespace<F> {
    fn from(other: ConstraintSystemRef<F>) -> Self {
        Self {
            inner: other,
            #[cfg(feature = "std")]
            id: None,
        }
    }
}

impl<F: Field> Namespace<F> {
    /// Construct a new `Namespace`.
    #[cfg(feature = "std")]
    pub fn new(inner: ConstraintSystemRef<F>, id: Option<tracing::Id>) -> Self {
        Self { inner, id }
    }

    /// Obtain the inner `ConstraintSystemRef<F>`.
    pub fn cs(&self) -> ConstraintSystemRef<F> {
        self.inner.clone()
    }

    /// Manually leave the namespace.
    pub fn leave_namespace(mut self) {
        drop(&mut self)
    }
}

#[cfg(feature = "std")]
impl<F: Field> Drop for Namespace<F> {
    fn drop(&mut self) {
        if let Some(id) = self.id.as_ref() {
            tracing::dispatcher::get_default(|dispatch| dispatch.exit(id))
        }
        drop(&mut self.inner)
    }
}

impl<F: Field> ConstraintSystemRef<F> {
    /// Returns `self` if `!self.is_none()`, otherwise returns `other`.
    pub fn or(self, other: Self) -> Self {
        match self {
            ConstraintSystemRef::None => other,
            _ => self,
        }
    }

    /// Returns `true` is `self == ConstraintSystemRef::None`.
    pub fn is_none(&self) -> bool {
        match self {
            ConstraintSystemRef::None => true,
            _ => false,
        }
    }

    /// Construct a `ConstraintSystemRef` from a `ConstraintSystem`.
    #[inline]
    pub fn new(inner: ConstraintSystem<F>) -> Self {
        Self::CS(Rc::new(RefCell::new(inner)))
    }

    fn inner(&self) -> Option<&Rc<RefCell<ConstraintSystem<F>>>> {
        match self {
            Self::CS(a) => Some(a),
            Self::None => None,
        }
    }

    /// Obtain an immutable reference to the underlying `ConstraintSystem`.
    ///
    /// # Panics
    /// This method panics if `self` is already mutably borrowed.
    #[inline]
    pub fn borrow(&self) -> Option<Ref<ConstraintSystem<F>>> {
        self.inner().map(|cs| cs.borrow())
    }

    /// Obtain a mutable reference to the underlying `ConstraintSystem`.
    ///
    /// # Panics
    /// This method panics if `self` is already mutably borrowed.
    #[inline]
    pub fn borrow_mut(&self) -> Option<RefMut<ConstraintSystem<F>>> {
        self.inner().map(|cs| cs.borrow_mut())
    }

    /// Set `self.mode` to `mode`.
    pub fn set_mode(&self, mode: SynthesisMode) {
        self.inner().map_or((), |cs| cs.borrow_mut().set_mode(mode))
    }

    /// Check whether `self.mode == SynthesisMode::Setup`.
    #[inline]
    pub fn is_in_setup_mode(&self) -> bool {
        self.inner()
            .map_or(false, |cs| cs.borrow().is_in_setup_mode())
    }

    /// Returns the number of constraints.
    #[inline]
    pub fn num_constraints(&self) -> usize {
        self.inner().map_or(0, |cs| cs.borrow().num_constraints)
    }

    /// Returns the number of instance variables.
    #[inline]
    pub fn num_instance_variables(&self) -> usize {
        self.inner()
            .map_or(0, |cs| cs.borrow().num_instance_variables)
    }

    /// Returns the number of witness variables.
    #[inline]
    pub fn num_witness_variables(&self) -> usize {
        self.inner()
            .map_or(0, |cs| cs.borrow().num_witness_variables)
    }

    /// Check whether or not `self` will construct matrices.
    #[inline]
    pub fn should_construct_matrices(&self) -> bool {
        self.inner()
            .map_or(false, |cs| cs.borrow().should_construct_matrices())
    }

    /// Obtain a variable representing a new public instance input.
    #[inline]
    pub fn new_input_variable<Func>(&self, f: Func) -> Result<Variable, SynthesisError>
    where
        Func: FnOnce() -> Result<F, SynthesisError>,
    {
        self.inner()
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
            })
    }

    /// Obtain a variable representing a new private witness input.
    #[inline]
    pub fn new_witness_variable<Func>(&self, f: Func) -> Result<Variable, SynthesisError>
    where
        Func: FnOnce() -> Result<F, SynthesisError>,
    {
        self.inner()
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
            })
    }

    /// Obtain a variable representing a linear combination.
    #[inline]
    pub fn new_lc(&self, lc: LinearCombination<F>) -> Result<Variable, SynthesisError> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| cs.borrow_mut().new_lc(lc))
    }

    /// Enforce a R1CS constraint with the name `name`.
    #[inline]
    pub fn enforce_constraint(
        &self,
        a: LinearCombination<F>,
        b: LinearCombination<F>,
        c: LinearCombination<F>,
    ) -> Result<(), SynthesisError> {
        self.inner()
            .ok_or(SynthesisError::MissingCS)
            .and_then(|cs| cs.borrow_mut().enforce_constraint(a, b, c))
    }

    /// Naively inlines symbolic linear combinations into the linear
    /// combinations that use them.
    ///
    /// Useful for standard pairing-based SNARKs where addition gates are free,
    /// such as the SNARKs in [[Groth16]](https://eprint.iacr.org/2016/260) and
    /// [[Groth-Maller17]](https://eprint.iacr.org/2017/540).
    pub fn inline_all_lcs(&self) {
        if let Some(cs) = self.inner() {
            cs.borrow_mut().inline_all_lcs()
        }
    }

    /// If a `SymbolicLc` is used in more than one location, this method makes a
    /// new variable for that `SymbolicLc`, adds a constraint ensuring the
    /// equality of the variable and the linear combination, and then uses
    /// that variable in every location the `SymbolicLc` is used.
    ///
    /// Useful for SNARKs like `Marlin` or `Fractal`, where where addition gates
    /// are not (entirely) free.
    pub fn outline_lcs(&self) {
        if let Some(cs) = self.inner() {
            cs.borrow_mut().outline_lcs()
        }
    }

    /// This step must be called after constraint generation has completed, and
    /// after all symbolic LCs have been inlined into the places that they
    /// are used.
    #[inline]
    pub fn to_matrices(&self) -> Option<ConstraintMatrices<F>> {
        self.inner().map_or(None, |cs| cs.borrow().to_matrices())
    }

    /// If `self` is satisfied, outputs `Ok(true)`.
    /// If `self` is unsatisfied, outputs `Ok(false)`.
    /// If `self.is_in_setup_mode()` or if `self == None`, outputs `Err(())`.
    pub fn is_satisfied(&self) -> Result<bool, SynthesisError> {
        self.inner()
            .map_or(Err(SynthesisError::AssignmentMissing), |cs| {
                cs.borrow().is_satisfied()
            })
    }

    /// If `self` is satisfied, outputs `Ok(None)`.
    /// If `self` is unsatisfied, outputs `Some(i)`, where `i` is the index of
    /// the first unsatisfied constraint.
    /// If `self.is_in_setup_mode()` or `self == None`, outputs `Err(())`.
    pub fn which_is_unsatisfied(&self) -> Result<Option<String>, SynthesisError> {
        self.inner()
            .map_or(Err(SynthesisError::AssignmentMissing), |cs| {
                cs.borrow().which_is_unsatisfied()
            })
    }

    /// Obtain the assignment corresponding to the `Variable` `v`.
    pub fn assigned_value(&self, v: Variable) -> Option<F> {
        self.inner()
            .map_or(None, |cs| cs.borrow().assigned_value(v))
    }

    /// Get trace information about all constraints in the system
    pub fn constraint_names(&self) -> Option<Vec<String>> {
        #[cfg(feature = "std")]
        {
            self.inner().and_then(|cs| {
                cs.borrow()
                    .constraint_traces
                    .iter()
                    .map(|trace| {
                        let mut constraint_path = String::new();
                        let mut prev_module_path = "";
                        let mut prefixes = crate::BTreeSet::new();
                        for step in trace.as_ref()?.path() {
                            let module_path = if prev_module_path == step.module_path {
                                prefixes.insert(step.module_path.to_string());
                                String::new()
                            } else {
                                let mut parts = step
                                    .module_path
                                    .split("::")
                                    .filter(|&part| part != "r1cs_std" && part != "constraints");
                                let mut path_so_far = String::new();
                                for part in parts.by_ref() {
                                    if path_so_far.is_empty() {
                                        path_so_far += part;
                                    } else {
                                        path_so_far += &["::", part].join("");
                                    }
                                    if prefixes.contains(&path_so_far) {
                                        continue;
                                    } else {
                                        prefixes.insert(path_so_far.clone());
                                        break;
                                    }
                                }
                                parts.collect::<Vec<_>>().join("::") + "::"
                            };
                            prev_module_path = step.module_path;
                            constraint_path += &["/", &module_path, step.name].join("");
                        }
                        Some(constraint_path)
                    })
                    .collect::<Option<Vec<_>>>()
            })
        }
        #[cfg(not(feature = "std"))]
        {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    use algebra::{bls12_381::Fr, One};

    #[test]
    fn matrix_generation() -> Result<(), SynthesisError> {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let two = Fr::one() + Fr::one();
        let a = cs.new_input_variable(|| Ok(Fr::one()))?;
        let b = cs.new_witness_variable(|| Ok(Fr::one()))?;
        let c = cs.new_witness_variable(|| Ok(two))?;
        cs.enforce_constraint(lc!() + a, lc!() + (two, b), lc!() + c)?;
        let d = cs.new_lc(lc!() + a + b)?;
        cs.enforce_constraint(lc!() + a, lc!() + d, lc!() + d)?;
        let e = cs.new_lc(lc!() + d + d)?;
        cs.enforce_constraint(lc!() + Variable::One, lc!() + e, lc!() + e)?;
        cs.inline_all_lcs();
        let matrices = cs.to_matrices().unwrap();
        assert_eq!(matrices.a[0], vec![(Fr::one(), 1)]);
        assert_eq!(matrices.b[0], vec![(two, 2)]);
        assert_eq!(matrices.c[0], vec![(Fr::one(), 3)]);

        assert_eq!(matrices.a[1], vec![(Fr::one(), 1)]);
        assert_eq!(matrices.b[1], vec![(Fr::one(), 1), (Fr::one(), 2)]);
        assert_eq!(matrices.c[1], vec![(Fr::one(), 1), (Fr::one(), 2)]);

        assert_eq!(matrices.a[2], vec![(Fr::one(), 0)]);
        assert_eq!(matrices.b[2], vec![(two, 1), (two, 2)]);
        assert_eq!(matrices.c[2], vec![(two, 1), (two, 2)]);
        Ok(())
    }
}
