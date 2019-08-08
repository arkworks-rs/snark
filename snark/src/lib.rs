#![deny(
    unused_import_braces,
    unused_qualifications,
    trivial_casts,
    trivial_numeric_casts
)]
#![deny(
    unused_qualifications,
    variant_size_differences,
    stable_features,
    unreachable_pub
)]
#![deny(
    non_shorthand_field_patterns,
    unused_attributes,
    unused_imports,
    unused_extern_crates
)]
#![deny(
    renamed_and_removed_lints,
    stable_features,
    unused_allocation,
    unused_comparisons,
    bare_trait_objects
)]
#![deny(const_err, unused_must_use, unused_mut, unused_unsafe, private_in_public)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate bench_utils;

pub mod gm17;

use algebra::{Field, PairingEngine};

use smallvec::{smallvec, SmallVec as StackVec};
use std::{
    cmp::Ordering,
    error::Error,
    fmt, io,
    marker::PhantomData,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub},
};

type SmallVec<E> = StackVec<[(Variable, <E as PairingEngine>::Fr); 16]>;

/// Computations are expressed in terms of arithmetic circuits, in particular
/// rank-1 quadratic constraint systems. The `Circuit` trait represents a
/// circuit that can be synthesized. The `synthesize` method is called during
/// CRS generation and during proving.
pub trait Circuit<E: PairingEngine> {
    /// Synthesize the circuit into a rank-1 quadratic constraint system
    fn synthesize<CS: ConstraintSystem<E>>(self, cs: &mut CS) -> Result<(), SynthesisError>;
}

/// Represents a variable in our constraint system.
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct Variable(Index);

impl PartialOrd for Variable {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Variable {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self.0, other.0) {
            (Index::Input(ref idx1), Index::Input(ref idx2))
            | (Index::Aux(ref idx1), Index::Aux(ref idx2)) => idx1.cmp(idx2),
            (Index::Input(_), Index::Aux(_)) => Ordering::Less,
            (Index::Aux(_), Index::Input(_)) => Ordering::Greater,
        }
    }
}

impl Variable {
    /// This constructs a variable with an arbitrary index.
    /// Circuit implementations are not recommended to use this.
    pub fn new_unchecked(idx: Index) -> Variable {
        Variable(idx)
    }

    /// This returns the index underlying the variable.
    /// Circuit implementations are not recommended to use this.
    pub fn get_unchecked(&self) -> Index {
        self.0
    }
}

/// Represents the index of either an input variable or
/// auxiliary variable.
#[derive(Copy, Clone, PartialEq, Debug, Eq)]
pub enum Index {
    Input(usize),
    Aux(usize),
}

/// This represents a linear combination of some variables, with coefficients
/// in the scalar field of a pairing-friendly elliptic curve group.
/// The `(coeff, var)` pairs in a `LinearCombination` are kept sorted according
/// to the index of the variable.
#[derive(Debug, Clone)]
pub struct LinearCombination<E: PairingEngine>(pub SmallVec<E>);

impl<E: PairingEngine> AsRef<[(Variable, E::Fr)]> for LinearCombination<E> {
    #[inline]
    fn as_ref(&self) -> &[(Variable, E::Fr)] {
        &self.0
    }
}

impl<E: PairingEngine> From<(E::Fr, Variable)> for LinearCombination<E> {
    #[inline]
    fn from((coeff, var): (E::Fr, Variable)) -> Self {
        LinearCombination(smallvec![(var, coeff)])
    }
}

impl<E: PairingEngine> From<Variable> for LinearCombination<E> {
    #[inline]
    fn from(var: Variable) -> Self {
        LinearCombination(smallvec![(var, E::Fr::one())])
    }
}

impl<E: PairingEngine> LinearCombination<E> {
    #[inline]
    pub fn zero() -> LinearCombination<E> {
        LinearCombination(SmallVec::<E>::new())
    }

    #[inline]
    pub fn replace_in_place(&mut self, other: Self) {
        self.0.truncate(0);
        self.0.extend_from_slice(&other.0)
    }

    #[inline]
    pub fn negate_in_place(&mut self) {
        self.0.iter_mut().for_each(|(_, coeff)| *coeff = -(*coeff));
    }

    #[inline]
    pub fn double_in_place(&mut self) {
        self.0.iter_mut().for_each(|(_, coeff)| {
            coeff.double_in_place();
        });
    }

    #[inline]
    pub fn get_var_loc(&self, search_var: &Variable) -> Result<usize, usize> {
        if self.0.len() < 6 {
            let mut found_index = 0;
            for (i, (var, _)) in self.0.iter().enumerate() {
                if var >= search_var {
                    found_index = i;
                    break;
                } else {
                    found_index += 1;
                }
            }
            return Err(found_index);
        } else {
            self.0
                .binary_search_by_key(search_var, |&(cur_var, _)| cur_var)
        }
    }
}

impl<E: PairingEngine> Add<(E::Fr, Variable)> for LinearCombination<E> {
    type Output = Self;

    #[inline]
    fn add(mut self, coeff_var: (E::Fr, Variable)) -> Self {
        self += coeff_var;
        self
    }
}

impl<E: PairingEngine> AddAssign<(E::Fr, Variable)> for LinearCombination<E> {
    #[inline]
    fn add_assign(&mut self, (coeff, var): (E::Fr, Variable)) {
        match self.get_var_loc(&var) {
            Ok(found) => self.0[found].1 += &coeff,
            Err(not_found) => self.0.insert(not_found, (var, coeff)),
        }
    }
}

impl<E: PairingEngine> Sub<(E::Fr, Variable)> for LinearCombination<E> {
    type Output = Self;

    #[inline]
    fn sub(self, (coeff, var): (E::Fr, Variable)) -> Self {
        self + (-coeff, var)
    }
}

impl<E: PairingEngine> Neg for LinearCombination<E> {
    type Output = Self;

    #[inline]
    fn neg(mut self) -> Self {
        self.negate_in_place();
        self
    }
}

impl<E: PairingEngine> Mul<E::Fr> for LinearCombination<E> {
    type Output = Self;

    #[inline]
    fn mul(mut self, scalar: E::Fr) -> Self {
        self *= scalar;
        self
    }
}

impl<E: PairingEngine> MulAssign<E::Fr> for LinearCombination<E> {
    #[inline]
    fn mul_assign(&mut self, scalar: E::Fr) {
        self.0.iter_mut().for_each(|(_, coeff)| *coeff *= &scalar);
    }
}

impl<E: PairingEngine> Add<Variable> for LinearCombination<E> {
    type Output = Self;

    #[inline]
    fn add(self, other: Variable) -> LinearCombination<E> {
        self + (E::Fr::one(), other)
    }
}

impl<E: PairingEngine> Sub<Variable> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    #[inline]
    fn sub(self, other: Variable) -> LinearCombination<E> {
        self - (E::Fr::one(), other)
    }
}

fn op_impl<E: PairingEngine, F1, F2>(
    cur: &LinearCombination<E>,
    other: &LinearCombination<E>,
    push_fn: F1,
    combine_fn: F2,
) -> LinearCombination<E>
where
    F1: Fn(E::Fr) -> E::Fr,
    F2: Fn(E::Fr, E::Fr) -> E::Fr,
{
    let mut new_vec = SmallVec::<E>::new(); // with_capacity($self.0.len() + $other.0.len());
    let mut i = 0;
    let mut j = 0;
    while i < cur.0.len() && j < other.0.len() {
        let self_cur = &cur.0[i];
        let other_cur = &other.0[j];
        if self_cur.0 > other_cur.0 {
            new_vec.push((other.0[j].0, push_fn(other.0[j].1)));
            j += 1;
        } else if self_cur.0 < other_cur.0 {
            new_vec.push(*self_cur);
            i += 1;
        } else {
            new_vec.push((self_cur.0, combine_fn(self_cur.1, other_cur.1)));
            i += 1;
            j += 1;
        }
    }
    new_vec.extend_from_slice(&cur.0[i..]);
    while j < other.0.len() {
        new_vec.push((other.0[j].0, push_fn(other.0[j].1)));
        j += 1;
    }
    LinearCombination(new_vec)
}

impl<E: PairingEngine> Add<&LinearCombination<E>> for &LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(self, other: &LinearCombination<E>) -> LinearCombination<E> {
        if other.0.is_empty() {
            return self.clone();
        } else if self.0.is_empty() {
            return other.clone();
        }
        op_impl(
            self,
            other,
            |coeff| coeff,
            |cur_coeff, other_coeff| cur_coeff + &other_coeff,
        )
    }
}

impl<E: PairingEngine> Add<LinearCombination<E>> for &LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(self, other: LinearCombination<E>) -> LinearCombination<E> {
        if self.0.is_empty() {
            return other;
        } else if other.0.is_empty() {
            return self.clone();
        }
        op_impl(
            self,
            &other,
            |coeff| coeff,
            |cur_coeff, other_coeff| cur_coeff + &other_coeff,
        )
    }
}

impl<'a, E: PairingEngine> Add<&'a LinearCombination<E>> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(self, other: &'a LinearCombination<E>) -> LinearCombination<E> {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            return other.clone();
        }
        op_impl(
            &self,
            other,
            |coeff| coeff,
            |cur_coeff, other_coeff| cur_coeff + &other_coeff,
        )
    }
}

impl<E: PairingEngine> Add<LinearCombination<E>> for LinearCombination<E> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            return other;
        }
        op_impl(
            &self,
            &other,
            |coeff| coeff,
            |cur_coeff, other_coeff| cur_coeff + &other_coeff,
        )
    }
}

impl<E: PairingEngine> Sub<&LinearCombination<E>> for &LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, other: &LinearCombination<E>) -> LinearCombination<E> {
        if other.0.is_empty() {
            let cur = self.clone();
            return cur;
        } else if self.0.is_empty() {
            let mut other = other.clone();
            other.negate_in_place();
            return other;
        }

        op_impl(
            self,
            other,
            |coeff| -coeff,
            |cur_coeff, other_coeff| cur_coeff - &other_coeff,
        )
    }
}

impl<'a, E: PairingEngine> Sub<&'a LinearCombination<E>> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, other: &'a LinearCombination<E>) -> LinearCombination<E> {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            let mut other = other.clone();
            other.negate_in_place();
            return other;
        }
        op_impl(
            &self,
            other,
            |coeff| -coeff,
            |cur_coeff, other_coeff| cur_coeff - &other_coeff,
        )
    }
}

impl<E: PairingEngine> Sub<LinearCombination<E>> for &LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, mut other: LinearCombination<E>) -> LinearCombination<E> {
        if self.0.is_empty() {
            other.negate_in_place();
            return other;
        } else if other.0.is_empty() {
            return self.clone();
        }

        op_impl(
            self,
            &other,
            |coeff| -coeff,
            |cur_coeff, other_coeff| cur_coeff - &other_coeff,
        )
    }
}

impl<E: PairingEngine> Sub<LinearCombination<E>> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, mut other: LinearCombination<E>) -> LinearCombination<E> {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            other.negate_in_place();
            return other;
        }
        op_impl(
            &self,
            &other,
            |coeff| -coeff,
            |cur_coeff, other_coeff| cur_coeff - &other_coeff,
        )
    }
}

impl<E: PairingEngine> Add<(E::Fr, &LinearCombination<E>)> for &LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(self, (mul_coeff, other): (E::Fr, &LinearCombination<E>)) -> LinearCombination<E> {
        if other.0.is_empty() {
            return self.clone();
        } else if self.0.is_empty() {
            let mut other = other.clone();
            other.mul_assign(mul_coeff);
            return other;
        }
        op_impl(
            self,
            other,
            |coeff| mul_coeff * &coeff,
            |cur_coeff, other_coeff| cur_coeff + &(mul_coeff * &other_coeff),
        )
    }
}

impl<'a, E: PairingEngine> Add<(E::Fr, &'a LinearCombination<E>)> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(self, (mul_coeff, other): (E::Fr, &'a LinearCombination<E>)) -> LinearCombination<E> {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            let mut other = other.clone();
            other.mul_assign(mul_coeff);
            return other;
        }
        op_impl(
            &self,
            other,
            |coeff| mul_coeff * &coeff,
            |cur_coeff, other_coeff| cur_coeff + &(mul_coeff * &other_coeff),
        )
    }
}

impl<E: PairingEngine> Add<(E::Fr, LinearCombination<E>)> for &LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn add(self, (mul_coeff, mut other): (E::Fr, LinearCombination<E>)) -> LinearCombination<E> {
        if other.0.is_empty() {
            return self.clone();
        } else if self.0.is_empty() {
            other.mul_assign(mul_coeff);
            return other;
        }
        op_impl(
            self,
            &other,
            |coeff| mul_coeff * &coeff,
            |cur_coeff, other_coeff| cur_coeff + &(mul_coeff * &other_coeff),
        )
    }
}

impl<E: PairingEngine> Add<(E::Fr, Self)> for LinearCombination<E> {
    type Output = Self;

    fn add(self, (mul_coeff, other): (E::Fr, Self)) -> Self {
        if other.0.is_empty() {
            return self;
        } else if self.0.is_empty() {
            let mut other = other;
            other.mul_assign(mul_coeff);
            return other;
        }
        op_impl(
            &self,
            &other,
            |coeff| mul_coeff * &coeff,
            |cur_coeff, other_coeff| cur_coeff + &(mul_coeff * &other_coeff),
        )
    }
}

impl<E: PairingEngine> Sub<(E::Fr, &LinearCombination<E>)> for &LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, (coeff, other): (E::Fr, &LinearCombination<E>)) -> LinearCombination<E> {
        self + (-coeff, other)
    }
}

impl<'a, E: PairingEngine> Sub<(E::Fr, &'a LinearCombination<E>)> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, (coeff, other): (E::Fr, &'a LinearCombination<E>)) -> LinearCombination<E> {
        self + (-coeff, other)
    }
}

impl<E: PairingEngine> Sub<(E::Fr, LinearCombination<E>)> for &LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, (coeff, other): (E::Fr, LinearCombination<E>)) -> LinearCombination<E> {
        self + (-coeff, other)
    }
}

impl<'a, E: PairingEngine> Sub<(E::Fr, LinearCombination<E>)> for LinearCombination<E> {
    type Output = LinearCombination<E>;

    fn sub(self, (coeff, other): (E::Fr, LinearCombination<E>)) -> LinearCombination<E> {
        self + (-coeff, other)
    }
}

/// This is an error that could occur during circuit synthesis contexts,
/// such as CRS generation, proving or verification.
#[derive(Debug)]
pub enum SynthesisError {
    /// During synthesis, we lacked knowledge of a variable assignment.
    AssignmentMissing,
    /// During synthesis, we divided by zero.
    DivisionByZero,
    /// During synthesis, we constructed an unsatisfiable constraint system.
    Unsatisfiable,
    /// During synthesis, our polynomials ended up being too high of degree
    PolynomialDegreeTooLarge,
    /// During proof generation, we encountered an identity in the CRS
    UnexpectedIdentity,
    /// During proof generation, we encountered an I/O error with the CRS
    IoError(io::Error),
    /// During verification, our verifying key was malformed.
    MalformedVerifyingKey,
    /// During CRS generation, we observed an unconstrained auxiliary variable
    UnconstrainedVariable,
}

impl From<io::Error> for SynthesisError {
    fn from(e: io::Error) -> SynthesisError {
        SynthesisError::IoError(e)
    }
}

impl Error for SynthesisError {
    fn description(&self) -> &str {
        match *self {
            SynthesisError::AssignmentMissing => {
                "an assignment for a variable could not be computed"
            },
            SynthesisError::DivisionByZero => "division by zero",
            SynthesisError::Unsatisfiable => "unsatisfiable constraint system",
            SynthesisError::PolynomialDegreeTooLarge => "polynomial degree is too large",
            SynthesisError::UnexpectedIdentity => "encountered an identity element in the CRS",
            SynthesisError::IoError(_) => "encountered an I/O error",
            SynthesisError::MalformedVerifyingKey => "malformed verifying key",
            SynthesisError::UnconstrainedVariable => "auxiliary variable was unconstrained",
        }
    }
}

impl fmt::Display for SynthesisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        if let &SynthesisError::IoError(ref e) = self {
            write!(f, "I/O error: ")?;
            e.fmt(f)
        } else {
            write!(f, "{}", self.description())
        }
    }
}

/// Represents a constraint system which can have new variables
/// allocated and constrains between them formed.
pub trait ConstraintSystem<E: PairingEngine>: Sized {
    /// Represents the type of the "root" of this constraint system
    /// so that nested namespaces can minimize indirection.
    type Root: ConstraintSystem<E>;

    /// Return the "one" input variable
    fn one() -> Variable {
        Variable::new_unchecked(Index::Input(0))
    }

    /// Allocate a private variable in the constraint system. The provided
    /// function is used to determine the assignment of the variable. The
    /// given `annotation` function is invoked in testing contexts in order
    /// to derive a unique name for this variable in the current namespace.
    fn alloc<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>;

    /// Allocate a public variable in the constraint system. The provided
    /// function is used to determine the assignment of the variable.
    fn alloc_input<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>;

    /// Enforce that `A` * `B` = `C`. The `annotation` function is invoked in
    /// testing contexts in order to derive a unique name for the constraint
    /// in the current namespace.
    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>;

    /// Create a new (sub)namespace and enter into it. Not intended
    /// for downstream use; use `namespace` instead.
    fn push_namespace<NR, N>(&mut self, name_fn: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR;

    /// Exit out of the existing namespace. Not intended for
    /// downstream use; use `namespace` instead.
    fn pop_namespace(&mut self);

    /// Gets the "root" constraint system, bypassing the namespacing.
    /// Not intended for downstream use; use `namespace` instead.
    fn get_root(&mut self) -> &mut Self::Root;

    /// Begin a namespace for this constraint system.
    fn ns<'a, NR, N>(&'a mut self, name_fn: N) -> Namespace<'a, E, Self::Root>
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        self.get_root().push_namespace(name_fn);

        Namespace(self.get_root(), PhantomData)
    }

    // Output the number of constraints in the system.
    fn num_constraints(&self) -> usize;
}

/// This is a "namespaced" constraint system which borrows a constraint system
/// (pushing a namespace context) and, when dropped, pops out of the namespace
/// context.
pub struct Namespace<'a, E: PairingEngine, CS: ConstraintSystem<E>>(&'a mut CS, PhantomData<E>);

impl<E: PairingEngine, CS: ConstraintSystem<E>> ConstraintSystem<E> for Namespace<'_, E, CS> {
    type Root = CS::Root;

    #[inline]
    fn one() -> Variable {
        CS::one()
    }

    #[inline]
    fn alloc<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.0.alloc(annotation, f)
    }

    #[inline]
    fn alloc_input<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        self.0.alloc_input(annotation, f)
    }

    #[inline]
    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
    {
        self.0.enforce(annotation, a, b, c)
    }

    // Downstream users who use `namespace` will never interact with these
    // functions and they will never be invoked because the namespace is
    // never a root constraint system.

    #[inline]
    fn push_namespace<NR, N>(&mut self, _: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        panic!("only the root's push_namespace should be called");
    }

    #[inline]
    fn pop_namespace(&mut self) {
        panic!("only the root's pop_namespace should be called");
    }

    #[inline]
    fn get_root(&mut self) -> &mut Self::Root {
        self.0.get_root()
    }

    #[inline]
    fn num_constraints(&self) -> usize {
        self.0.num_constraints()
    }
}

impl<E: PairingEngine, CS: ConstraintSystem<E>> Drop for Namespace<'_, E, CS> {
    #[inline]
    fn drop(&mut self) {
        self.get_root().pop_namespace()
    }
}

/// Convenience implementation of ConstraintSystem<E> for mutable references to
/// constraint systems.
impl<E: PairingEngine, CS: ConstraintSystem<E>> ConstraintSystem<E> for &mut CS {
    type Root = CS::Root;

    #[inline]
    fn one() -> Variable {
        CS::one()
    }

    #[inline]
    fn alloc<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        (**self).alloc(annotation, f)
    }

    #[inline]
    fn alloc_input<F, A, AR>(&mut self, annotation: A, f: F) -> Result<Variable, SynthesisError>
    where
        F: FnOnce() -> Result<E::Fr, SynthesisError>,
        A: FnOnce() -> AR,
        AR: Into<String>,
    {
        (**self).alloc_input(annotation, f)
    }

    #[inline]
    fn enforce<A, AR, LA, LB, LC>(&mut self, annotation: A, a: LA, b: LB, c: LC)
    where
        A: FnOnce() -> AR,
        AR: Into<String>,
        LA: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LB: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
        LC: FnOnce(LinearCombination<E>) -> LinearCombination<E>,
    {
        (**self).enforce(annotation, a, b, c)
    }

    #[inline]
    fn push_namespace<NR, N>(&mut self, name_fn: N)
    where
        NR: Into<String>,
        N: FnOnce() -> NR,
    {
        (**self).push_namespace(name_fn)
    }

    #[inline]
    fn pop_namespace(&mut self) {
        (**self).pop_namespace()
    }

    #[inline]
    fn get_root(&mut self) -> &mut Self::Root {
        (**self).get_root()
    }

    #[inline]
    fn num_constraints(&self) -> usize {
        (**self).num_constraints()
    }
}
