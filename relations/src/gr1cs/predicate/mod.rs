//! This module contains the implementation of a general  predicate, defined in https://eprint.iacr.org/2024/1245
//! A  predicate is a function from t (arity) variables to a boolean
//! variable A predicate can be as simple as f(a,b,c)=a.b-c=0 or as
//! complex as a lookup table

pub mod polynomial_constraint;

use super::{Constraint, ConstraintSystem, LcIndex, Matrix};
use crate::utils::{error::SynthesisError::ArityMismatch, variable::Variable::SymbolicLc};
use ark_ff::Field;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, SerializationError};
use ark_std::{io::Write, vec::Vec};
use polynomial_constraint::PolynomialPredicate;

/// GR1CS can potentially support different types of predicates
/// For now, we only support polynomial predicates
/// In the future, we can add other types of predicates, e.g. lookup table
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum PredicateType<F: Field> {
    /// A polynomial local predicate. This is the most common predicate that
    /// captures high-degree custome gates
    Polynomial(PolynomialPredicate<F>),
    // Add other predicates in the future, e.g. lookup table
}

impl<F: Field> ark_serialize::Valid for PredicateType<F> {
    fn check(&self) -> Result<(), SerializationError> {
        match self {
            PredicateType::Polynomial(p) => p.check(),
        }
    }
}
impl<F: Field> CanonicalDeserialize for PredicateType<F> {
    fn deserialize_with_mode<R: ark_serialize::Read>(
        reader: R,
        compress: Compress,
        should_validate: ark_serialize::Validate,
    ) -> Result<Self, SerializationError> {
        let predicate_type =
            PolynomialPredicate::<F>::deserialize_with_mode(reader, compress, should_validate)?;
        Ok(PredicateType::Polynomial(predicate_type))
    }
}

impl<F: Field> CanonicalSerialize for PredicateType<F> {
    fn serialize_with_mode<W: Write>(
        &self,
        writer: W,
        compress: Compress,
    ) -> Result<(), SerializationError> {
        match self {
            PredicateType::Polynomial(p) => p.serialize_with_mode(writer, compress),
        }
    }
    fn serialized_size(&self, compress: Compress) -> usize {
        match self {
            PredicateType::Polynomial(p) => p.serialized_size(compress),
        }
    }
}

impl<F: Field> PredicateType<F> {
    fn is_satisfied(&self, variables: &[F]) -> bool {
        match self {
            PredicateType::Polynomial(p) => p.is_satisfied(variables),
            // TODO: Add other predicates in the future, e.g. lookup table
        }
    }

    fn arity(&self) -> usize {
        match self {
            PredicateType::Polynomial(p) => p.arity(),
            // TODO: Add other predicates in the future, e.g. lookup table
        }
    }
}

/// A constraint system that enforces a predicate
#[derive(Debug, Clone)]
pub struct PredicateConstraintSystem<F: Field> {
    /// The list of linear combinations for each arguments of the predicate
    /// The length of this list is equal to the arity of the predicate
    /// Each element in the list has size equal to the number of constraints
    argument_lcs: Vec<Vec<LcIndex>>,

    /// The number of constraints enforced by this predicate
    num_constraints: usize,

    /// The type of the predicate enforced by this constraint system.  
    predicate_type: PredicateType<F>,
}

impl<F: Field> PredicateConstraintSystem<F> {
    /// Create a new predicate constraint system with a specific predicate
    fn new(predicate_type: PredicateType<F>) -> Self {
        Self {
            argument_lcs: vec![Vec::new(); predicate_type.arity()],
            predicate_type,
            num_constraints: 0,
        }
    }

    /// Create new polynomial predicate constraint system
    pub fn new_polynomial_predicate(arity: usize, terms: Vec<(F, Vec<(usize, usize)>)>) -> Self {
        Self::new(PredicateType::Polynomial(PolynomialPredicate::new(
            arity, terms,
        )))
    }

    /// creates an R1CS predicate which is a special kind of polynomial
    /// predicate
    pub fn new_r1cs_predicate() -> crate::utils::Result<Self> {
        Ok(Self::new_polynomial_predicate(
            3,
            vec![
                (F::from(1u8), vec![(0, 1), (1, 1)]),
                (F::from(-1i8), vec![(2, 1)]),
            ],
        ))
    }

    /// Get the arity of the  predicate in this predicate constraint system
    pub fn get_arity(&self) -> usize {
        self.predicate_type.arity()
    }

    /// Get the number of constraints enforced by this predicate
    pub fn num_constraints(&self) -> usize {
        self.num_constraints
    }

    /// Get the vector of constrints enforced by this predicate
    /// Each constraint is a list of linear combinations with size equal to the
    /// arity
    pub fn get_constraints(&self) -> &Vec<Constraint> {
        &self.argument_lcs
    }

    /// Get a reference to the  predicate in this predicate constraint
    /// system
    pub fn get_predicate_type(&self) -> &PredicateType<F> {
        &self.predicate_type
    }

    /// Enforce a constraint in this predicate constraint system
    /// The constraint is a list of linear combinations with size equal to the
    /// arity
    pub fn enforce_constraint(
        &mut self,
        constraint: impl IntoIterator<Item = LcIndex>,
    ) -> crate::utils::Result<()> {
        let mut arity = 0;
        constraint
            .into_iter()
            .zip(&mut self.argument_lcs)
            .for_each(|(lc_index, arg_lc)| {
                arity += 1;
                arg_lc.push(lc_index);
            });
        if arity != self.get_arity() {
            return Err(ArityMismatch);
        }

        self.num_constraints += 1;
        Ok(())
    }

    fn iter_constraints(&self) -> impl Iterator<Item = Constraint> + '_ {
        // Transpose the `argument_lcs` to iterate over constraints
        let num_constraints = self.num_constraints;

        (0..num_constraints).map(move |i| {
            self.argument_lcs
                .iter()
                .map(|lc_s| lc_s[i])
                .collect::<Vec<LcIndex>>()
        })
    }

    /// Check if the constraints enforced by this predicate are satisfied
    /// i.e. L(x_1, x_2, ..., x_n) = 0
    pub fn which_constraint_is_unsatisfied(&self, cs: &ConstraintSystem<F>) -> Option<usize> {
        for (i, constraint) in self.iter_constraints().enumerate() {
            let variables: Vec<F> = constraint
                .iter()
                .map(|lc_index| cs.assigned_value(SymbolicLc(*lc_index)).unwrap())
                .collect();
            let result = self.predicate_type.is_satisfied(&variables);
            if result {
                return Some(i);
            }
        }
        None
    }

    /// Create the set of matrices for this predicate constraint system
    pub fn to_matrices(&self, cs: &ConstraintSystem<F>) -> Vec<Matrix<F>> {
        let mut matrices: Vec<Matrix<F>> = vec![Vec::new(); self.get_arity()];
        for constraint in self.iter_constraints() {
            for (matrix_ind, lc_index) in constraint.iter().enumerate() {
                let lc = cs.get_lc(*lc_index).unwrap();
                let row = cs.make_row(&lc);
                matrices[matrix_ind].push(row);
            }
        }
        matrices
    }
}
