use core::fmt;

/// This is an error that could occur during circuit synthesis.
#[derive(PartialEq, Eq, Clone, Copy, Debug)]
pub enum SynthesisError {
    /// During synthesis, we tried to allocate a variable when
    /// `ConstraintSystemRef` was `None`.
    MissingCS,
    /// During synthesis, we lacked knowledge of a variable assignment.
    AssignmentMissing,
    /// During synthesis, we divided by zero.
    DivisionByZero,
    /// During synthesis, we constructed an unsatisfiable constraint system.
    Unsatisfiable,
    /// During synthesis, our polynomials ended up being too high of degree
    PolynomialDegreeTooLarge,
    /// The string does not match to any predicate
    PredicateNotFound,
    /// The predicate expects a different arity
    ArityMismatch,
}

impl ark_std::error::Error for SynthesisError {}

impl fmt::Display for SynthesisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        match self {
            SynthesisError::MissingCS => write!(f, "the constraint system was `None`"),
            SynthesisError::AssignmentMissing => write!(f, "assignment couldn't be computed"),
            SynthesisError::DivisionByZero => write!(f, "division by zero"),
            SynthesisError::Unsatisfiable => write!(f, "unsatisfiable constraint system"),
            SynthesisError::PolynomialDegreeTooLarge => write!(f, "polynomial degree too large"),
            SynthesisError::ArityMismatch => write!(f, "predicate arity doesn't match input"),
            SynthesisError::PredicateNotFound => write!(f, "predicate was not found"),
        }
    }
}
