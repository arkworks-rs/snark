use std::io;

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
    /// A generic error, not classifiable in any of the above categories
    Other(String),
}

impl From<io::Error> for SynthesisError {
    fn from(e: io::Error) -> SynthesisError {
        SynthesisError::IoError(e)
    }
}

impl From<Box<dyn std::error::Error>> for SynthesisError {
    fn from(e: Box<dyn std::error::Error>) -> SynthesisError {
        SynthesisError::Other(e.to_string())
    }
}

impl std::fmt::Display for SynthesisError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SynthesisError::AssignmentMissing => {
                write!(f, "an assignment for a variable could not be computed")
            }
            SynthesisError::DivisionByZero => write!(f, "division by zero"),
            SynthesisError::Unsatisfiable => write!(f, "unsatisfiable constraint system"),
            SynthesisError::PolynomialDegreeTooLarge => write!(f, "polynomial degree is too large"),
            SynthesisError::UnexpectedIdentity => {
                write!(f, "encountered an identity element in the CRS")
            }
            SynthesisError::IoError(e) => write!(f, "{:?}", e),
            SynthesisError::MalformedVerifyingKey => write!(f, "malformed verifying key"),
            SynthesisError::UnconstrainedVariable => {
                write!(f, "auxiliary variable was unconstrained")
            }
            SynthesisError::Other(e) => write!(f, "{:?}", e),
        }
    }
}

impl std::error::Error for SynthesisError {}
