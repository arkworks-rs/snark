use error::SynthesisError;

/// Different errors and error messages used in the library.
pub mod error;

/// Tools to work with linear combinations.
pub mod linear_combination;

/// Tools to work with matrices.
pub mod matrix;

/// Tools to work with variables.
pub mod variable;

/// A result type specialized to `SynthesisError`.
pub type Result<T> = core::result::Result<T, SynthesisError>;




