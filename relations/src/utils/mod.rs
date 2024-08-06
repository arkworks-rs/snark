use ark_std::vec::Vec;
use error::SynthesisError;
use linear_combination::LcIndex;

pub mod error;
pub mod linear_combination;
pub mod matrix;
pub mod variable;

/// A result type specialized to `SynthesisError`.
pub type Result<T> = core::result::Result<T, SynthesisError>;
