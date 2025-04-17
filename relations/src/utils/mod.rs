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

/// A type alias for a hash builder that uses the foldhash algorithm.
pub type HashBuilder = foldhash::fast::SeedableRandomState;
/// A type alias for an index map that uses the foldhash algorithm.
pub type IndexMap<K, V> = indexmap::IndexMap<K, V, HashBuilder>;
