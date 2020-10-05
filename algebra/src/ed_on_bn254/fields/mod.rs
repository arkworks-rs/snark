pub mod fq;
pub mod fr;

pub use fq::*;
pub use fr::*;

#[cfg(all(feature = "ed_on_bn254", test, feature = "prime_fields"))]
mod tests;
