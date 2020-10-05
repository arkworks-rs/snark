pub mod fq;
pub mod fr;

pub use fq::*;
pub use fr::*;

#[cfg(all(feature = "ed_on_mnt4_753", test, feature = "prime_fields"))]
mod tests;
