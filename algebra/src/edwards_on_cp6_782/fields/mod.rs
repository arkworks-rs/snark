pub mod fq;
pub mod fr;

pub use fq::*;
pub use fr::*;

#[cfg(all(feature = "edwards_on_cp6_782", test))]
mod tests;
