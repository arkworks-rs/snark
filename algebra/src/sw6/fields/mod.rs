#[cfg(any(feature = "sw6", feature = "edwards_sw6"))]
pub mod fr;
#[cfg(any(feature = "sw6", feature = "edwards_sw6"))]
pub use self::fr::*;

#[cfg(feature = "sw6")]
pub mod fq;
#[cfg(feature = "sw6")]
pub use self::fq::*;

#[cfg(feature = "sw6")]
pub mod fq3;
#[cfg(feature = "sw6")]
pub use self::fq3::*;

#[cfg(feature = "sw6")]
pub mod fq6;
#[cfg(feature = "sw6")]
pub use self::fq6::*;

#[cfg(all(feature = "sw6", test))]
mod tests;
