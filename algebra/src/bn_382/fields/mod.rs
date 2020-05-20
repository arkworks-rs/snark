#[cfg(feature = "bn_382")]
pub mod fp;
#[cfg(feature = "bn_382")]
pub use self::fp::*;

#[cfg(feature = "bn_382")]
pub mod fq;
#[cfg(feature = "bn_382")]
pub use self::fq::*;

#[cfg(feature = "bn_382")]
pub mod fq2;
#[cfg(feature = "bn_382")]
pub use self::fq2::*;

#[cfg(feature = "bn_382")]
pub mod fq6;
#[cfg(feature = "bn_382")]
pub use self::fq6::*;

#[cfg(feature = "bn_382")]
pub mod fq12;
#[cfg(feature = "bn_382")]
pub use self::fq12::*;

#[cfg(all(feature = "bn_382", test))]
mod tests;

