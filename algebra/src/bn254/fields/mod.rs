#[cfg(any(feature = "bn254", feature = "ed_on_bn254"))]
pub mod fr;
#[cfg(any(feature = "bn254", feature = "ed_on_bn254"))]
pub use self::fr::*;

#[cfg(feature = "bn254")]
pub mod fq;
#[cfg(feature = "bn254")]
pub use self::fq::*;

#[cfg(feature = "bn254")]
pub mod fq2;
#[cfg(feature = "bn254")]
pub use self::fq2::*;

#[cfg(feature = "bn254")]
pub mod fq6;
#[cfg(feature = "bn254")]
pub use self::fq6::*;

#[cfg(feature = "bn254")]
pub mod fq12;
#[cfg(feature = "bn254")]
pub use self::fq12::*;

#[cfg(all(feature = "bn254", test))]
mod tests;
