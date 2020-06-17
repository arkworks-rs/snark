#[cfg(feature = "bn256")]
pub mod fr;
#[cfg(feature = "bn256")]
pub use self::fr::*;

#[cfg(feature = "bn256")]
pub mod fq;
#[cfg(feature = "bn256")]
pub use self::fq::*;

#[cfg(feature = "bn256")]
pub mod fq2;
#[cfg(feature = "bn256")]
pub use self::fq2::*;

#[cfg(feature = "bn256")]
pub mod fq6;
#[cfg(feature = "bn256")]
pub use self::fq6::*;

#[cfg(feature = "bn256")]
pub mod fq12;
#[cfg(feature = "bn256")]
pub use self::fq12::*;

#[cfg(all(feature = "bn256", test))]
mod tests;
