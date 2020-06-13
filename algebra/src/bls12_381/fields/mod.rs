#[cfg(any(feature = "bls12_381", feature = "ed_on_bls12_381"))]
pub mod fr;
#[cfg(any(feature = "bls12_381", feature = "ed_on_bls12_381"))]
pub use self::fr::*;

#[cfg(feature = "bls12_381")]
pub mod fq;
#[cfg(feature = "bls12_381")]
pub use self::fq::*;

#[cfg(feature = "bls12_381")]
pub mod fq2;
#[cfg(feature = "bls12_381")]
pub use self::fq2::*;

#[cfg(feature = "bls12_381")]
pub mod fq6;
#[cfg(feature = "bls12_381")]
pub use self::fq6::*;

#[cfg(feature = "bls12_381")]
pub mod fq12;
#[cfg(feature = "bls12_381")]
pub use self::fq12::*;

#[cfg(all(feature = "bls12_381", test))]
mod tests;
