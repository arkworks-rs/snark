#[cfg(any(feature = "bls12_377", feature = "edwards_bls12"))]
pub mod fr;
#[cfg(any(feature = "bls12_377", feature = "edwards_bls12"))]
pub use self::fr::*;

#[cfg(any(feature = "bls12_377", feature = "sw6", feature = "edwards_sw6", feature = "bw6_761"))]
pub mod fq;
#[cfg(any(feature = "bls12_377", feature = "sw6", feature = "edwards_sw6", feature = "bw6_761"))]
pub use self::fq::*;

#[cfg(feature = "bls12_377")]
pub mod fq2;
#[cfg(feature = "bls12_377")]
pub use self::fq2::*;

#[cfg(feature = "bls12_377")]
pub mod fq6;
#[cfg(feature = "bls12_377")]
pub use self::fq6::*;

#[cfg(feature = "bls12_377")]
pub mod fq12;
#[cfg(feature = "bls12_377")]
pub use self::fq12::*;

#[cfg(all(feature = "bls12_377", test))]
mod tests;
