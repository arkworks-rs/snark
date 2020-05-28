#[cfg(any(feature = "bw6_761", feature = "edwards_sw6"))]
pub mod fr;
#[cfg(any(feature = "bw6_761", feature = "edwards_sw6"))]
pub use self::fr::*;

#[cfg(feature = "bw6_761")]
pub mod fq;
#[cfg(feature = "bw6_761")]
pub use self::fq::*;

#[cfg(feature = "bw6_761")]
pub mod fq3;
#[cfg(feature = "bw6_761")]
pub use self::fq3::*;

#[cfg(feature = "bw6_761")]
pub mod fq6;
#[cfg(feature = "bw6_761")]
pub use self::fq6::*;

#[cfg(all(feature = "bw6_761", test))]
mod tests;
