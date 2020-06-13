#[cfg(any(feature = "mnt6_298", feature = "mnt4_298", feature = "ed_on_mnt4_298"))]
pub mod fr;
#[cfg(any(feature = "mnt6_298", feature = "mnt4_298", feature = "ed_on_mnt4_298"))]
pub use self::fr::*;

#[cfg(any(feature = "mnt6_298", feature = "mnt4_298", feature = "ed_on_mnt4_298"))]
pub mod fq;
#[cfg(any(feature = "mnt6_298", feature = "mnt4_298", feature = "ed_on_mnt4_298"))]
pub use self::fq::*;

#[cfg(feature = "mnt6_298")]
pub mod fq3;
#[cfg(feature = "mnt6_298")]
pub use self::fq3::*;

#[cfg(feature = "mnt6_298")]
pub mod fq6;
#[cfg(feature = "mnt6_298")]
pub use self::fq6::*;

#[cfg(all(feature = "mnt6_298", test))]
#[cfg(test)]
mod tests;
