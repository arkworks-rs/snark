#[cfg(feature = "mnt6_298")]
mod curves;
#[cfg(any(feature = "mnt6_298", feature = "mnt4_298"))]
mod fields;

#[cfg(feature = "mnt6_298")]
pub use curves::*;
#[cfg(any(feature = "mnt6_298", feature = "mnt4_298"))]
pub use fields::*;
