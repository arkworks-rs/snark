#[cfg(feature = "mnt6_753")]
mod curves;
#[cfg(any(feature = "mnt6_753", feature = "mnt4_753"))]
mod fields;

#[cfg(feature = "mnt6_753")]
pub use curves::*;
#[cfg(any(feature = "mnt6_753", feature = "mnt4_753"))]
pub use fields::*;
