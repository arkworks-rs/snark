#[macro_use]
pub mod scalar_mul;
pub use scalar_mul::*;

#[cfg(not(feature = "cuda"))]
pub mod accel_dummy;
