#[cfg(not(feature = "std"))]
use alloc::vec::Vec;
pub mod error {
    pub type Result<T> = T;
}

pub struct Context {}

pub type DeviceMemory<T> = Vec<T>;
