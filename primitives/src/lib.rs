#[macro_use]
extern crate bench_utils;

#[macro_use]
extern crate derivative;

pub mod crh;
pub use self::crh::*;

#[cfg(feature = "commitment")]
pub mod commitment;
#[cfg(feature = "commitment")]
pub use self::commitment::*;

#[cfg(feature = "merkle_tree")]
pub mod merkle_tree;
#[cfg(feature = "merkle_tree")]
pub use self::merkle_tree::*;

#[cfg(feature = "prf")]
pub mod prf;
#[cfg(feature = "prf")]
pub use self::prf::*;

#[cfg(feature = "signature")]
pub mod signature;
#[cfg(feature = "signature")]
pub use self::signature::*;

#[cfg(feature = "vrf")]
pub mod vrf;
#[cfg(feature = "vrf")]
pub use self::vrf::*;


pub type Error = Box<dyn std::error::Error>;

#[derive(Debug)]
pub enum CryptoError {
    IncorrectInputLength(String, usize),
    InvalidElement(String),
    NotPrimeOrder(String),
    FailedVerification,
    InitializationError(String),
    HashingError(String),
    Other(String),
}

impl std::fmt::Display for CryptoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let msg = match self {
            CryptoError::IncorrectInputLength(elem, len) => format!("{} length is wrong: {}", elem, len),
            CryptoError::InvalidElement(elem) => format!("{} is invalid", elem),
            CryptoError::NotPrimeOrder(elem) => format!("element {} is not prime order", elem),
            CryptoError::FailedVerification => "verification failed".to_owned(),
            CryptoError::InitializationError(message) => format!("{}", message),
            CryptoError::HashingError(message) => format!("Failed to compute the hash: {}", message),
            CryptoError::Other(message) => format!("{}", message),
        };
        write!(f, "{}", msg)
    }
}

impl std::error::Error for CryptoError {
    #[inline]
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

/// Return the number of leading bits to skip in a field element belonging to a field
/// 'from' having `modulus_from` bits in order to safely convert it into a field element
/// belonging to a field 'to' having `modulus_to` bits.
pub fn compute_truncation_size(modulus_from: i32, modulus_to: i32) -> usize {
    (match modulus_from - modulus_to {
        moduli_diff if moduli_diff > 0 => moduli_diff + 1,
        moduli_diff if moduli_diff == 0 => 1,
        moduli_diff if moduli_diff < 0 => 0,
        _ => unreachable!(),
    }) as usize
}

use algebra::{
    PrimeField, FpParameters,
};

/// Return the number of bytes to skip in a little-endian byte order representation
/// of a field element belonging to field `F`.
#[allow(dead_code)]
pub fn compute_bytes_truncation_size<F: PrimeField>() -> usize {
    let bigint_bytes = (F::Params::MODULUS_BITS + F::Params::REPR_SHAVE_BITS)/8;
    let safe_bytes = F::Params::CAPACITY/8;
    (bigint_bytes - safe_bytes) as usize
}

pub fn bytes_to_bits(bytes: &[u8]) -> Vec<bool> {
    let mut bits = Vec::with_capacity(bytes.len() * 8);
    for byte in bytes {
        for i in 0..8 {
            let bit = (*byte >> i) & 1;
            bits.push(bit == 1)
        }
    }
    bits
}