#[macro_use]
extern crate bench_utils;

#[macro_use]
extern crate derivative;

pub mod commitment;
pub mod crh;
pub mod merkle_tree;
pub mod nizk;
pub mod prf;
pub mod signature;
pub mod vrf;

pub use self::{
    commitment::CommitmentScheme,
    crh::FixedLengthCRH,
    merkle_tree::{MerkleHashTree, MerkleTreePath},
    nizk::NIZK,
    prf::PRF,
    signature::SignatureScheme,
};

#[cfg(feature = "r1cs")]
pub use self::{
    commitment::CommitmentGadget, crh::FixedLengthCRHGadget,
    merkle_tree::constraints::MerkleTreePathGadget, nizk::NIZKVerifierGadget, prf::PRFGadget,
    signature::SigRandomizePkGadget,
};

pub type Error = Box<dyn std::error::Error>;

#[derive(Debug)]
pub enum CryptoError {
    IncorrectInputLength(usize),
    NotPrimeOrder,
    FailedVerification,
}

impl std::fmt::Display for CryptoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let msg = match self {
            CryptoError::IncorrectInputLength(len) => format!("input length is wrong: {}", len),
            CryptoError::NotPrimeOrder => "element is not prime order".to_owned(),
            CryptoError::FailedVerification => "verification failed".to_owned(),
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
pub(crate) fn compute_truncation_size(modulus_from: i32, modulus_to: i32) -> usize {
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
pub(crate) fn compute_bytes_truncation_size<F: PrimeField>() -> usize {
    let bigint_bytes = (F::Params::MODULUS_BITS + F::Params::REPR_SHAVE_BITS)/8;
    let safe_bytes = F::Params::CAPACITY/8;
    (bigint_bytes - safe_bytes) as usize
}