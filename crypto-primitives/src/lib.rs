#![cfg_attr(not(feature = "std"), no_std)]

#[macro_use]
extern crate bench_utils;

#[macro_use]
extern crate derivative;

#[macro_use]
extern crate alloc;

#[cfg(not(feature = "std"))]
pub(crate) use alloc::{borrow::ToOwned, boxed::Box, vec::Vec};

#[cfg(feature = "std")]
pub(crate) use std::{borrow::ToOwned, boxed::Box, vec::Vec};

pub mod commitment;
pub mod crh;
pub mod merkle_tree;
pub mod nizk;
pub mod prf;
pub mod signature;

pub use self::{
    commitment::CommitmentScheme,
    crh::FixedLengthCRH,
    merkle_tree::{MerkleTree, Path},
    nizk::NIZK,
    prf::PRF,
    signature::SignatureScheme,
};

#[cfg(feature = "r1cs")]
#[cfg(feature = "r1cs")]
pub use self::{
    commitment::CommitmentGadget, crh::FixedLengthCRHGadget, merkle_tree::constraints::PathVar,
    nizk::NIZKVerifierGadget, prf::PRFGadget, signature::SigRandomizePkGadget,
};

pub type Error = Box<dyn algebra_core::Error>;

#[derive(Debug)]
pub enum CryptoError {
    IncorrectInputLength(usize),
    NotPrimeOrder,
}

impl core::fmt::Display for CryptoError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let msg = match self {
            CryptoError::IncorrectInputLength(len) => format!("input length is wrong: {}", len),
            CryptoError::NotPrimeOrder => "element is not prime order".to_owned(),
        };
        write!(f, "{}", msg)
    }
}

impl algebra_core::Error for CryptoError {}
