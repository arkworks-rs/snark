pub mod commitment;
pub mod crh;
pub mod mht;
pub mod nizk;
pub mod prf;
pub mod signature;

pub use self::{
    commitment::CommitmentScheme,
    crh::FixedLengthCRH,
    mht::{HashMembershipProof, MerkleHashTree},
    nizk::NIZK,
    prf::PRF,
    signature::SignatureScheme,
};

#[derive(Debug)]
pub enum CryptoError {
    IncorrectInputLength(usize),
    NotPrimeOrder,
}

impl std::fmt::Display for CryptoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let msg = match self {
            CryptoError::IncorrectInputLength(len) => format!("input length is wrong: {}", len),
            CryptoError::NotPrimeOrder => "element is not prime order".to_owned(),
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
