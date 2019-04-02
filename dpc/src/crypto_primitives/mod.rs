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

#[derive(Debug, Fail)]
pub enum CryptoError {
    #[fail(display = "input length is wrong: {}", _0)]
    IncorrectInputLength(usize),
    #[fail(display = "Element is not prime order")]
    NotPrimeOrder,
}
