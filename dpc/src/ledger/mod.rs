use crypto_primitives::{FixedLengthCRH, mht::{MHTParameters, HashMembershipProof}};
use crate::dpc::Transaction;
use crate::Error;
use rand::Rng;

pub mod ideal_ledger;
pub use self::ideal_ledger::*;

pub type MHTParams<P> = <<P as MHTParameters>::H as FixedLengthCRH>::Parameters;
pub type MHTDigest<P> = <<P as MHTParameters>::H as FixedLengthCRH>::Output;

pub trait Ledger {
    type Parameters: MHTParameters;

    type Commitment;
    type SerialNumber;
    type Memo;

    type Transaction: Transaction;

    fn setup<R: Rng>(rng: &mut R) -> Result<MHTParams<Self::Parameters>, Error>;

    /// Creates an empty ledger
    fn new(
        parameters: MHTParams<Self::Parameters>,
        dummy_cm: Self::Commitment,
        dummy_sn: Self::SerialNumber,
        dummy_memo: Self::Memo,
    ) -> Self;

    /// Return the current number of transactions on the ledger.
    fn len(&self) -> usize;

    /// Return the parameters used to construct the ledger data structure.
    fn parameters(&self) -> &MHTParams<Self::Parameters>;

    /// Append a (valid) transaction tx to the ledger.
    fn push(&mut self, transaction: Self::Transaction) -> Result<(), Error>;

    /// Return a short digest of the current state of the transaction set data
    /// structure.
    fn digest(&self) -> Option<MHTDigest<Self::Parameters>>;

    /// Check that st_{ts} is a valid digest for some (past) ledger state.
    fn validate_digest(&self, digest: &MHTDigest<Self::Parameters>) -> bool;

    fn contains_cm(&self, cm: &Self::Commitment) -> bool;
    fn contains_sn(&self, sn: &Self::SerialNumber) -> bool;
    fn contains_memo(&self, memo: &Self::Memo) -> bool;

    fn prove_cm(&self, cm: &Self::Commitment) -> Result<HashMembershipProof<Self::Parameters>, Error>;
    fn prove_sn(&self, sn: &Self::SerialNumber) -> Result<HashMembershipProof<Self::Parameters>, Error>;
    fn prove_memo(&self, memo: &Self::Memo) -> Result<HashMembershipProof<Self::Parameters>, Error>;

    fn verify_cm(
        parameters: &MHTParams<Self::Parameters>,
        digest: &MHTDigest<Self::Parameters>,
        cm: &Self::Commitment,
        witness: &HashMembershipProof<Self::Parameters>,
    ) -> bool;

    fn verify_sn(
        parameters: &MHTParams<Self::Parameters>,
        digest: &MHTDigest<Self::Parameters>,
        sn: &Self::SerialNumber,
        witness: &HashMembershipProof<Self::Parameters>,
    ) -> bool;

    fn verify_memo(
        parameters: &MHTParams<Self::Parameters>,
        digest: &MHTDigest<Self::Parameters>,
        memo: &Self::Memo,
        witness: &HashMembershipProof<Self::Parameters>,
    ) -> bool;
}
