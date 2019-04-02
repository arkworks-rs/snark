use crate::dpc::Transaction;
use algebra::bytes::ToBytes;
use failure::Error;
use rand::Rng;

pub mod ideal_ledger;
pub use self::ideal_ledger::*;

pub trait LedgerDigest: Clone + ToBytes + Default + Eq {
    type Parameters: Clone + Default;
}
pub trait LedgerWitness<D: LedgerDigest>: Clone + Default {
    fn dummy_witness() -> Self;
}

pub trait Ledger {
    type Parameters: Clone + Default;
    type LedgerStateDigest: LedgerDigest<Parameters = Self::Parameters>;

    type Commitment;
    type CommWitness: Clone + LedgerWitness<Self::LedgerStateDigest>;

    type SerialNumber;
    type SnWitness: LedgerWitness<Self::LedgerStateDigest>;

    type Memo;
    type MemoWitness: LedgerWitness<Self::LedgerStateDigest>;

    type Transaction: Transaction;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Parameters, Error>;

    /// Creates an empty ledger
    fn new(
        parameters: Self::Parameters,
        dummy_cm: Self::Commitment,
        dummy_sn: Self::SerialNumber,
        dummy_memo: Self::Memo,
    ) -> Self;

    /// Return the current number of transactions on the ledger.
    fn len(&self) -> usize;

    /// Return the parameters used to construct the ledger data structure.
    fn parameters(&self) -> &Self::Parameters;

    /// Append a (valid) transaction tx to the ledger.
    fn push(&mut self, transaction: Self::Transaction) -> Result<(), Error>;

    /// Return a short digest of the current state of the transaction set data
    /// structure.
    fn digest(&self) -> Option<Self::LedgerStateDigest>;

    /// Check that st_{ts} is a valid digest for some (past) ledger state.
    fn validate_digest(&self, digest: &Self::LedgerStateDigest) -> bool;

    fn contains_cm(&self, cm: &Self::Commitment) -> bool;
    fn contains_sn(&self, sn: &Self::SerialNumber) -> bool;
    fn contains_memo(&self, memo: &Self::Memo) -> bool;

    fn prove_cm(&self, cm: &Self::Commitment) -> Result<Self::CommWitness, Error>;
    fn prove_sn(&self, sn: &Self::SerialNumber) -> Result<Self::SnWitness, Error>;
    fn prove_memo(&self, memo: &Self::Memo) -> Result<Self::MemoWitness, Error>;

    fn verify_cm(
        params: &Self::Parameters,
        digest: &Self::LedgerStateDigest,
        cm: &Self::Commitment,
        witness: &Self::CommWitness,
    ) -> bool;

    fn verify_sn(
        params: &Self::Parameters,
        digest: &Self::LedgerStateDigest,
        sn: &Self::SerialNumber,
        witness: &Self::SnWitness,
    ) -> bool;

    fn verify_memo(
        params: &Self::Parameters,
        digest: &Self::LedgerStateDigest,
        memo: &Self::Memo,
        witness: &Self::MemoWitness,
    ) -> bool;
}
