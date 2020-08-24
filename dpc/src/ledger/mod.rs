use crate::{dpc::Transaction, Error};
pub use crypto_primitives::merkle_tree;
use rand::Rng;

pub mod ideal_ledger;
pub use self::ideal_ledger::*;

pub trait Ledger {
    type Parameters: merkle_tree::Config;

    type Commitment;
    type SerialNumber;
    type Memo;

    type Transaction: Transaction;

    fn setup<R: Rng>(rng: &mut R) -> Result<merkle_tree::Parameters<Self::Parameters>, Error>;

    /// Creates an empty ledger
    fn new(
        parameters: merkle_tree::Parameters<Self::Parameters>,
        dummy_cm: Self::Commitment,
        dummy_sn: Self::SerialNumber,
        dummy_memo: Self::Memo,
    ) -> Self;

    /// Return the current number of transactions on the ledger.
    fn len(&self) -> usize;

    /// Return the parameters used to construct the ledger data structure.
    fn parameters(&self) -> &merkle_tree::Parameters<Self::Parameters>;

    /// Append a (valid) transaction tx to the ledger.
    fn push(&mut self, transaction: Self::Transaction) -> Result<(), Error>;

    /// Return a short digest of the current state of the transaction set data
    /// structure.
    fn digest(&self) -> Option<merkle_tree::Digest<Self::Parameters>>;

    /// Check that st_{ts} is a valid digest for some (past) ledger state.
    fn validate_digest(&self, digest: &merkle_tree::Digest<Self::Parameters>) -> bool;

    fn contains_cm(&self, cm: &Self::Commitment) -> bool;
    fn contains_sn(&self, sn: &Self::SerialNumber) -> bool;
    fn contains_memo(&self, memo: &Self::Memo) -> bool;

    fn prove_cm(&self, cm: &Self::Commitment)
        -> Result<merkle_tree::Path<Self::Parameters>, Error>;
    fn prove_sn(
        &self,
        sn: &Self::SerialNumber,
    ) -> Result<merkle_tree::Path<Self::Parameters>, Error>;
    fn prove_memo(&self, memo: &Self::Memo) -> Result<merkle_tree::Path<Self::Parameters>, Error>;

    fn verify_cm(
        parameters: &merkle_tree::Parameters<Self::Parameters>,
        digest: &merkle_tree::Digest<Self::Parameters>,
        cm: &Self::Commitment,
        witness: &merkle_tree::Path<Self::Parameters>,
    ) -> bool;

    fn verify_sn(
        parameters: &merkle_tree::Parameters<Self::Parameters>,
        digest: &merkle_tree::Digest<Self::Parameters>,
        sn: &Self::SerialNumber,
        witness: &merkle_tree::Path<Self::Parameters>,
    ) -> bool;

    fn verify_memo(
        parameters: &merkle_tree::Parameters<Self::Parameters>,
        digest: &merkle_tree::Digest<Self::Parameters>,
        memo: &Self::Memo,
        witness: &merkle_tree::Path<Self::Parameters>,
    ) -> bool;
}
