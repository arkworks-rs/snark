pub use crypto_primitives::merkle_tree::*;
use crate::dpc::Transaction;
use crate::Error;
use rand::Rng;

pub mod ideal_ledger;
pub use self::ideal_ledger::*;

pub trait Ledger {
    type Parameters: MerkleTreeConfig;

    type Commitment;
    type SerialNumber;
    type Memo;

    type Transaction: Transaction;

    fn setup<R: Rng>(rng: &mut R) -> Result<MerkleTreeParams<Self::Parameters>, Error>;

    /// Creates an empty ledger
    fn new(
        parameters: MerkleTreeParams<Self::Parameters>,
        dummy_cm: Self::Commitment,
        dummy_sn: Self::SerialNumber,
        dummy_memo: Self::Memo,
    ) -> Self;

    /// Return the current number of transactions on the ledger.
    fn len(&self) -> usize;

    /// Return the parameters used to construct the ledger data structure.
    fn parameters(&self) -> &MerkleTreeParams<Self::Parameters>;

    /// Append a (valid) transaction tx to the ledger.
    fn push(&mut self, transaction: Self::Transaction) -> Result<(), Error>;

    /// Return a short digest of the current state of the transaction set data
    /// structure.
    fn digest(&self) -> Option<MerkleTreeDigest<Self::Parameters>>;

    /// Check that st_{ts} is a valid digest for some (past) ledger state.
    fn validate_digest(&self, digest: &MerkleTreeDigest<Self::Parameters>) -> bool;

    fn contains_cm(&self, cm: &Self::Commitment) -> bool;
    fn contains_sn(&self, sn: &Self::SerialNumber) -> bool;
    fn contains_memo(&self, memo: &Self::Memo) -> bool;

    fn prove_cm(&self, cm: &Self::Commitment) -> Result<MerkleTreePath<Self::Parameters>, Error>;
    fn prove_sn(&self, sn: &Self::SerialNumber) -> Result<MerkleTreePath<Self::Parameters>, Error>;
    fn prove_memo(&self, memo: &Self::Memo) -> Result<MerkleTreePath<Self::Parameters>, Error>;

    fn verify_cm(
        parameters: &MerkleTreeParams<Self::Parameters>,
        digest: &MerkleTreeDigest<Self::Parameters>,
        cm: &Self::Commitment,
        witness: &MerkleTreePath<Self::Parameters>,
    ) -> bool;

    fn verify_sn(
        parameters: &MerkleTreeParams<Self::Parameters>,
        digest: &MerkleTreeDigest<Self::Parameters>,
        sn: &Self::SerialNumber,
        witness: &MerkleTreePath<Self::Parameters>,
    ) -> bool;

    fn verify_memo(
        parameters: &MerkleTreeParams<Self::Parameters>,
        digest: &MerkleTreeDigest<Self::Parameters>,
        memo: &Self::Memo,
        witness: &MerkleTreePath<Self::Parameters>,
    ) -> bool;
}
