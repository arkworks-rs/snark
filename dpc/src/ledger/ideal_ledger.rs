use crate::Error;
use rand::Rng;
use std::{
    collections::{HashMap, HashSet},
    fmt,
    hash::Hash,
    io::{Result as IoResult, Write},
    rc::Rc,
};

use crate::{
    crypto_primitives::{FixedLengthCRH, HashMembershipProof, MerkleHashTree},
    dpc::Transaction,
    ledger::{Ledger, LedgerDigest, LedgerWitness},
};
use algebra::bytes::ToBytes;

#[derive(Debug)]
pub enum LedgerError {
    DuplicateSn,
    DuplicateMemo,
    InvalidCm,
    InvalidCmIndex,
}

impl std::fmt::Display for LedgerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let msg = match self {
            LedgerError::DuplicateSn => "duplicate sn pushed to ledger",
            LedgerError::DuplicateMemo => "duplicate memo pushed to ledger",
            LedgerError::InvalidCm => "invalid cm pushed to ledger",
            LedgerError::InvalidCmIndex => "invalid cm index during proving",

        };
        write!(f, "{}", msg)
    }
}

impl std::error::Error for LedgerError {
    #[inline]
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

#[derive(Derivative)]
#[derivative(
    Default(bound = "H: FixedLengthCRH"),
    Clone(bound = "H: FixedLengthCRH"),
    PartialEq(bound = "H: FixedLengthCRH"),
    Eq(bound = "H: FixedLengthCRH"),
    Hash(bound = "H: FixedLengthCRH"),
    Debug(bound = "H: FixedLengthCRH, H::Output: fmt::Debug")
)]
pub struct Digest<H: FixedLengthCRH>(pub(crate) H::Output);

impl<H: FixedLengthCRH> ToBytes for Digest<H> {
    #[inline]
    fn write<W: Write>(&self, writer: W) -> IoResult<()> {
        self.0.write(writer)
    }
}

impl<H: FixedLengthCRH> LedgerDigest for Digest<H> {
    type Parameters = H::Parameters;
}

#[derive(Derivative)]
#[derivative(
    Default(bound = "H: FixedLengthCRH, L: ToBytes + Eq"),
    Clone(bound = "H: FixedLengthCRH, L: ToBytes + Eq"),
    Debug(bound = "H: FixedLengthCRH, L: ToBytes + Eq, H::Output: fmt::Debug")
)]
pub struct CommPath<H: FixedLengthCRH, L: ToBytes + Eq>(pub(crate) HashMembershipProof<H, L>);

#[derive(Default, Clone)]
pub struct SnPath;

#[derive(Default, Clone)]
pub struct MemoPath;

impl<H: FixedLengthCRH> LedgerWitness<Digest<H>> for SnPath {
    fn dummy_witness() -> Self {
        SnPath
    }
}

impl<H: FixedLengthCRH> LedgerWitness<Digest<H>> for MemoPath {
    fn dummy_witness() -> Self {
        MemoPath
    }
}

impl<H: FixedLengthCRH, L: ToBytes + Eq> LedgerWitness<Digest<H>> for CommPath<H, L> {
    fn dummy_witness() -> Self {
        CommPath(HashMembershipProof::default())
    }
}

pub struct IdealLedger<T: Transaction, H: FixedLengthCRH>
where
    T::Commitment: ToBytes,
{
    crh_params:     Rc<H::Parameters>,
    transactions:   Vec<T>,
    cm_merkle_tree: MerkleHashTree<H, T::Commitment>,
    cur_cm_index:   usize,
    cur_sn_index:   usize,
    cur_memo_index: usize,
    comm_to_index:  HashMap<T::Commitment, usize>,
    sn_to_index:    HashMap<T::SerialNumber, usize>,
    memo_to_index:  HashMap<T::Memorandum, usize>,
    current_digest: Option<Digest<H>>,
    past_digests:   HashSet<Digest<H>>,
    genesis_cm:     T::Commitment,
    genesis_sn:     T::SerialNumber,
    genesis_memo:   T::Memorandum,
}

impl<T: Transaction, H: FixedLengthCRH> Ledger for IdealLedger<T, H>
where
    T: Eq,
    T::Commitment: ToBytes + Clone,
    T::SerialNumber: ToBytes + Clone,
    T::Memorandum: Hash + Clone,
{
    type Parameters = H::Parameters;
    type LedgerStateDigest = Digest<H>;
    type Commitment = T::Commitment;
    type CommWitness = CommPath<H, T::Commitment>;

    type SerialNumber = T::SerialNumber;
    type SnWitness = SnPath;

    type Memo = T::Memorandum;
    type MemoWitness = MemoPath;
    type Transaction = T;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Parameters, Error> {
        H::setup(rng)
    }

    fn new(
        parameters: Self::Parameters,
        genesis_cm: Self::Commitment,
        genesis_sn: Self::SerialNumber,
        genesis_memo: Self::Memo,
    ) -> Self {
        let params = Rc::new(parameters);
        let cm_merkle_tree = MerkleHashTree::new(params.clone(), &[genesis_cm.clone()]).unwrap();

        let mut cur_cm_index = 0;
        let mut comm_to_index = HashMap::new();
        comm_to_index.insert(genesis_cm.clone(), cur_cm_index);
        cur_cm_index += 1;

        let root = Digest(cm_merkle_tree.root());
        let mut past_digests = HashSet::new();
        past_digests.insert(root.clone());

        IdealLedger {
            crh_params: params,
            transactions: Vec::new(),
            cm_merkle_tree,
            cur_cm_index,
            cur_sn_index: 0,
            cur_memo_index: 0,

            comm_to_index,
            sn_to_index: HashMap::new(),
            memo_to_index: HashMap::new(),
            current_digest: Some(root),
            past_digests,
            genesis_cm,
            genesis_sn,
            genesis_memo,
        }
    }

    fn len(&self) -> usize {
        self.transactions.len()
    }

    fn parameters(&self) -> &Self::Parameters {
        &self.crh_params
    }

    fn push(&mut self, transaction: Self::Transaction) -> Result<(), Error> {
        let push_time = timer_start!(|| "IdealLedger::PushTx");

        let mut cur_sn_index = self.cur_sn_index;
        for sn in transaction.old_serial_numbers() {
            if sn != &self.genesis_sn {
                if self.sn_to_index.contains_key(sn) {
                    Err(LedgerError::DuplicateSn)?;
                }
                self.sn_to_index.insert(sn.clone(), cur_sn_index);
                cur_sn_index += 1;
            }
        }
        self.cur_sn_index = cur_sn_index;

        let mut cur_cm_index = self.cur_cm_index;
        for cm in transaction.new_commitments() {
            if cm == &self.genesis_cm || self.comm_to_index.contains_key(cm) {
                Err(LedgerError::InvalidCm)?;
            }
            self.comm_to_index.insert(cm.clone(), cur_cm_index);
            cur_cm_index += 1;
        }
        self.cur_cm_index = cur_cm_index;

        if transaction.memorandum() != &self.genesis_memo {
            if self.memo_to_index.contains_key(transaction.memorandum()) {
                Err(LedgerError::DuplicateMemo)?;
            } else {
                self.memo_to_index
                    .insert(transaction.memorandum().clone(), self.cur_memo_index);
                self.cur_memo_index += 1;
            }
        }

        // Rebuild the tree.
        let mut cm_and_indices = self.comm_to_index.iter().collect::<Vec<_>>();
        cm_and_indices.sort_by(|&(_, i), &(_, j)| i.cmp(j));
        let commitments = cm_and_indices
            .into_iter()
            .map(|(cm, _)| cm)
            .cloned()
            .collect::<Vec<_>>();
        assert!(commitments[0] == self.genesis_cm);
        self.cm_merkle_tree = MerkleHashTree::new(self.crh_params.clone(), &commitments)?;

        let new_digest = Digest(self.cm_merkle_tree.root());
        self.past_digests.insert(new_digest.clone());
        self.current_digest = Some(new_digest);

        self.transactions.push(transaction);

        timer_end!(push_time);
        Ok(())
    }

    fn digest(&self) -> Option<Self::LedgerStateDigest> {
        self.current_digest.clone()
    }

    fn validate_digest(&self, digest: &Self::LedgerStateDigest) -> bool {
        self.past_digests.contains(digest)
    }

    fn contains_cm(&self, cm: &Self::Commitment) -> bool {
        self.comm_to_index.contains_key(cm)
    }

    fn contains_sn(&self, sn: &Self::SerialNumber) -> bool {
        self.sn_to_index.contains_key(sn) && sn != &self.genesis_sn
    }

    fn contains_memo(&self, memo: &Self::Memo) -> bool {
        self.memo_to_index.contains_key(memo)
    }

    fn prove_cm(&self, cm: &Self::Commitment) -> Result<Self::CommWitness, Error> {
        let witness_time = timer_start!(|| "Generate membership witness");

        let cm_index = self
            .comm_to_index
            .get(cm)
            .ok_or(LedgerError::InvalidCmIndex)?;

        let result = CommPath(self.cm_merkle_tree.generate_proof(*cm_index, cm)?);

        timer_end!(witness_time);
        Ok(result)
    }

    fn prove_sn(&self, _sn: &Self::SerialNumber) -> Result<Self::SnWitness, Error> {
        Ok(SnPath)
    }

    fn prove_memo(&self, _memo: &Self::Memo) -> Result<Self::MemoWitness, Error> {
        Ok(MemoPath)
    }

    fn verify_cm(
        params: &Self::Parameters,
        digest: &Self::LedgerStateDigest,
        cm: &Self::Commitment,
        witness: &Self::CommWitness,
    ) -> bool {
        witness.0.verify(params, &digest.0, cm).unwrap()
    }

    fn verify_sn(
        _params: &Self::Parameters,
        _digest: &Self::LedgerStateDigest,
        _sn: &Self::SerialNumber,
        _witness: &Self::SnWitness,
    ) -> bool {
        true
    }

    fn verify_memo(
        _params: &Self::Parameters,
        _digest: &Self::LedgerStateDigest,
        _memo: &Self::Memo,
        _witness: &Self::MemoWitness,
    ) -> bool {
        true
    }
}
