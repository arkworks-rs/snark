use crate::{
    dpc::delegable_dpc::{DelegableDPCComponents, Transaction},
    ledger::*,
};
use crypto_primitives::{CommitmentScheme, SignatureScheme, NIZK};

#[derive(Derivative)]
#[derivative(
    Clone(bound = "C: DelegableDPCComponents"),
    PartialEq(bound = "C: DelegableDPCComponents"),
    Eq(bound = "C: DelegableDPCComponents")
)]
pub struct DPCTransaction<C: DelegableDPCComponents> {
    old_serial_numbers: Vec<<C::S as SignatureScheme>::PublicKey>,
    new_commitments: Vec<<C::RecC as CommitmentScheme>::Output>,
    memorandum: [u8; 32],
    pub stuff: DPCStuff<C>,
}

#[derive(Derivative)]
#[derivative(
    Clone(bound = "C: DelegableDPCComponents"),
    PartialEq(bound = "C: DelegableDPCComponents"),
    Eq(bound = "C: DelegableDPCComponents")
)]
pub struct DPCStuff<C: DelegableDPCComponents> {
    pub digest: merkle_tree::Digest<C::MerkleTreeConfig>,

    #[derivative(PartialEq = "ignore")]
    pub core_proof: <C::MainNIZK as NIZK>::Proof,

    #[derivative(PartialEq = "ignore")]
    pub predicate_proof: <C::ProofCheckNIZK as NIZK>::Proof,

    #[derivative(PartialEq = "ignore")]
    pub predicate_comm: <C::PredVkComm as CommitmentScheme>::Output,

    #[derivative(PartialEq = "ignore")]
    pub local_data_comm: <C::LocalDataComm as CommitmentScheme>::Output,

    #[derivative(PartialEq = "ignore")]
    pub signatures: Vec<<C::S as SignatureScheme>::Signature>,
}

impl<C: DelegableDPCComponents> DPCTransaction<C> {
    pub fn new(
        old_serial_numbers: Vec<<Self as Transaction>::SerialNumber>,
        new_commitments: Vec<<Self as Transaction>::Commitment>,
        memorandum: <Self as Transaction>::Memorandum,
        digest: merkle_tree::Digest<C::MerkleTreeConfig>,
        core_proof: <C::MainNIZK as NIZK>::Proof,
        predicate_proof: <C::ProofCheckNIZK as NIZK>::Proof,
        predicate_comm: <C::PredVkComm as CommitmentScheme>::Output,
        local_data_comm: <C::LocalDataComm as CommitmentScheme>::Output,
        signatures: Vec<<C::S as SignatureScheme>::Signature>,
    ) -> Self {
        let stuff = DPCStuff {
            digest,

            core_proof,
            predicate_proof,

            predicate_comm,
            local_data_comm,

            signatures,
        };
        DPCTransaction {
            old_serial_numbers,
            new_commitments,
            memorandum,
            stuff,
        }
    }
}

impl<C: DelegableDPCComponents> Transaction for DPCTransaction<C> {
    type Stuff = DPCStuff<C>;
    type SerialNumber = <C::S as SignatureScheme>::PublicKey;
    type Commitment = <C::RecC as CommitmentScheme>::Output;
    type Memorandum = [u8; 32];

    fn old_serial_numbers(&self) -> &[Self::SerialNumber] {
        self.old_serial_numbers.as_slice()
    }

    fn new_commitments(&self) -> &[Self::Commitment] {
        self.new_commitments.as_slice()
    }

    fn memorandum(&self) -> &Self::Memorandum {
        &self.memorandum
    }

    fn stuff(&self) -> &Self::Stuff {
        &self.stuff
    }
}
