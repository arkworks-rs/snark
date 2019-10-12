use algebra::curves::{
    bls12_377::Bls12_377, edwards_bls12::EdwardsProjective as EdwardsBls,
    edwards_sw6::EdwardsProjective as EdwardsSW, sw6::SW6,
};

use algebra::fields::bls12_377::{fr::Fr as Bls12_377Fr, fq::Fq as Bls12_377Fq};

use crypto_primitives::{
    commitment::{blake2s::Blake2sCommitment, injective_map::PedersenCommCompressor},
    crh::{
        injective_map::{PedersenCRHCompressor, TECompressor},
        pedersen::PedersenWindow,
    },
    nizk::Gm17,
    merkle_tree::MerkleTreeConfig,
    prf::blake2s::Blake2s,
};

use crypto_primitives::{
    commitment::{
        blake2s::constraints::Blake2sCommitmentGadget, 
        injective_map::constraints::PedersenCommitmentCompressorGadget,
    },
    crh::injective_map::constraints::{PedersenCRHCompressorGadget, TECompressorGadget},
    prf::blake2s::constraints::Blake2sGadget,
    nizk::gm17::constraints::Gm17VerifierGadget,
};
use r1cs_std::{
    groups::curves::twisted_edwards::{
        edwards_bls12::EdwardsBlsGadget, edwards_sw6::EdwardsSWGadget,
    },
    pairing::bls12_377::PairingGadget,
};

use crate::dpc::plain_dpc::{
    core_checks_circuit::*, predicate::DPCPredicate, predicate_circuit::*, proof_check_circuit::*,
    transaction::DPCTransaction, LocalData as DPCLocalData, PlainDPCComponents, DPC,
};

use crate::ledger::*;

pub const NUM_INPUT_RECORDS: usize = 2;
pub const NUM_OUTPUT_RECORDS: usize = 2;

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct SnNonceWindow;

// `WINDOW_SIZE * NUM_WINDOWS` = 2 * 256 + 8 + 256 bits
const SN_NONCE_SIZE_BITS: usize = NUM_INPUT_RECORDS * 2 * 256 + 8 + 256;
impl PedersenWindow for SnNonceWindow {
    const WINDOW_SIZE: usize = SN_NONCE_SIZE_BITS / 8;
    const NUM_WINDOWS: usize = 8;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct PredVkHashWindow;

impl PedersenWindow for PredVkHashWindow {
    const WINDOW_SIZE: usize = 248;
    const NUM_WINDOWS: usize = 38;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct LocalDataWindow;

impl PedersenWindow for LocalDataWindow {
    const WINDOW_SIZE: usize = 248;
    const NUM_WINDOWS: usize = 30;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct TwoToOneWindow;
// `WINDOW_SIZE * NUM_WINDOWS` = 2 * 256 bits
impl PedersenWindow for TwoToOneWindow {
    const WINDOW_SIZE: usize = 128;
    const NUM_WINDOWS: usize = 4;
}

pub struct CommitmentMerkleTreeConfig;
impl MerkleTreeConfig for CommitmentMerkleTreeConfig {
    const HEIGHT: usize = 32;
    type H = MerkleTreeCRH;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct RecordWindow;
impl PedersenWindow for RecordWindow {
    const WINDOW_SIZE: usize = 225;
    const NUM_WINDOWS: usize = 8;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct AddressWindow;
impl PedersenWindow for AddressWindow {
    const WINDOW_SIZE: usize = 128;
    const NUM_WINDOWS: usize = 4;
}

pub struct Components;

impl PlainDPCComponents for Components {
    const NUM_INPUT_RECORDS: usize = NUM_INPUT_RECORDS;
    const NUM_OUTPUT_RECORDS: usize = NUM_OUTPUT_RECORDS;

    type CoreCheckF = CoreCheckF;
    type ProofCheckF = ProofCheckF;

    type MerkleTreeConfig = CommitmentMerkleTreeConfig;
    type MerkleTree_HGadget = MerkleTreeCRHGadget;

    type AddrC = AddressComm;
    type RecC = RecordComm;

    type AddrCGadget = AddressCommGadget;
    type RecCGadget = RecordCommGadget;

    type SnNonceH = SnNonceCRH;
    type SnNonceHGadget = SnNonceCRHGadget;
    type MainNIZK = CoreCheckNIZK;
    type ProofCheckNIZK = ProofCheckNIZK;
    type P = PRF;
    type PGadget = PRFGadget;

    type PredicateNIZK = PredicateNIZK<Self>;
    type PredicateNIZKGadget = PredicateNIZKGadget;

    type PredVkH = PredVkCRH;
    type PredVkHGadget = PredVkCRHGadget;
    type PredVkComm = PredicateComm;
    type PredVkCommGadget = PredicateCommGadget;
    type LocalDataComm = LocalDataComm;
    type LocalDataCommGadget = LocalDataCommGadget;
}

// Native primitives
pub type EdwardsCompressor = TECompressor;
pub type CoreCheckPairing = Bls12_377;
pub type ProofCheckPairing = SW6;
pub type CoreCheckF = Bls12_377Fr;
pub type ProofCheckF = Bls12_377Fq;

pub type AddressComm = PedersenCommCompressor<EdwardsBls, EdwardsCompressor, AddressWindow>;
pub type RecordComm = PedersenCommCompressor<EdwardsBls, EdwardsCompressor, RecordWindow>;
pub type PredicateComm = Blake2sCommitment;
pub type LocalDataComm = PedersenCommCompressor<EdwardsBls, EdwardsCompressor, LocalDataWindow>;

pub type MerkleTreeCRH = PedersenCRHCompressor<EdwardsBls, EdwardsCompressor, TwoToOneWindow>;
pub type SnNonceCRH = PedersenCRHCompressor<EdwardsBls, EdwardsCompressor, SnNonceWindow>;
pub type PredVkCRH = PedersenCRHCompressor<EdwardsSW, EdwardsCompressor, PredVkHashWindow>;

pub type Predicate = DPCPredicate<Components>;
pub type CoreCheckNIZK =
    Gm17<CoreCheckPairing, CoreChecksCircuit<Components>, CoreChecksVerifierInput<Components>>;
pub type ProofCheckNIZK =
    Gm17<ProofCheckPairing, ProofCheckCircuit<Components>, ProofCheckVerifierInput<Components>>;
pub type PredicateNIZK<C> = Gm17<CoreCheckPairing, EmptyPredicateCircuit<C>, PredicateLocalData<C>>;
pub type PRF = Blake2s;

// Gadgets
pub type EdwardsCompressorGadget = TECompressorGadget;

pub type RecordCommGadget = PedersenCommitmentCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreCheckF,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;
pub type AddressCommGadget = PedersenCommitmentCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreCheckF,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;
pub type PredicateCommGadget = Blake2sCommitmentGadget;
pub type LocalDataCommGadget = PedersenCommitmentCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreCheckF,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;

pub type SnNonceCRHGadget = PedersenCRHCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreCheckF,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;
pub type MerkleTreeCRHGadget = PedersenCRHCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreCheckF,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;
pub type PredVkCRHGadget = PedersenCRHCompressorGadget<
    EdwardsSW,
    EdwardsCompressor,
    ProofCheckF,
    EdwardsSWGadget,
    EdwardsCompressorGadget,
>;

pub type PRFGadget = Blake2sGadget;
pub type PredicateNIZKGadget = Gm17VerifierGadget<CoreCheckPairing, ProofCheckF, PairingGadget>;
//

pub type MerkleTreeIdealLedger = IdealLedger<Tx, CommitmentMerkleTreeConfig>;
pub type Tx = DPCTransaction<Components>;

pub type InstantiatedDPC = DPC<Components>;
pub type LocalData = DPCLocalData<Components>;
