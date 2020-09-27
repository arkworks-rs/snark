use algebra::{
    bls12_377, cp6_782::CP6_782, ed_on_bls12_377::EdwardsProjective as EdwardsBls,
    ed_on_cp6_782::EdwardsProjective as EdwardsCP6, Bls12_377,
};

use crypto_primitives::{
    commitment::{blake2s, injective_map::PedersenCommCompressor},
    crh::{
        injective_map::{PedersenCRHCompressor, TECompressor},
        pedersen,
    },
    merkle_tree,
    nizk::Groth16,
    prf::blake2s::Blake2s,
};

use crypto_primitives::{
    commitment::injective_map,
    crh::injective_map::constraints::{PedersenCRHCompressorGadget, TECompressorGadget},
    nizk::groth16::constraints::Groth16VerifierGadget,
    prf::blake2s::constraints::Blake2sGadget,
};
use r1cs_std::{
    bls12_377::PairingVar, ed_on_bls12_377::EdwardsVar as EdwardsBlsVar,
    ed_on_cp6_782::EdwardsVar as EdwardsCP6Var,
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
impl pedersen::Window for SnNonceWindow {
    const WINDOW_SIZE: usize = SN_NONCE_SIZE_BITS / 8;
    const NUM_WINDOWS: usize = 8;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct PredVkHashWindow;

impl pedersen::Window for PredVkHashWindow {
    const WINDOW_SIZE: usize = 248;
    const NUM_WINDOWS: usize = 38;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct LocalDataWindow;

impl pedersen::Window for LocalDataWindow {
    const WINDOW_SIZE: usize = 248;
    const NUM_WINDOWS: usize = 30;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct TwoToOneWindow;
// `WINDOW_SIZE * NUM_WINDOWS` = 2 * 256 bits
impl pedersen::Window for TwoToOneWindow {
    const WINDOW_SIZE: usize = 128;
    const NUM_WINDOWS: usize = 4;
}

pub struct CommitmentMerkleTreeConfig;
impl merkle_tree::Config for CommitmentMerkleTreeConfig {
    const HEIGHT: usize = 32;
    type H = MerkleTreeCRH;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct RecordWindow;
impl pedersen::Window for RecordWindow {
    const WINDOW_SIZE: usize = 225;
    const NUM_WINDOWS: usize = 8;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct AddressWindow;
impl pedersen::Window for AddressWindow {
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
    type MerkleTreeHGadget = MerkleTreeCRHGadget;

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
pub type ProofCheckPairing = CP6_782;
pub type CoreCheckF = bls12_377::Fr;
pub type ProofCheckF = bls12_377::Fq;

pub type AddressComm = PedersenCommCompressor<EdwardsBls, EdwardsCompressor, AddressWindow>;
pub type RecordComm = PedersenCommCompressor<EdwardsBls, EdwardsCompressor, RecordWindow>;
pub type PredicateComm = blake2s::Commitment;
pub type LocalDataComm = PedersenCommCompressor<EdwardsBls, EdwardsCompressor, LocalDataWindow>;

pub type MerkleTreeCRH = PedersenCRHCompressor<EdwardsBls, EdwardsCompressor, TwoToOneWindow>;
pub type SnNonceCRH = PedersenCRHCompressor<EdwardsBls, EdwardsCompressor, SnNonceWindow>;
pub type PredVkCRH = PedersenCRHCompressor<EdwardsCP6, EdwardsCompressor, PredVkHashWindow>;

pub type Predicate = DPCPredicate<Components>;
pub type CoreCheckNIZK =
    Groth16<CoreCheckPairing, CoreChecksCircuit<Components>, CoreChecksVerifierInput<Components>>;
pub type ProofCheckNIZK =
    Groth16<ProofCheckPairing, ProofCheckCircuit<Components>, ProofCheckVerifierInput<Components>>;
pub type PredicateNIZK<C> =
    Groth16<CoreCheckPairing, EmptyPredicateCircuit<C>, PredicateLocalData<C>>;
pub type PRF = Blake2s;

// Gadgets
pub type EdwardsCompressorGadget = TECompressorGadget;

pub type RecordCommGadget = injective_map::constraints::CommitmentCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    RecordWindow,
    EdwardsBlsVar,
    EdwardsCompressorGadget,
>;
pub type AddressCommGadget = injective_map::constraints::CommitmentCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    AddressWindow,
    EdwardsBlsVar,
    EdwardsCompressorGadget,
>;
pub type PredicateCommGadget = blake2s::constraints::CommGadget;
pub type LocalDataCommGadget = injective_map::constraints::CommitmentCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    LocalDataWindow,
    EdwardsBlsVar,
    EdwardsCompressorGadget,
>;

pub type SnNonceCRHGadget = PedersenCRHCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    SnNonceWindow,
    EdwardsBlsVar,
    EdwardsCompressorGadget,
>;
pub type MerkleTreeCRHGadget = PedersenCRHCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    TwoToOneWindow,
    EdwardsBlsVar,
    EdwardsCompressorGadget,
>;
pub type PredVkCRHGadget = PedersenCRHCompressorGadget<
    EdwardsCP6,
    EdwardsCompressor,
    PredVkHashWindow,
    EdwardsCP6Var,
    EdwardsCompressorGadget,
>;

pub type PRFGadget = Blake2sGadget;
pub type PredicateNIZKGadget = Groth16VerifierGadget<CoreCheckPairing, PairingVar>;
//

pub type MerkleTreeIdealLedger = IdealLedger<Tx, CommitmentMerkleTreeConfig>;
pub type Tx = DPCTransaction<Components>;

pub type InstantiatedDPC = DPC<Components>;
pub type LocalData = DPCLocalData<Components>;
