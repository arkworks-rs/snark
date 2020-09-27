use algebra::curves::bls12_377::Bls12_377;
use algebra::curves::cp6_782::CP6_782;
use algebra::curves::ed_on_bls12_377::EdwardsProjective as EdwardsBls;
use algebra::curves::ed_on_cp6_782::EdwardsProjective as E;

use algebra::{to_bytes, ToBytes};
use rand::thread_rng;

use crate::crypto_primitives::{
    commitment::{blake2s::Blake2sCommitment, injective_map::PedersenCommCompressor},
    crh::{
        injective_map::{PedersenCRHCompressor, TECompressor},
        pedersen::PedersenWindow,
    },
    nizk::gm17_unprep::Gm17 as Gm17UnPrep,
    nizk::Gm17,
    prf::blake2s::Blake2s,
    signature::schnorr::SchnorrSignature,
};
use crate::crypto_primitives::{CommitmentScheme, FixedLengthCRH};
use blake2::Blake2s as Blake2sHash;

use r1cs_core::ConstraintSystem;

use crate::constraints::commitment::{
    blake2s::Blake2sCommitmentGadget, injective_map::PedersenCommitmentCompressorGadget,
};
use crate::constraints::crh::injective_map::{PedersenCRHCompressorGadget, TECompressorGadget};
use crate::constraints::delegable_dpc::execute_core_checks_gadget;
use crate::constraints::delegable_dpc::execute_proof_check_gadget;
use crate::constraints::merkle_tree::IdealLedgerGadget;
use crate::constraints::prf::blake2s::Blake2sGadget;
use crate::constraints::verifier::gm17::Gm17VerifierGadget;
use r1cs_std::groups::curves::twisted_edwards::ed_on_bls12_377::EdwardsGadget as EdwardsBlsGadget;
use r1cs_std::groups::curves::twisted_edwards::ed_on_cp6_782::EdwardsGadget as EdwardsCP6Gadget;
use r1cs_std::pairing::bls12_377::PairingGadget;
use r1cs_std::test_constraint_system::TestConstraintSystem;

use crate::constraints::signature::schnorr::SchnorrRandomizePkGadget;

use crate::dpc::{DPCScheme, Record};

use crate::dpc::delegable_dpc::core_checks_circuit::*;
use crate::dpc::delegable_dpc::predicate::DPCPredicate;
use crate::dpc::delegable_dpc::predicate::PrivatePredInput;
use crate::dpc::delegable_dpc::predicate_circuit::*;
use crate::dpc::delegable_dpc::proof_check_circuit::*;
use crate::dpc::delegable_dpc::transaction::DPCTransaction;
use crate::dpc::delegable_dpc::DelegableDPCComponents;
use crate::dpc::delegable_dpc::ExecuteContext;
use crate::dpc::delegable_dpc::DPC;

use crate::ledger::{CommPath, Digest, IdealLedger, Ledger};

const NUM_INPUT_RECORDS: usize = 2;
const NUM_OUTPUT_RECORDS: usize = 2;

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct SnNonceWindow;

// `WINDOW_SIZE * NUM_WINDOWS` = NUM_INPUT_RECORDS * 64 + 1 + 32 = 225 bytes
const SN_NONCE_SIZE_BITS: usize = NUM_INPUT_RECORDS * 2 * 512 + 8 + 256;
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
    const NUM_WINDOWS: usize = 36;
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct TwoToOneWindow;
// `WINDOW_SIZE * NUM_WINDOWS` = 2 * 256 bits
impl PedersenWindow for TwoToOneWindow {
    const WINDOW_SIZE: usize = 128;
    const NUM_WINDOWS: usize = 4;
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
    const WINDOW_SIZE: usize = 192;
    const NUM_WINDOWS: usize = 8;
}

struct Components;

impl DelegableDPCComponents for Components {
    const NUM_INPUT_RECORDS: usize = NUM_INPUT_RECORDS;
    const NUM_OUTPUT_RECORDS: usize = NUM_OUTPUT_RECORDS;
    type E = CoreEngine;
    type ProofCheckE = ProofCheckEngine;

    type AddrC = AddressComm;
    type RecC = RecordComm;

    type AddrCGadget = AddressCommGadget;
    type RecCGadget = RecordCommGadget;

    type D = MerkleTreeDigest;
    type SnNonceH = SnNonceCRH;
    type SnNonceHGadget = SnNonceCRHGadget;
    type MainNIZK = CoreCheckNIZK;
    type ProofCheckNIZK = ProofCheckNIZK;

    type LCW = MerkleTreeWitness;
    type LCWGadget = MerkleTreeWitnessGadget;
    type P = PRF;
    type PGadget = PRFGadget;

    type PredicateNIZK = PredicateNIZK<Self>;
    type PredicateNIZKGadget = PredicateNIZKGadget;

    type S = AuthSignature;
    type SGadget = AuthSignatureGadget;

    type PredVkH = PredVkCRH;
    type PredVkHGadget = PredVkCRHGadget;
    type PredVkComm = PredicateComm;
    type PredVkCommGadget = PredicateCommGadget;
    type LocalDataComm = LocalDataComm;
    type LocalDataCommGadget = LocalDataCommGadget;
}

// Native primitives
type EdwardsCompressor = TECompressor;
type CoreEngine = Bls12_377;
type ProofCheckEngine = CP6_782;

type AddressComm = PedersenCommCompressor<EdwardsBls, EdwardsCompressor, AddressWindow>;
type RecordComm = PedersenCommCompressor<EdwardsBls, EdwardsCompressor, RecordWindow>;
type PredicateComm = Blake2sCommitment;
type LocalDataComm = PedersenCommCompressor<EdwardsBls, EdwardsCompressor, LocalDataWindow>;

type AuthSignature = SchnorrSignature<EdwardsBls, Blake2sHash>;
type MerkleTreeCRH = PedersenCRHCompressor<EdwardsBls, EdwardsCompressor, TwoToOneWindow>;
type SnNonceCRH = PedersenCRHCompressor<EdwardsBls, EdwardsCompressor, SnNonceWindow>;
type PredVkCRH = PedersenCRHCompressor<EdwardsCP6, EdwardsCompressor, PredVkHashWindow>;

type Predicate = DPCPredicate<Components>;
type CoreCheckNIZK =
    Gm17<CoreEngine, CoreChecksCircuit<Components>, CoreChecksVerifierInput<Components>>;
type ProofCheckNIZK =
    Gm17<ProofCheckEngine, ProofCheckCircuit<Components>, ProofCheckVerifierInput<Components>>;
type PredicateNIZK<C> = Gm17UnPrep<CoreEngine, EmptyPredicateCircuit<C>, PredicateLocalData<C>>;
type PRF = Blake2s;

type MerkleTreeDigest = Digest<MerkleTreeCRH>;
type MerkleTreeWitness = CommPath<MerkleTreeCRH, <RecordComm as CommitmentScheme>::Output>;
//

// Gadgets
type EdwardsCompressorGadget = TECompressorGadget;

type RecordCommGadget = PedersenCommitmentCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreEngine,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;
type AddressCommGadget = PedersenCommitmentCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreEngine,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;
type PredicateCommGadget = Blake2sCommitmentGadget;
type LocalDataCommGadget = PedersenCommitmentCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreEngine,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;

type SnNonceCRHGadget = PedersenCRHCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreEngine,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;
type MerkleTreeCRHGadget = PedersenCRHCompressorGadget<
    EdwardsBls,
    EdwardsCompressor,
    CoreEngine,
    EdwardsBlsGadget,
    EdwardsCompressorGadget,
>;
type PredVkCRHGadget = PedersenCRHCompressorGadget<
    EdwardsCP6,
    EdwardsCompressor,
    ProofCheckEngine,
    EdwardsCP6Gadget,
    EdwardsCompressorGadget,
>;

type AuthSignatureGadget = SchnorrRandomizePkGadget<EdwardsBls, Bls12_377, EdwardsBlsGadget>;
type MerkleTreeWitnessGadget =
    IdealLedgerGadget<RecordComm, MerkleTreeCRH, MerkleTreeCRHGadget, RecordCommGadget>;
type PRFGadget = Blake2sGadget;
type PredicateNIZKGadget = Gm17VerifierGadget<CoreEngine, ProofCheckEngine, PairingGadget>;
//

type MerkleTreeIdealLedger = IdealLedger<Tx, MerkleTreeCRH>;
type Tx = DPCTransaction<Components>;

type InstantiatedDPC = DPC<Components>;

#[test]
fn core_checks_gadget() {
    let mut rng = thread_rng();
    let ledger_parameters = MerkleTreeIdealLedger::setup(&mut rng).expect("Ledger setup failed");
    let comm_crh_sig_pp = InstantiatedDPC::generate_comm_crh_sig_parameters(&mut rng).unwrap();
    let pred_nizk_pp =
        InstantiatedDPC::generate_pred_nizk_parameters(&comm_crh_sig_pp, &mut rng).unwrap();
    let pred_nizk_vk_bytes = to_bytes![PredVkCRH::evaluate(
        &comm_crh_sig_pp.pred_vk_crh_pp,
        &to_bytes![pred_nizk_pp.vk].unwrap()
    )
    .unwrap()]
    .unwrap();

    // let dummy_pred_input_bytes = to_bytes![PredicateHashInput::<Components>::default()].unwrap();
    // println!("Predicate input length: {:?}", dummy_pred_input_bytes.len());
    // let _result = LocalDataComm::evaluate(
    //     &comm_crh_sig_pp.pred_input_crh_pp,
    //     &dummy_pred_input_bytes,
    // ).unwrap();

    // Create genesis record
    let genesis_metadata = [1u8; 32];
    let genesis_address =
        InstantiatedDPC::create_address_helper(&comm_crh_sig_pp, &genesis_metadata, &mut rng)
            .unwrap();
    let genesis_sn_nonce =
        SnNonceCRH::evaluate(&comm_crh_sig_pp.sn_nonce_crh_pp, &[0u8; 1]).unwrap();
    let genesis_record = InstantiatedDPC::generate_record(
        &comm_crh_sig_pp,
        &genesis_sn_nonce,
        &genesis_address.public_key,
        true,
        &[0u8; 32],
        &Predicate::new(pred_nizk_vk_bytes.clone()),
        &Predicate::new(pred_nizk_vk_bytes.clone()),
        &mut rng,
    )
    .unwrap();

    // Create ledger
    let (genesis_sn, _) = InstantiatedDPC::generate_sn(
        &comm_crh_sig_pp,
        &genesis_record,
        &genesis_address.secret_key,
    )
    .unwrap();
    println!(
        "genesis_record.len(): {:?}",
        to_bytes![genesis_record.commitment()].unwrap().len()
    );
    let genesis_memo = [0u8; 32];
    let ledger = MerkleTreeIdealLedger::new(
        ledger_parameters,
        genesis_record.commitment(),
        genesis_sn.clone(),
        genesis_memo,
    );

    // Create address 1
    let metadata1 = [1u8; 32];
    let address1 =
        InstantiatedDPC::create_address_helper(&comm_crh_sig_pp, &metadata1, &mut rng).unwrap();

    let test_payload = [1u8; 32];
    let test_predicate = Predicate::new(pred_nizk_vk_bytes.clone());

    println!("Execution 0");
    let old_records1 = vec![genesis_record.clone(); NUM_INPUT_RECORDS];
    let old_asks1 = vec![genesis_address.secret_key.clone(); NUM_INPUT_RECORDS];

    let new_apks1 = vec![address1.public_key.clone(); NUM_OUTPUT_RECORDS];
    let new_payloads1 = vec![test_payload.clone(); NUM_OUTPUT_RECORDS];
    let new_birth_predicates1 = vec![test_predicate.clone(); NUM_OUTPUT_RECORDS];
    let new_death_predicates1 = vec![test_predicate.clone(); NUM_OUTPUT_RECORDS];
    let auxiliary = [0u8; 32];
    let memo = [0u8; 32];

    let context = InstantiatedDPC::execute_helper(
        &comm_crh_sig_pp,
        &old_records1,
        &old_asks1,
        &new_apks1,
        &[false, false],
        &new_payloads1,
        &new_birth_predicates1,
        &new_death_predicates1,
        &memo,
        &auxiliary,
        &ledger,
        &mut rng,
    )
    .unwrap();

    let ExecuteContext {
        comm_crh_sig_pp: _comm_crh_sig_pp,
        ledger_digest,

        old_records,
        old_witnesses,
        old_address_secret_keys,
        old_serial_numbers,
        old_randomizers: _,

        new_records,
        new_sn_nonce_randomness,
        new_commitments,

        predicate_comm,
        predicate_rand,

        local_data_comm,
        local_data_rand,
    } = context;

    let mut core_cs = TestConstraintSystem::<Bls12_377>::new();

    execute_core_checks_gadget::<_, _>(
        &mut core_r1cs_core::ns!(cs, || "Core checks"),
        &comm_crh_sig_pp,
        ledger.parameters(),
        &ledger_digest,
        &old_records,
        &old_witnesses,
        &old_address_secret_keys,
        &old_serial_numbers,
        &new_records,
        &new_sn_nonce_randomness,
        &new_commitments,
        &predicate_comm,
        &predicate_rand,
        &local_data_comm,
        &local_data_rand,
        &[0u8; 32],
        &auxiliary,
    )
    .unwrap();

    if !core_cs.is_satisfied() {
        println!("=========================================================");
        println!("Unsatisfied constraints:");
        println!("{}", core_cs.which_is_unsatisfied().unwrap());
        println!("=========================================================");
    }
    println!("\n\n\n\nAll Core check constraints:");
    core_cs.print_named_objects();

    assert!(core_cs.is_satisfied());
}

#[test]
fn proof_check_gadget() {
    let mut rng = thread_rng();
    let ledger_parameters = MerkleTreeIdealLedger::setup(&mut rng).expect("Ledger setup failed");
    let comm_crh_sig_pp = InstantiatedDPC::generate_comm_crh_sig_parameters(&mut rng).unwrap();
    let pred_nizk_pp =
        InstantiatedDPC::generate_pred_nizk_parameters(&comm_crh_sig_pp, &mut rng).unwrap();
    let pred_nizk_vk_bytes = to_bytes![PredVkCRH::evaluate(
        &comm_crh_sig_pp.pred_vk_crh_pp,
        &to_bytes![pred_nizk_pp.vk].unwrap()
    )
    .unwrap()]
    .unwrap();

    // Create genesis record
    let genesis_metadata = [1u8; 32];
    let genesis_address =
        DPC::create_address_helper(&comm_crh_sig_pp, &genesis_metadata, &mut rng).unwrap();
    let genesis_sn_nonce =
        SnNonceCRH::evaluate(&comm_crh_sig_pp.sn_nonce_crh_pp, &[0u8; 1]).unwrap();
    let genesis_record = DPC::generate_record(
        &comm_crh_sig_pp,
        &genesis_sn_nonce,
        &genesis_address.public_key,
        true,
        &[0u8; 32],
        &Predicate::new(pred_nizk_vk_bytes.clone()),
        &Predicate::new(pred_nizk_vk_bytes.clone()),
        &mut rng,
    )
    .unwrap();

    // Create ledger
    let (genesis_sn, _) = InstantiatedDPC::generate_sn(
        &comm_crh_sig_pp,
        &genesis_record,
        &genesis_address.secret_key,
    )
    .unwrap();
    println!(
        "genesis_record.cm.len(): {:?}",
        to_bytes![genesis_record.commitment()].unwrap().len()
    );
    let genesis_memo = [0u8; 32];
    let ledger = MerkleTreeIdealLedger::new(
        ledger_parameters,
        genesis_record.commitment(),
        genesis_sn.clone(),
        genesis_memo,
    );

    // Create address 1
    let metadata1 = [1u8; 32];
    let address1 =
        InstantiatedDPC::create_address_helper(&comm_crh_sig_pp, &metadata1, &mut rng).unwrap();

    let test_payload = [1u8; 32];
    let test_predicate = Predicate::new(pred_nizk_vk_bytes.clone());

    println!("Execution 0");
    let old_records1 = vec![genesis_record.clone(); NUM_INPUT_RECORDS];
    let old_asks1 = vec![genesis_address.secret_key.clone(); NUM_INPUT_RECORDS];

    let new_apks1 = vec![address1.public_key.clone(); NUM_OUTPUT_RECORDS];
    let new_payloads1 = vec![test_payload.clone(); NUM_OUTPUT_RECORDS];
    let new_birth_predicates1 = vec![test_predicate.clone(); NUM_OUTPUT_RECORDS];
    let new_death_predicates1 = vec![test_predicate.clone(); NUM_OUTPUT_RECORDS];
    let auxiliary = [0u8; 32];
    let memo = [0u8; 32];

    let context = InstantiatedDPC::execute_helper(
        &comm_crh_sig_pp,
        &old_records1,
        &old_asks1,
        &new_apks1,
        &[false, false],
        &new_payloads1,
        &new_birth_predicates1,
        &new_death_predicates1,
        &memo,
        &auxiliary,
        &ledger,
        &mut rng,
    )
    .unwrap();

    let private_input: PrivatePredInput<Components> = PrivatePredInput {
        vk: pred_nizk_pp.vk.clone(),
        proof: pred_nizk_pp.proof.clone(),
    };
    let old_proof_and_vk = vec![private_input.clone(); NUM_INPUT_RECORDS];
    let new_proof_and_vk = vec![private_input.clone(); NUM_OUTPUT_RECORDS];

    let ExecuteContext {
        comm_crh_sig_pp: _comm_crh_sig_pp,
        ledger_digest: _ledger_digest,

        old_records: _,
        old_witnesses: _,
        old_address_secret_keys: _,
        old_serial_numbers: _,
        old_randomizers: _,

        new_records: _,
        new_sn_nonce_randomness: _,
        new_commitments: _,

        predicate_comm,
        predicate_rand,

        local_data_comm,
        local_data_rand: _,
    } = context;

    let mut pf_check_cs = TestConstraintSystem::<CP6_782>::new();

    execute_proof_check_gadget::<_, _>(
        &mut pf_check_r1cs_core::ns!(cs, || "Check predicate proofs"),
        &comm_crh_sig_pp,
        &old_proof_and_vk,
        &new_proof_and_vk,
        &predicate_comm,
        &predicate_rand,
        &local_data_comm,
    )
    .unwrap();
    if !pf_check_cs.is_satisfied() {
        println!("=========================================================");
        println!("Unsatisfied constraints:");
        println!("{}", pf_check_cs.which_is_unsatisfied().unwrap());
        println!("=========================================================");
    }
    println!("\n\n\n\nAll Proof check constraints:");
    pf_check_cs.print_named_objects();
    assert!(pf_check_cs.is_satisfied());
}

#[test]
fn execution() {
    let rng = &mut thread_rng();

    // DPC Setup
    let ledger_parameters = MerkleTreeIdealLedger::setup(rng).expect("Ledger setup failed");
    let parameters =
        <InstantiatedDPC as DPCScheme<MerkleTreeIdealLedger>>::setup(&ledger_parameters, rng)
            .expect("DPC setup failed");

    // Create genesis record
    let genesis_metadata = [1u8; 32];
    let genesis_address = <InstantiatedDPC as DPCScheme<MerkleTreeIdealLedger>>::create_address(
        &parameters,
        &genesis_metadata,
        rng,
    )
    .unwrap();

    let genesis_sn_nonce = SnNonceCRH::evaluate(&parameters.sn_nonce_crh_pp(), &[0u8; 1]).unwrap();
    let genesis_pred_vk_bytes = to_bytes![PredVkCRH::evaluate(
        &parameters.comm_crh_sig_pp.pred_vk_crh_pp,
        &to_bytes![parameters.pred_nizk_pp.vk].unwrap()
    )
    .unwrap()]
    .unwrap();
    let genesis_record = InstantiatedDPC::generate_record(
        &parameters.comm_crh_sig_pp,
        &genesis_sn_nonce,
        &genesis_address.public_key,
        true,
        &[0u8; 32],
        &Predicate::new(genesis_pred_vk_bytes.clone()),
        &Predicate::new(genesis_pred_vk_bytes.clone()),
        rng,
    )
    .unwrap();

    // Create ledger
    let (genesis_sn, _) = InstantiatedDPC::generate_sn(
        &parameters.comm_crh_sig_pp,
        &genesis_record,
        &genesis_address.secret_key,
    )
    .unwrap();
    let genesis_memo = [0u8; 32];
    let mut ledger = MerkleTreeIdealLedger::new(
        ledger_parameters,
        genesis_record.commitment(),
        genesis_sn.clone(),
        genesis_memo,
    );

    // Create address 1
    let metadata1 = [1u8; 32];
    let address1 = <InstantiatedDPC as DPCScheme<MerkleTreeIdealLedger>>::create_address(
        &parameters,
        &metadata1,
        rng,
    )
    .unwrap();

    let test_payload = [1u8; 32];
    let test_predicate = Predicate::new(genesis_pred_vk_bytes.clone());
    let test_predicate_private_input = PrivatePredInput {
        vk: parameters.pred_nizk_pp.vk.clone(),
        proof: parameters.pred_nizk_pp.proof.clone(),
    };

    println!("Execution 0");
    let old_records1 = vec![genesis_record.clone(); NUM_INPUT_RECORDS];
    let old_asks1 = vec![genesis_address.secret_key.clone(); NUM_INPUT_RECORDS];
    let old_private_predicate_input = vec![test_predicate_private_input.clone(); NUM_INPUT_RECORDS];

    let new_apks1 = vec![address1.public_key.clone(); NUM_OUTPUT_RECORDS];
    let new_payloads1 = vec![test_payload.clone(); NUM_OUTPUT_RECORDS];
    let new_dummy_flags_1 = vec![false; NUM_OUTPUT_RECORDS];
    let new_birth_predicates1 = vec![test_predicate.clone(); NUM_OUTPUT_RECORDS];
    let new_death_predicates1 = vec![test_predicate.clone(); NUM_OUTPUT_RECORDS];
    let new_private_predicate_input =
        vec![test_predicate_private_input.clone(); NUM_OUTPUT_RECORDS];
    let auxiliary1 = [1u8; 32];

    let (_records1, transaction1) = InstantiatedDPC::execute(
        &parameters,
        &old_records1,
        &old_asks1,
        &old_private_predicate_input,
        &new_apks1,
        &new_dummy_flags_1,
        &new_payloads1,
        &new_birth_predicates1,
        &new_death_predicates1,
        &new_private_predicate_input,
        &auxiliary1,
        &[0u8; 32],
        &ledger,
        rng,
    )
    .unwrap();

    let _is_valid1 = InstantiatedDPC::verify(&parameters, &transaction1, &ledger).unwrap();

    ledger.push(transaction1).unwrap();

    assert_eq!(ledger.len(), 1);
}
