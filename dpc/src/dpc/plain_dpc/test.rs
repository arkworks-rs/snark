use super::instantiated::*;
use algebra::{
    bls12_377::{Fq, Fr},
    to_bytes,
};
use tracing_subscriber::layer::SubscriberExt;

#[cfg(debug_assertions)]
use groth16::PreparedVerifyingKey;
use rand::SeedableRng;
use rand_xorshift::XorShiftRng;

use crypto_primitives::FixedLengthCRH;

use r1cs_core::{ConstraintLayer, ConstraintSystem};

use crate::constraints::plain_dpc::{execute_core_checks_gadget, execute_proof_check_gadget};

use crate::dpc::{
    plain_dpc::{predicate::PrivatePredInput, predicate_circuit::*, ExecuteContext, DPC},
    Record,
};

use crate::ledger::Ledger;

#[test]
fn test_execute_constraint_systems() {
    let mut layer = ConstraintLayer::default();
    layer.mode = r1cs_core::TracingMode::OnlyConstraints;
    let subscriber = tracing_subscriber::Registry::default().with(layer);
    tracing::subscriber::set_global_default(subscriber).unwrap();
    let mut rng = XorShiftRng::seed_from_u64(1231275789u64);
    // Generate parameters for the ledger, commitment schemes, CRH, and the
    // "always-accept" predicate.
    let ledger_parameters = MerkleTreeIdealLedger::setup(&mut rng).expect("Ledger setup failed");
    let comm_and_crh_pp = InstantiatedDPC::generate_comm_and_crh_parameters(&mut rng).unwrap();
    let pred_nizk_pp =
        InstantiatedDPC::generate_pred_nizk_parameters(&comm_and_crh_pp, &mut rng).unwrap();
    #[cfg(debug_assertions)]
    let pred_nizk_pvk: PreparedVerifyingKey<_> = pred_nizk_pp.vk.clone().into();

    let pred_nizk_vk_bytes = to_bytes![PredVkCRH::evaluate(
        &comm_and_crh_pp.pred_vk_crh_pp,
        &to_bytes![pred_nizk_pp.vk].unwrap()
    )
    .unwrap()]
    .unwrap();

    // Generate metadata and an address for a dummy initial, or "genesis", record.
    let genesis_metadata = [1u8; 32];
    let genesis_address =
        DPC::create_address_helper(&comm_and_crh_pp, &genesis_metadata, &mut rng).unwrap();
    let genesis_sn_nonce =
        SnNonceCRH::evaluate(&comm_and_crh_pp.sn_nonce_crh_pp, &[0u8; 1]).unwrap();
    let genesis_record = DPC::generate_record(
        &comm_and_crh_pp,
        &genesis_sn_nonce,
        &genesis_address.public_key,
        true,
        &[0u8; 32],
        &Predicate::new(pred_nizk_vk_bytes.clone()),
        &Predicate::new(pred_nizk_vk_bytes.clone()),
        &mut rng,
    )
    .unwrap();

    // Generate serial number for the genesis record.
    let genesis_sn = DPC::generate_sn(&genesis_record, &genesis_address.secret_key).unwrap();
    let genesis_memo = [0u8; 32];

    // Use genesis record, serial number, and memo to initialize the ledger.
    let ledger = MerkleTreeIdealLedger::new(
        ledger_parameters,
        genesis_record.commitment(),
        genesis_sn.clone(),
        genesis_memo,
    );

    // Set the input records for our transaction to be the initial dummy records.
    let old_records = vec![genesis_record; NUM_INPUT_RECORDS];
    let old_asks = vec![genesis_address.secret_key; NUM_INPUT_RECORDS];

    // Construct new records.

    // Create an address for an actual new record.
    let new_metadata = [1u8; 32];
    let new_address =
        DPC::create_address_helper(&comm_and_crh_pp, &new_metadata, &mut rng).unwrap();

    // Create a payload.
    let new_payload = [1u8; 32];
    // Set the new record's predicate to be the "always-accept" predicate.
    let new_predicate = Predicate::new(pred_nizk_vk_bytes);

    let new_apks = vec![new_address.public_key; NUM_OUTPUT_RECORDS];
    let new_payloads = vec![new_payload.clone(); NUM_OUTPUT_RECORDS];
    let new_birth_predicates = vec![new_predicate.clone(); NUM_OUTPUT_RECORDS];
    let new_death_predicates = vec![new_predicate; NUM_OUTPUT_RECORDS];
    let new_dummy_flags = vec![false; NUM_OUTPUT_RECORDS];
    let auxiliary = [122u8; 32];
    let memo = [238u8; 32];

    let context = DPC::execute_helper(
        &comm_and_crh_pp,
        &old_records,
        &old_asks,
        &new_apks,
        &new_dummy_flags,
        &new_payloads,
        &new_birth_predicates,
        &new_death_predicates,
        &memo,
        &auxiliary,
        &ledger,
        &mut rng,
    )
    .unwrap();

    let ExecuteContext {
        comm_and_crh_pp: _comm_and_crh_pp,
        ledger_digest,

        old_records,
        old_witnesses,
        old_address_secret_keys,
        old_serial_numbers,

        new_records,
        new_sn_nonce_randomness,
        new_commitments,

        predicate_comm,
        predicate_rand,

        local_data_comm,
        local_data_rand,
    } = context;

    //////////////////////////////////////////////////////////////////////////
    // Check that the core check constraint system was satisfied.
    let core_cs = ConstraintSystem::<Fr>::new_ref();

    let core_ns = r1cs_core::ns!(core_cs, "Core checks");
    let cs = core_ns.cs();
    execute_core_checks_gadget::<_>(
        cs.clone(),
        &comm_and_crh_pp,
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
        &memo,
        &auxiliary,
    )
    .unwrap();

    if !cs.is_satisfied().unwrap() {
        println!("=========================================================");
        println!("Unsatisfied constraints:");
        println!("{}", core_cs.which_is_unsatisfied().unwrap().unwrap());
        println!("=========================================================");
    }

    if cs.is_satisfied().unwrap() {
        println!("\n\n\n\nAll Core check constraints:");
        for constraint in cs.constraint_names().unwrap() {
            println!("{}", constraint)
        }
    }

    println!("=========================================================");
    println!("=========================================================");
    println!("=========================================================\n\n\n");

    assert!(cs.is_satisfied().unwrap());
    drop(core_ns);

    // Check that the proof check constraint system was satisfied.
    let pf_check_cs = ConstraintSystem::<Fq>::new_ref();

    let mut old_proof_and_vk = vec![];
    for i in 0..NUM_INPUT_RECORDS {
        use crypto_primitives::nizk::NIZK;
        let proof = PredicateNIZK::prove(
            &pred_nizk_pp.pk,
            EmptyPredicateCircuit::new(&comm_and_crh_pp, &local_data_comm, i as u8),
            &mut rng,
        )
        .expect("Proof should work");
        #[cfg(debug_assertions)]
        {
            let pred_pub_input: PredicateLocalData<Components> = PredicateLocalData {
                local_data_comm_pp: comm_and_crh_pp.local_data_comm_pp.clone(),
                local_data_comm: local_data_comm.clone(),
                position: i as u8,
            };
            assert!(
                PredicateNIZK::verify(&pred_nizk_pvk, &pred_pub_input, &proof)
                    .expect("Proof should verify")
            );
        }
        let private_input: PrivatePredInput<Components> = PrivatePredInput {
            vk: pred_nizk_pp.vk.clone(),
            proof,
        };
        old_proof_and_vk.push(private_input);
    }

    let mut new_proof_and_vk = vec![];
    for i in 0..NUM_OUTPUT_RECORDS {
        use crypto_primitives::nizk::NIZK;
        let proof = PredicateNIZK::prove(
            &pred_nizk_pp.pk,
            EmptyPredicateCircuit::new(&comm_and_crh_pp, &local_data_comm, i as u8),
            &mut rng,
        )
        .expect("Proof should work");
        let private_input: PrivatePredInput<Components> = PrivatePredInput {
            vk: pred_nizk_pp.vk.clone(),
            proof,
        };
        new_proof_and_vk.push(private_input);
    }

    let pf_check_ns = r1cs_core::ns!(pf_check_cs, "Check predicate proofs");
    let cs = pf_check_ns.cs();
    execute_proof_check_gadget::<_>(
        cs,
        &comm_and_crh_pp,
        &old_proof_and_vk,
        &new_proof_and_vk,
        &predicate_comm,
        &predicate_rand,
        &local_data_comm,
    )
    .unwrap();

    if !pf_check_cs.is_satisfied().unwrap() {
        println!("=========================================================");
        println!("Unsatisfied constraints:");
        println!("{}", pf_check_cs.which_is_unsatisfied().unwrap().unwrap());
        println!("=========================================================");
    }
    println!("\n\n\n\nAll Proof check constraints:");
    if pf_check_cs.is_satisfied().unwrap() {
        for constraint in pf_check_cs.constraint_names().unwrap() {
            println!("{}", constraint)
        }
    }
    println!("=========================================================");
    println!("=========================================================");
    println!("=========================================================");

    assert!(pf_check_cs.is_satisfied().unwrap());
}
