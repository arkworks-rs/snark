use crate::dpc::{
    plain_dpc::{
        address::AddressSecretKey, parameters::CommAndCRHPublicParameters,
        predicate::PrivatePredInput, record::DPCRecord, PlainDPCComponents,
    },
    Record,
};
use algebra::{to_bytes, ToConstraintField};
use crypto_primitives::{
    merkle_tree, CommitmentGadget, CommitmentScheme, FixedLengthCRH, FixedLengthCRHGadget,
    NIZKVerifierGadget, PRFGadget, PRF,
};
use r1cs_core::{ConstraintSystemRef, SynthesisError};
use r1cs_std::prelude::*;

use r1cs_std::boolean::Boolean;

pub fn execute_core_checks_gadget<C: PlainDPCComponents>(
    cs: ConstraintSystemRef<C::CoreCheckF>,
    // Parameters
    comm_crh_parameters: &CommAndCRHPublicParameters<C>,
    ledger_parameters: &merkle_tree::Parameters<C::MerkleTreeConfig>,

    // Digest
    ledger_digest: &merkle_tree::Digest<C::MerkleTreeConfig>,

    // Old record stuff
    old_records: &[DPCRecord<C>],
    old_witnesses: &[merkle_tree::Path<C::MerkleTreeConfig>],
    old_address_secret_keys: &[AddressSecretKey<C>],
    old_serial_numbers: &[<C::P as PRF>::Output],

    // New record stuff
    new_records: &[DPCRecord<C>],
    new_sn_nonce_randomness: &[[u8; 32]],
    new_commitments: &[<C::RecC as CommitmentScheme>::Output],

    // Rest
    predicate_comm: &<C::PredVkComm as CommitmentScheme>::Output,
    predicate_rand: &<C::PredVkComm as CommitmentScheme>::Randomness,
    local_data_comm: &<C::LocalDataComm as CommitmentScheme>::Output,
    local_data_rand: &<C::LocalDataComm as CommitmentScheme>::Randomness,
    memo: &[u8; 32],
    auxiliary: &[u8; 32],
) -> Result<(), SynthesisError> {
    execute_core_checks_gadget_helper::<
        C,
        C::AddrC,
        C::RecC,
        C::SnNonceH,
        C::P,
        C::AddrCGadget,
        C::RecCGadget,
        C::SnNonceHGadget,
        C::PGadget,
    >(
        cs,
        //
        comm_crh_parameters,
        ledger_parameters,
        //
        ledger_digest,
        //
        old_records,
        old_witnesses,
        old_address_secret_keys,
        old_serial_numbers,
        //
        new_records,
        new_sn_nonce_randomness,
        new_commitments,
        //
        predicate_comm,
        predicate_rand,
        local_data_comm,
        local_data_rand,
        memo,
        auxiliary,
    )
}

fn execute_core_checks_gadget_helper<
    C,
    AddrC,
    RecC,
    SnNonceH,
    P,
    AddrCGadget,
    RecCGadget,
    SnNonceHGadget,
    PGadget,
>(
    cs: ConstraintSystemRef<C::CoreCheckF>,

    //
    comm_crh_parameters: &CommAndCRHPublicParameters<C>,
    ledger_parameters: &merkle_tree::Parameters<C::MerkleTreeConfig>,

    //
    ledger_digest: &merkle_tree::Digest<C::MerkleTreeConfig>,

    //
    old_records: &[DPCRecord<C>],
    old_witnesses: &[merkle_tree::Path<C::MerkleTreeConfig>],
    old_address_secret_keys: &[AddressSecretKey<C>],
    old_serial_numbers: &[P::Output],

    //
    new_records: &[DPCRecord<C>],
    new_sn_nonce_randomness: &[[u8; 32]],
    new_commitments: &[RecC::Output],

    //
    predicate_comm: &<C::PredVkComm as CommitmentScheme>::Output,
    predicate_rand: &<C::PredVkComm as CommitmentScheme>::Randomness,
    local_data_comm: &<C::LocalDataComm as CommitmentScheme>::Output,
    local_data_rand: &<C::LocalDataComm as CommitmentScheme>::Randomness,
    memo: &[u8; 32],
    auxiliary: &[u8; 32],
) -> Result<(), SynthesisError>
where
    C: PlainDPCComponents<
        AddrC = AddrC,
        RecC = RecC,
        SnNonceH = SnNonceH,
        P = P,
        AddrCGadget = AddrCGadget,
        SnNonceHGadget = SnNonceHGadget,
        RecCGadget = RecCGadget,
        PGadget = PGadget,
    >,
    AddrC: CommitmentScheme,
    RecC: CommitmentScheme,
    SnNonceH: FixedLengthCRH,
    P: PRF,
    RecC::Output: Eq,
    AddrCGadget: CommitmentGadget<AddrC, C::CoreCheckF>,
    RecCGadget: CommitmentGadget<RecC, C::CoreCheckF>,
    SnNonceHGadget: FixedLengthCRHGadget<SnNonceH, C::CoreCheckF>,
    PGadget: PRFGadget<P, C::CoreCheckF>,
{
    let mut old_sns = Vec::with_capacity(old_records.len());
    let mut old_rec_comms = Vec::with_capacity(old_records.len());
    let mut old_apks = Vec::with_capacity(old_records.len());
    let mut old_dummy_flags = Vec::with_capacity(old_records.len());
    let mut old_payloads = Vec::with_capacity(old_records.len());
    let mut old_birth_pred_hashes = Vec::with_capacity(old_records.len());
    let mut old_death_pred_hashes = Vec::with_capacity(old_records.len());

    let mut new_rec_comms = Vec::with_capacity(new_records.len());
    let mut new_apks = Vec::with_capacity(new_records.len());
    let mut new_dummy_flags = Vec::with_capacity(new_records.len());
    let mut new_payloads = Vec::with_capacity(new_records.len());
    let mut new_death_pred_hashes = Vec::with_capacity(new_records.len());
    let mut new_birth_pred_hashes = Vec::with_capacity(new_records.len());

    // Order for new_witnessation of input:
    // 1. addr_comm_pp.
    // 2. rec_comm_pp.
    // 3. crh_pp.
    // 4. ledger_parameters.
    // 5. ledger_digest.
    // 6. for i in 0..NUM_INPUT_RECORDS: old_serial_numbers[i].
    // 7. for j in 0..NUM_OUTPUT_RECORDS: new_commitments[i].
    let (
        addr_comm_pp,
        rec_comm_pp,
        pred_vk_comm_pp,
        local_data_comm_pp,
        sn_nonce_crh_pp,
        ledger_pp,
    ) = {
        let _ns = r1cs_core::ns!(cs, "Declare Comm and CRH parameters");
        let addr_comm_pp = AddrCGadget::ParametersVar::new_constant(
            r1cs_core::ns!(cs, "Declare Addr Comm parameters"),
            &comm_crh_parameters.addr_comm_pp,
        )?;

        let rec_comm_pp = RecCGadget::ParametersVar::new_constant(
            r1cs_core::ns!(cs, "Declare Rec Comm parameters"),
            &comm_crh_parameters.rec_comm_pp,
        )?;

        let local_data_comm_pp =
            <C::LocalDataCommGadget as CommitmentGadget<_, _>>::ParametersVar::new_constant(
                r1cs_core::ns!(cs, "Declare Local Data Comm parameters"),
                &comm_crh_parameters.local_data_comm_pp,
            )?;

        let pred_vk_comm_pp =
            <C::PredVkCommGadget as CommitmentGadget<_, C::CoreCheckF>>::ParametersVar::new_constant(
                r1cs_core::ns!(cs, "Declare Pred Vk COMM parameters"),
                &comm_crh_parameters.pred_vk_comm_pp,
            )?;

        let sn_nonce_crh_pp = SnNonceHGadget::ParametersVar::new_constant(
            r1cs_core::ns!(cs, "Declare SN Nonce CRH parameters"),
            &comm_crh_parameters.sn_nonce_crh_pp,
        )?;

        let ledger_pp =
            <C::MerkleTreeHGadget as FixedLengthCRHGadget<_, _>>::ParametersVar::new_constant(
                r1cs_core::ns!(cs, "Declare Ledger Parameters"),
                ledger_parameters,
            )?;
        (
            addr_comm_pp,
            rec_comm_pp,
            pred_vk_comm_pp,
            local_data_comm_pp,
            sn_nonce_crh_pp,
            ledger_pp,
        )
    };

    let digest_gadget = <C::MerkleTreeHGadget as FixedLengthCRHGadget<_, _>>::OutputVar::new_input(
        r1cs_core::ns!(cs, "Declare ledger digest"),
        || Ok(ledger_digest),
    )?;

    for (((record, witness), secret_key), given_serial_number) in old_records
        .iter()
        .zip(old_witnesses)
        .zip(old_address_secret_keys)
        .zip(old_serial_numbers)
    {
        let _input_ns = r1cs_core::ns!(cs, "Process input record");
        // Declare record contents
        let (
            given_apk,
            given_commitment,
            given_is_dummy,
            given_payload,
            given_birth_pred_hash,
            given_death_pred_hash,
            given_comm_rand,
            sn_nonce,
        ) = {
            let _declare_ns = r1cs_core::ns!(cs, "Declare input record");
            // No need to check that commitments, public keys and hashes are in
            // prime order subgroup because the commitment and CRH parameters
            // are trusted, and so when we recompute these, the newly computed
            // values will always be in correct subgroup. If the input cm, pk
            // or hash is incorrect, then it will not match the computed equivalent.
            let given_apk =
                AddrCGadget::OutputVar::new_witness(r1cs_core::ns!(cs, "Addr PubKey"), || {
                    Ok(&record.address_public_key().public_key)
                })?;
            old_apks.push(given_apk.clone());

            let given_commitment =
                RecCGadget::OutputVar::new_witness(r1cs_core::ns!(cs, "Commitment"), || {
                    Ok(record.commitment())
                })?;
            old_rec_comms.push(given_commitment.clone());

            let given_is_dummy =
                Boolean::new_witness(r1cs_core::ns!(cs, "is_dummy"), || Ok(record.is_dummy()))?;
            old_dummy_flags.push(given_is_dummy.clone());

            let given_payload =
                UInt8::new_witness_vec(r1cs_core::ns!(cs, "Payload"), record.payload())?;
            old_payloads.push(given_payload.clone());

            let given_birth_pred_hash = UInt8::new_witness_vec(
                r1cs_core::ns!(cs, "Birth predicate"),
                &record.birth_predicate_repr(),
            )?;
            old_birth_pred_hashes.push(given_birth_pred_hash.clone());

            let given_death_pred_hash = UInt8::new_witness_vec(
                r1cs_core::ns!(cs, "Death predicate"),
                &record.death_predicate_repr(),
            )?;
            old_death_pred_hashes.push(given_death_pred_hash.clone());

            let given_comm_rand = RecCGadget::RandomnessVar::new_witness(
                r1cs_core::ns!(cs, "Commitment randomness"),
                || Ok(record.commitment_randomness()),
            )?;

            let sn_nonce =
                SnNonceHGadget::OutputVar::new_witness(r1cs_core::ns!(cs, "Sn nonce"), || {
                    Ok(record.serial_number_nonce())
                })?;
            (
                given_apk,
                given_commitment,
                given_is_dummy,
                given_payload,
                given_birth_pred_hash,
                given_death_pred_hash,
                given_comm_rand,
                sn_nonce,
            )
        };

        // ********************************************************************
        // Check that the commitment appears on the ledger,
        // i.e., the membership witness is valid with respect to the
        // transaction set digest.
        // ********************************************************************
        {
            let _witness_ns = r1cs_core::ns!(cs, "Check membership witness");

            let witness =
                merkle_tree::constraints::PathVar::<_, C::MerkleTreeHGadget, _>::new_witness(
                    r1cs_core::ns!(cs, "Declare witness"),
                    || Ok(witness),
                )?;

            witness
                .check_membership(&ledger_pp, &digest_gadget, &given_commitment)?
                .enforce_equal(&given_is_dummy.not())?;
            // if it's not dummy, then the membership check can return false.
        }
        // ********************************************************************

        // ********************************************************************
        // Check that the address public key and secret key form a valid key
        // pair.
        // ********************************************************************

        let sk_prf = {
            // Declare variables for addr_sk contents.
            let _address_ns = r1cs_core::ns!(cs, "Check address keypair");
            let sk_prf = PGadget::new_seed(cs.clone(), &secret_key.sk_prf);
            let metadata = UInt8::new_witness_vec(
                r1cs_core::ns!(cs, "Declare metadata"),
                &secret_key.metadata,
            )?;
            let r_pk = AddrCGadget::RandomnessVar::new_witness(
                r1cs_core::ns!(cs, "Declare r_pk"),
                || Ok(&secret_key.r_pk),
            )?;

            let mut apk_input = sk_prf.clone();
            apk_input.extend_from_slice(&metadata);

            let candidate_apk = AddrCGadget::commit(&addr_comm_pp, &apk_input, &r_pk)?;

            candidate_apk.enforce_equal(&given_apk)?;
            sk_prf
        };
        // ********************************************************************

        // ********************************************************************
        // Check that the serial number is derived correctly.
        // ********************************************************************
        let sn_nonce_bytes = {
            let _sn_ns = r1cs_core::ns!(cs, "Check that sn is derived correctly");

            let sn_nonce_bytes = sn_nonce.to_bytes()?;

            let prf_seed = sk_prf;
            let candidate_serial_number = PGadget::evaluate(&prf_seed, &sn_nonce_bytes)?;

            let given_sn = PGadget::OutputVar::new_input(
                r1cs_core::ns!(cs, "Declare given serial number"),
                || Ok(given_serial_number),
            )?;

            candidate_serial_number.enforce_equal(&given_sn)?;
            old_sns.push(candidate_serial_number);
            sn_nonce_bytes
        };
        // ********************************************************************

        // Check that the record is well-formed.
        {
            let _comm_ns = r1cs_core::ns!(cs, "Check that record is well-formed");
            let apk_bytes = given_apk.to_bytes()?;
            let is_dummy_bytes = given_is_dummy.to_bytes()?;

            let mut comm_input = Vec::new();
            comm_input.extend_from_slice(&apk_bytes);
            comm_input.extend_from_slice(&is_dummy_bytes);
            comm_input.extend_from_slice(&given_payload);
            comm_input.extend_from_slice(&given_birth_pred_hash);
            comm_input.extend_from_slice(&given_death_pred_hash);
            comm_input.extend_from_slice(&sn_nonce_bytes);
            let candidate_commitment =
                RecCGadget::commit(&rec_comm_pp, &comm_input, &given_comm_rand)?;
            candidate_commitment.enforce_equal(&given_commitment)?;
        }
    }

    let sn_nonce_input = {
        let _ns = r1cs_core::ns!(cs, "Convert input serial numbers to bytes");
        let mut sn_nonce_input = Vec::new();
        for (_i, old_sn) in old_sns.iter().enumerate() {
            let bytes = old_sn.to_bytes()?;
            sn_nonce_input.extend_from_slice(&bytes);
        }
        sn_nonce_input
    };

    for (j, ((record, sn_nonce_rand), commitment)) in new_records
        .iter()
        .zip(new_sn_nonce_randomness)
        .zip(new_commitments)
        .enumerate()
    {
        let _ns = r1cs_core::ns!(cs, "Process output record");
        let j = j as u8;

        let (
            given_apk,
            given_record_comm,
            given_comm,
            given_is_dummy,
            given_payload,
            given_birth_pred_hash,
            given_death_pred_hash,
            given_comm_rand,
            sn_nonce,
        ) = {
            let _declare_ns = r1cs_core::ns!(cs, "Declare output record");
            let given_apk =
                AddrCGadget::OutputVar::new_witness(r1cs_core::ns!(cs, "Addr PubKey"), || {
                    Ok(&record.address_public_key().public_key)
                })?;
            new_apks.push(given_apk.clone());
            let given_record_comm = RecCGadget::OutputVar::new_witness(
                r1cs_core::ns!(cs, "Record Commitment"),
                || Ok(record.commitment()),
            )?;
            new_rec_comms.push(given_record_comm.clone());
            let given_comm =
                RecCGadget::OutputVar::new_input(r1cs_core::ns!(cs, "Given Commitment"), || {
                    Ok(commitment)
                })?;

            let given_is_dummy =
                Boolean::new_witness(r1cs_core::ns!(cs, "is_dummy"), || Ok(record.is_dummy()))?;
            new_dummy_flags.push(given_is_dummy.clone());

            let given_payload =
                UInt8::new_witness_vec(r1cs_core::ns!(cs, "Payload"), record.payload())?;
            new_payloads.push(given_payload.clone());

            let given_birth_pred_hash = UInt8::new_witness_vec(
                r1cs_core::ns!(cs, "Birth predicate"),
                &record.birth_predicate_repr(),
            )?;
            new_birth_pred_hashes.push(given_birth_pred_hash.clone());
            let given_death_pred_hash = UInt8::new_witness_vec(
                r1cs_core::ns!(cs, "Death predicate"),
                &record.death_predicate_repr(),
            )?;
            new_death_pred_hashes.push(given_death_pred_hash.clone());

            let given_comm_rand = RecCGadget::RandomnessVar::new_witness(
                r1cs_core::ns!(cs, "Commitment randomness"),
                || Ok(record.commitment_randomness()),
            )?;

            let sn_nonce =
                SnNonceHGadget::OutputVar::new_witness(r1cs_core::ns!(cs, "Sn nonce"), || {
                    Ok(record.serial_number_nonce())
                })?;

            (
                given_apk,
                given_record_comm,
                given_comm,
                given_is_dummy,
                given_payload,
                given_birth_pred_hash,
                given_death_pred_hash,
                given_comm_rand,
                sn_nonce,
            )
        };

        // *******************************************************************
        // Check that the serial number nonce is computed correctly.
        // *******************************************************************
        {
            let _sn_ns = r1cs_core::ns!(cs, "Check that serial number nonce is computed correctly");

            let cur_record_num = UInt8::constant(j);
            let mut cur_record_num_bytes_le = vec![cur_record_num];

            let sn_nonce_randomness = UInt8::new_witness_vec(
                r1cs_core::ns!(cs, "Allocate serial number nonce randomness"),
                sn_nonce_rand,
            )?;
            cur_record_num_bytes_le.extend_from_slice(&sn_nonce_randomness);
            cur_record_num_bytes_le.extend_from_slice(&sn_nonce_input);

            let sn_nonce_input = cur_record_num_bytes_le;

            let candidate_sn_nonce = SnNonceHGadget::evaluate(&sn_nonce_crh_pp, &sn_nonce_input)?;
            candidate_sn_nonce.enforce_equal(&sn_nonce)?;
        }
        // *******************************************************************

        // *******************************************************************
        // Check that the record is well-formed.
        // *******************************************************************
        {
            let _comm_ns = r1cs_core::ns!(cs, "Check that record is well-formed");
            let apk_bytes = given_apk.to_bytes()?;
            let is_dummy_bytes = given_is_dummy.to_bytes()?;
            let sn_nonce_bytes = sn_nonce.to_bytes()?;

            let mut comm_input = Vec::new();
            comm_input.extend_from_slice(&apk_bytes);
            comm_input.extend_from_slice(&is_dummy_bytes);
            comm_input.extend_from_slice(&given_payload);
            comm_input.extend_from_slice(&given_birth_pred_hash);
            comm_input.extend_from_slice(&given_death_pred_hash);
            comm_input.extend_from_slice(&sn_nonce_bytes);

            let candidate_commitment =
                RecCGadget::commit(&rec_comm_pp, &comm_input, &given_comm_rand)?;
            candidate_commitment.enforce_equal(&given_comm)?;
            candidate_commitment.enforce_equal(&given_record_comm)?;
        }
    }
    // *******************************************************************
    // Check that predicate commitment is well formed.
    // *******************************************************************
    {
        let _comm_ns = r1cs_core::ns!(cs, "Check that predicate commitment is well-formed");

        let mut input = Vec::new();
        for i in 0..C::NUM_INPUT_RECORDS {
            input.extend_from_slice(&old_death_pred_hashes[i]);
        }

        for j in 0..C::NUM_OUTPUT_RECORDS {
            input.extend_from_slice(&new_birth_pred_hashes[j]);
        }

        let given_comm_rand =
        <C::PredVkCommGadget as CommitmentGadget<_, C::CoreCheckF>>::RandomnessVar::new_witness(
            r1cs_core::ns!(cs, "Commitment randomness"),
            || Ok(predicate_rand),
        )?;

        let given_comm =
            <C::PredVkCommGadget as CommitmentGadget<_, C::CoreCheckF>>::OutputVar::new_input(
                r1cs_core::ns!(cs, "Commitment output"),
                || Ok(predicate_comm),
            )?;

        let candidate_commitment =
            <C::PredVkCommGadget as CommitmentGadget<_, C::CoreCheckF>>::commit(
                &pred_vk_comm_pp,
                &input,
                &given_comm_rand,
            )?;

        candidate_commitment.enforce_equal(&given_comm)?;
    }
    {
        let _ns = r1cs_core::ns!(cs, "Check that local data commitment is valid.");

        let mut local_data_bytes = Vec::new();
        for i in 0..C::NUM_INPUT_RECORDS {
            let _ns = r1cs_core::ns!(cs, "Construct local data with Input Record");
            local_data_bytes.extend_from_slice(&old_rec_comms[i].to_bytes()?);
            local_data_bytes.extend_from_slice(&old_apks[i].to_bytes()?);
            local_data_bytes.extend_from_slice(&old_dummy_flags[i].to_bytes()?);
            local_data_bytes.extend_from_slice(&old_payloads[i]);
            local_data_bytes.extend_from_slice(&old_birth_pred_hashes[i]);
            local_data_bytes.extend_from_slice(&old_death_pred_hashes[i]);
            local_data_bytes.extend_from_slice(&old_sns[i].to_bytes()?);
        }

        for j in 0..C::NUM_OUTPUT_RECORDS {
            let _ns = r1cs_core::ns!(cs, "Construct local data with Output Record");
            local_data_bytes.extend_from_slice(&new_rec_comms[j].to_bytes()?);
            local_data_bytes.extend_from_slice(&new_apks[j].to_bytes()?);
            local_data_bytes.extend_from_slice(&new_dummy_flags[j].to_bytes()?);
            local_data_bytes.extend_from_slice(&new_payloads[j]);
            local_data_bytes.extend_from_slice(&new_birth_pred_hashes[j]);
            local_data_bytes.extend_from_slice(&new_death_pred_hashes[j]);
        }
        let memo = UInt8::new_input_vec(r1cs_core::ns!(cs, "Allocate memorandum"), memo)?;
        local_data_bytes.extend_from_slice(&memo);

        let auxiliary =
            UInt8::new_witness_vec(r1cs_core::ns!(cs, "Allocate auxiliary input"), auxiliary)?;
        local_data_bytes.extend_from_slice(&auxiliary);

        let local_data_comm_rand =
            <C::LocalDataCommGadget as CommitmentGadget<_, _>>::RandomnessVar::new_witness(
                r1cs_core::ns!(cs, "Allocate local data commitment randomness"),
                || Ok(local_data_rand),
            )?;

        let declared_local_data_comm =
            <C::LocalDataCommGadget as CommitmentGadget<_, _>>::OutputVar::new_input(
                r1cs_core::ns!(cs, "Allocate local data commitment"),
                || Ok(local_data_comm),
            )?;

        let comm = C::LocalDataCommGadget::commit(
            &local_data_comm_pp,
            &local_data_bytes,
            &local_data_comm_rand,
        )?;

        comm.enforce_equal(&declared_local_data_comm)?;
    }
    Ok(())
}

pub fn execute_proof_check_gadget<C: PlainDPCComponents>(
    cs: ConstraintSystemRef<C::ProofCheckF>,
    // Parameters
    comm_crh_parameters: &CommAndCRHPublicParameters<C>,

    // Old record death predicate verif. keys and proofs
    old_death_pred_vk_and_pf: &[PrivatePredInput<C>],

    // New record birth predicate verif. keys and proofs
    new_birth_pred_vk_and_pf: &[PrivatePredInput<C>],

    // Rest
    predicate_comm: &<C::PredVkComm as CommitmentScheme>::Output,
    predicate_rand: &<C::PredVkComm as CommitmentScheme>::Randomness,

    local_data_comm: &<C::LocalDataComm as CommitmentScheme>::Output,
) -> Result<(), SynthesisError>
where
    <C::LocalDataComm as CommitmentScheme>::Output: ToConstraintField<C::CoreCheckF>,
    <C::LocalDataComm as CommitmentScheme>::Parameters: ToConstraintField<C::CoreCheckF>,
{
    // Declare public parameters.
    let (pred_vk_comm_pp, pred_vk_crh_pp) = {
        let _ns = r1cs_core::ns!(cs, "Declare Comm and CRH parameters");

        let pred_vk_comm_pp =
            <C::PredVkCommGadget as CommitmentGadget<_, C::ProofCheckF>>::ParametersVar::new_constant(
                r1cs_core::ns!(cs, "Declare Pred Vk COMM parameters"),
                &comm_crh_parameters.pred_vk_comm_pp,
            )?;

        let pred_vk_crh_pp = <C::PredVkHGadget as FixedLengthCRHGadget<_, C::ProofCheckF>>::ParametersVar::new_constant(
            r1cs_core::ns!(cs, "Declare Pred Vk CRH parameters"),
            &comm_crh_parameters.pred_vk_crh_pp,
        )?;

        (pred_vk_comm_pp, pred_vk_crh_pp)
    };

    // ************************************************************************
    // Construct predicate input
    // ************************************************************************

    // First we convert the input for the predicates into `CoreCheckF` field
    // elements
    let local_data_comm_pp_fe = ToConstraintField::<C::CoreCheckF>::to_field_elements(
        &comm_crh_parameters.local_data_comm_pp,
    )
    .map_err(|_| SynthesisError::AssignmentMissing)?;

    let local_data_comm_fe = ToConstraintField::<C::CoreCheckF>::to_field_elements(local_data_comm)
        .map_err(|_| SynthesisError::AssignmentMissing)?;

    // Then we convert these field elements into bytes
    let pred_input = [
        to_bytes![local_data_comm_pp_fe].map_err(|_| SynthesisError::AssignmentMissing)?,
        to_bytes![local_data_comm_fe].map_err(|_| SynthesisError::AssignmentMissing)?,
    ];

    let pred_input_bytes = [
        UInt8::new_input_vec(
            r1cs_core::ns!(cs, "Allocate local data pp "),
            &pred_input[0],
        )?,
        UInt8::new_input_vec(
            r1cs_core::ns!(cs, "Allocate local data comm"),
            &pred_input[1],
        )?,
    ];

    let pred_input_bits = [
        pred_input_bytes[0].to_bits_le()?,
        pred_input_bytes[1].to_bits_le()?,
    ];
    // ************************************************************************
    // ************************************************************************

    let mut old_death_pred_hashes = Vec::new();
    let mut new_birth_pred_hashes = Vec::new();
    for i in 0..C::NUM_INPUT_RECORDS {
        let _ns = r1cs_core::ns!(cs, "Check death predicate for input record");

        let death_pred_proof =
            <C::PredicateNIZKGadget as NIZKVerifierGadget<_, _>>::ProofVar::new_witness(
                r1cs_core::ns!(cs, "Allocate proof"),
                || Ok(&old_death_pred_vk_and_pf[i].proof),
            )?;

        let death_pred_vk =
            <C::PredicateNIZKGadget as NIZKVerifierGadget<_, _>>::new_verification_key_unchecked(
                r1cs_core::ns!(cs, "Allocate verification key"),
                || Ok(&old_death_pred_vk_and_pf[i].vk),
                AllocationMode::Witness,
            )?;

        let _ns = r1cs_core::ns!(cs, "Convert vk to bytes");
        let death_pred_vk_bytes = death_pred_vk.to_bytes()?;
        drop(_ns);

        let _ns = r1cs_core::ns!(cs, "Evaluate death predicate hash");
        let claimed_death_pred_hash =
            C::PredVkHGadget::evaluate(&pred_vk_crh_pp, &death_pred_vk_bytes)?;
        drop(_ns);

        let _ns = r1cs_core::ns!(cs, "Convert death predicate hash to bytes");
        let claimed_death_pred_hash_bytes = claimed_death_pred_hash.to_bytes()?;
        drop(_ns);

        old_death_pred_hashes.push(claimed_death_pred_hash_bytes);

        let _ns = r1cs_core::ns!(cs, "Verify death predicate");
        let position = UInt8::constant(i as u8).to_bits_le()?;
        C::PredicateNIZKGadget::verify(
            &death_pred_vk,
            ([position].iter())
                .chain(&pred_input_bits)
                .filter(|v| !v.is_empty()),
            &death_pred_proof,
        )?
        .enforce_equal(&Boolean::TRUE)?;
    }

    for j in 0..C::NUM_OUTPUT_RECORDS {
        let _ns = r1cs_core::ns!(cs, "Check birth predicate for output record");

        let birth_pred_proof =
            <C::PredicateNIZKGadget as NIZKVerifierGadget<_, _>>::ProofVar::new_witness(
                r1cs_core::ns!(cs, "Allocate proof"),
                || Ok(&new_birth_pred_vk_and_pf[j].proof),
            )?;

        let birth_pred_vk =
            <C::PredicateNIZKGadget as NIZKVerifierGadget<_, _>>::new_verification_key_unchecked(
                r1cs_core::ns!(cs, "Allocate verification key"),
                || Ok(&new_birth_pred_vk_and_pf[j].vk),
                AllocationMode::Witness,
            )?;

        let _ns = r1cs_core::ns!(cs, "Convert vk to bytes");
        let birth_pred_vk_bytes = birth_pred_vk.to_bytes()?;
        drop(_ns);

        let _ns = r1cs_core::ns!(cs, "Evaluate birth predicate hash");
        let claimed_birth_pred_hash =
            C::PredVkHGadget::evaluate(&pred_vk_crh_pp, &birth_pred_vk_bytes)?;
        drop(_ns);

        let _ns = r1cs_core::ns!(cs, "Convert birth predicate hash to bytes");
        let claimed_birth_pred_hash_bytes = claimed_birth_pred_hash.to_bytes()?;
        drop(_ns);

        new_birth_pred_hashes.push(claimed_birth_pred_hash_bytes);

        let _ns = r1cs_core::ns!(cs, "Verify birth predicate");
        let position = UInt8::constant(j as u8).to_bits_le()?;
        C::PredicateNIZKGadget::verify(
            &birth_pred_vk,
            ([position].iter())
                .chain(&pred_input_bits)
                .filter(|v| !v.is_empty()),
            &birth_pred_proof,
        )?
        .enforce_equal(&Boolean::TRUE)?;
    }
    {
        let _comm_ns = r1cs_core::ns!(cs, "Check that predicate commitment is well-formed");

        let mut input = Vec::new();
        for i in 0..C::NUM_INPUT_RECORDS {
            input.extend_from_slice(&old_death_pred_hashes[i]);
        }

        for j in 0..C::NUM_OUTPUT_RECORDS {
            input.extend_from_slice(&new_birth_pred_hashes[j]);
        }

        let given_comm_rand =
        <C::PredVkCommGadget as CommitmentGadget<_, C::ProofCheckF>>::RandomnessVar::new_witness(
            r1cs_core::ns!(cs, "Commitment randomness"),
            || Ok(predicate_rand),
        )?;

        let given_comm =
            <C::PredVkCommGadget as CommitmentGadget<_, C::ProofCheckF>>::OutputVar::new_input(
                r1cs_core::ns!(cs, "Commitment output"),
                || Ok(predicate_comm),
            )?;

        let candidate_commitment =
            <C::PredVkCommGadget as CommitmentGadget<_, C::ProofCheckF>>::commit(
                &pred_vk_comm_pp,
                &input,
                &given_comm_rand,
            )?;

        candidate_commitment.enforce_equal(&given_comm)?;
    }
    Ok(())
}
