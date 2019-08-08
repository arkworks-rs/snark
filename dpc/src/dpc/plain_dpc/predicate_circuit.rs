use crate::{
    crypto_primitives::{CommitmentScheme, PRF},
    dpc::{plain_dpc::DPCRecord, Record},
    gadgets::Assignment,
    plain_dpc::*,
};
use snark_gadgets::{uint8::UInt8, utils::AllocGadget};
use std::io::{Result as IoResult, Write};

use algebra::{bytes::ToBytes, utils::ToEngineFr, PairingEngine};

use snark::{Circuit, ConstraintSystem, SynthesisError};

use crate::Error;

pub struct PredicateHashInput<C: PlainDPCComponents> {
    pub old_rec_comms:      Vec<<C::RecC as CommitmentScheme>::Output>,
    pub old_apks:           Vec<<C::AddrC as CommitmentScheme>::Output>,
    pub old_dummy_flags:    Vec<bool>,
    pub old_payloads:       Vec<<DPCRecord<C> as Record>::Payload>,
    pub old_death_pred_ids: Vec<Vec<u8>>,
    pub old_birth_pred_ids: Vec<Vec<u8>>,
    pub old_serial_numbers: Vec<<C::P as PRF>::Output>,

    pub new_rec_comms:      Vec<<C::RecC as CommitmentScheme>::Output>,
    pub new_apks:           Vec<<C::AddrC as CommitmentScheme>::Output>,
    pub new_dummy_flags:    Vec<bool>,
    pub new_payloads:       Vec<<DPCRecord<C> as Record>::Payload>,
    pub new_death_pred_ids: Vec<Vec<u8>>,
    pub new_birth_pred_ids: Vec<Vec<u8>>,

    pub memo:      [u8; 32],
    pub auxiliary: [u8; 32],
}

impl<C: PlainDPCComponents> Default for PredicateHashInput<C> {
    fn default() -> Self {
        Self {
            old_rec_comms:      vec![
                <C::RecC as CommitmentScheme>::Output::default();
                C::NUM_INPUT_RECORDS
            ],
            old_apks:           vec![
                <C::AddrC as CommitmentScheme>::Output::default();
                C::NUM_INPUT_RECORDS
            ],
            old_dummy_flags:    vec![false; C::NUM_INPUT_RECORDS],
            old_payloads:       vec![
                <DPCRecord<C> as Record>::Payload::default();
                C::NUM_INPUT_RECORDS
            ],
            old_death_pred_ids: vec![vec![0u8; 48]; C::NUM_INPUT_RECORDS],
            old_birth_pred_ids: vec![vec![0u8; 48]; C::NUM_INPUT_RECORDS],
            old_serial_numbers: vec![<C::P as PRF>::Output::default(); C::NUM_INPUT_RECORDS],

            new_rec_comms:      vec![
                <C::RecC as CommitmentScheme>::Output::default();
                C::NUM_OUTPUT_RECORDS
            ],
            new_apks:           vec![
                <C::AddrC as CommitmentScheme>::Output::default();
                C::NUM_OUTPUT_RECORDS
            ],
            new_dummy_flags:    vec![false; C::NUM_OUTPUT_RECORDS],
            new_payloads:       vec![
                <DPCRecord<C> as Record>::Payload::default();
                C::NUM_OUTPUT_RECORDS
            ],
            new_death_pred_ids: vec![vec![0u8; 48]; C::NUM_OUTPUT_RECORDS],
            new_birth_pred_ids: vec![vec![0u8; 48]; C::NUM_OUTPUT_RECORDS],

            memo:      [0u8; 32],
            auxiliary: [0u8; 32],
        }
    }
}

impl<C: PlainDPCComponents> ToBytes for PredicateHashInput<C> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        for i in 0..C::NUM_INPUT_RECORDS {
            self.old_rec_comms[i].write(&mut writer)?;
            self.old_apks[i].write(&mut writer)?;
            self.old_dummy_flags[i].write(&mut writer)?;
            self.old_payloads[i].write(&mut writer)?;
            self.old_death_pred_ids[i].write(&mut writer)?;
            self.old_birth_pred_ids[i].write(&mut writer)?;
            self.old_serial_numbers[i].write(&mut writer)?;
        }

        for i in 0..C::NUM_OUTPUT_RECORDS {
            self.new_rec_comms[i].write(&mut writer)?;
            self.new_apks[i].write(&mut writer)?;
            self.new_dummy_flags[i].write(&mut writer)?;
            self.new_payloads[i].write(&mut writer)?;
            self.new_death_pred_ids[i].write(&mut writer)?;
            self.new_birth_pred_ids[i].write(&mut writer)?;
        }
        self.memo.write(&mut writer)?;
        self.auxiliary.write(&mut writer)?;
        Ok(())
    }
}

pub struct PredicateLocalData<C: PlainDPCComponents> {
    pub local_data_comm_pp: <C::LocalDataComm as CommitmentScheme>::Parameters,
    pub local_data_comm:    <C::LocalDataComm as CommitmentScheme>::Output,
    pub position:           u8,
}

// Convert each component to bytes and pack into field elements.
impl<C: PlainDPCComponents> ToEngineFr<C::E> for PredicateLocalData<C>
where
    <C::LocalDataComm as CommitmentScheme>::Output: ToEngineFr<C::E>,
    <C::LocalDataComm as CommitmentScheme>::Parameters: ToEngineFr<C::E>,
{
    fn to_engine_fr(&self) -> Result<Vec<<C::E as PairingEngine>::Fr>, Error> {
        let mut v = ToEngineFr::<C::E>::to_engine_fr([self.position].as_ref())?;
        v.extend_from_slice(&self.local_data_comm_pp.to_engine_fr()?);
        v.extend_from_slice(&self.local_data_comm.to_engine_fr()?);
        Ok(v)
    }
}

pub struct EmptyPredicateCircuit<C: PlainDPCComponents> {
    // Parameters
    comm_and_crh_parameters: Option<CommAndCRHPublicParameters<C>>,

    // Commitment to Predicate input.
    local_data_comm: Option<<C::LocalDataComm as CommitmentScheme>::Output>,
    position:        u8,
}

impl<C: PlainDPCComponents> EmptyPredicateCircuit<C> {
    pub fn blank(comm_and_crh_parameters: &CommAndCRHPublicParameters<C>) -> Self {
        let local_data_comm = <C::LocalDataComm as CommitmentScheme>::Output::default();

        Self {
            comm_and_crh_parameters: Some(comm_and_crh_parameters.clone()),
            local_data_comm:         Some(local_data_comm),
            position:                0u8,
        }
    }

    pub fn new(
        comm_amd_crh_parameters: &CommAndCRHPublicParameters<C>,
        local_data_comm: &<C::LocalDataComm as CommitmentScheme>::Output,
        position: u8,
    ) -> Self {
        Self {
            // Parameters
            comm_and_crh_parameters: Some(comm_amd_crh_parameters.clone()),

            // Other stuff
            local_data_comm: Some(local_data_comm.clone()),
            position,
        }
    }
}

impl<C: PlainDPCComponents> Circuit<C::E> for EmptyPredicateCircuit<C> {
    fn synthesize<CS: ConstraintSystem<C::E>>(self, cs: &mut CS) -> Result<(), SynthesisError> {
        let _position = UInt8::alloc_input_vec(cs.ns(|| "Alloc position"), &[self.position])?;

        let _local_data_comm_pp =
            <C::LocalDataCommGadget as CommitmentGadget<_, _>>::ParametersGadget::alloc_input(
                &mut cs.ns(|| "Declare Pred Input Comm parameters"),
                || {
                    self.comm_and_crh_parameters
                        .get()
                        .map(|pp| &pp.local_data_comm_pp)
                },
            )?;

        let _local_data_comm =
            <C::LocalDataCommGadget as CommitmentGadget<_, _>>::OutputGadget::alloc_input(
                cs.ns(|| "Allocate predicate commitment"),
                || self.local_data_comm.get(),
            )?;

        Ok(())
    }
}
