use crate::dpc::delegable_dpc::DelegableDPCComponents;
use crypto_primitives::{CommitmentScheme, FixedLengthCRH, SignatureScheme, NIZK};

#[derive(Derivative)]
#[derivative(Clone(bound = "C: DelegableDPCComponents"))]
pub struct CommCRHSigPublicParameters<C: DelegableDPCComponents> {
    pub addr_comm_pp: <C::AddrC as CommitmentScheme>::Parameters,
    pub rec_comm_pp: <C::RecC as CommitmentScheme>::Parameters,
    pub pred_vk_comm_pp: <C::PredVkComm as CommitmentScheme>::Parameters,
    pub local_data_comm_pp: <C::LocalDataComm as CommitmentScheme>::Parameters,

    pub sn_nonce_crh_pp: <C::SnNonceH as FixedLengthCRH>::Parameters,
    pub pred_vk_crh_pp: <C::PredVkH as FixedLengthCRH>::Parameters,

    pub sig_pp: <C::S as SignatureScheme>::Parameters,
}

#[derive(Derivative)]
#[derivative(Clone(bound = "C: DelegableDPCComponents"))]
pub struct PredNIZKParameters<C: DelegableDPCComponents> {
    pub pk: <C::PredicateNIZK as NIZK>::ProvingParameters,
    pub vk: <C::PredicateNIZK as NIZK>::VerificationParameters,
    pub proof: <C::PredicateNIZK as NIZK>::Proof,
}

pub struct PublicParameters<C: DelegableDPCComponents> {
    pub comm_crh_sig_pp: CommCRHSigPublicParameters<C>,
    pub pred_nizk_pp: PredNIZKParameters<C>,
    pub proof_check_nizk_pp: (
        <C::ProofCheckNIZK as NIZK>::ProvingParameters,
        <C::ProofCheckNIZK as NIZK>::PreparedVerificationParameters,
    ),
    pub core_nizk_pp: (
        <C::MainNIZK as NIZK>::ProvingParameters,
        <C::MainNIZK as NIZK>::PreparedVerificationParameters,
    ),
}

impl<C: DelegableDPCComponents> PublicParameters<C> {
    pub fn core_check_nizk_pp(
        &self,
    ) -> &(
        <C::MainNIZK as NIZK>::ProvingParameters,
        <C::MainNIZK as NIZK>::PreparedVerificationParameters,
    ) {
        &self.core_nizk_pp
    }

    pub fn proof_check_nizk_pp(
        &self,
    ) -> &(
        <C::ProofCheckNIZK as NIZK>::ProvingParameters,
        <C::ProofCheckNIZK as NIZK>::PreparedVerificationParameters,
    ) {
        &self.proof_check_nizk_pp
    }

    pub fn pred_nizk_pp(&self) -> &PredNIZKParameters<C> {
        &self.pred_nizk_pp
    }

    pub fn sn_nonce_crh_pp(&self) -> &<C::SnNonceH as FixedLengthCRH>::Parameters {
        &self.comm_crh_sig_pp.sn_nonce_crh_pp
    }

    pub fn pred_vk_crh_pp(&self) -> &<C::PredVkH as FixedLengthCRH>::Parameters {
        &self.comm_crh_sig_pp.pred_vk_crh_pp
    }

    pub fn local_data_comm_pp(&self) -> &<C::LocalDataComm as CommitmentScheme>::Parameters {
        &self.comm_crh_sig_pp.local_data_comm_pp
    }

    pub fn addr_comm_pp(&self) -> &<C::AddrC as CommitmentScheme>::Parameters {
        &self.comm_crh_sig_pp.addr_comm_pp
    }

    pub fn rec_comm_pp(&self) -> &<C::RecC as CommitmentScheme>::Parameters {
        &self.comm_crh_sig_pp.rec_comm_pp
    }

    pub fn pred_vk_comm_pp(&self) -> &<C::PredVkComm as CommitmentScheme>::Parameters {
        &self.comm_crh_sig_pp.pred_vk_comm_pp
    }

    pub fn sig_pp(&self) -> &<C::S as SignatureScheme>::Parameters {
        &self.comm_crh_sig_pp.sig_pp
    }
}
