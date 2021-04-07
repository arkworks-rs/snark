use algebra::{AffineCurve, ToConstraintField, UniformRand};
use r1cs_core::ConstraintSynthesizer;
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC, UniversalParams,
        CommitterKey as DLogCommitterKey, VerifierKey as DLogVerifierKey,
    },
    PolynomialCommitment, Error as PCError
};
use crate::darlin::{
    accumulators::{
        ItemAccumulator,
        dlog::{DualDLogItem, DualDLogItemAccumulator},
    },
    pcd::{
        final_darlin::{FinalDarlinPCD, FinalDarlinPCDVerifierKey},
        simple_marlin::{SimpleMarlinPCD, SimpleMarlinPCDVerifierKey},
        error::PCDError,
    },
    data_structures::FinalDarlinDeferredData,
};
use rand::RngCore;
use digest::Digest;
use std::fmt::Debug;

pub mod simple_marlin;
pub mod final_darlin;
pub mod error;

/// Configuration parameters for the PCD scheme: for now, just the size of the
/// committer key to be used throughout the PCD scheme.
pub struct PCDParameters {
    pub segment_size: usize
}

impl PCDParameters {
    /// We assume the DLOG keys to be generated outside the PCD scheme,
    /// so this function actually just trim them to the segment size
    /// specified in the config.
    pub fn universal_setup<G: AffineCurve, D: Digest>(
        &self,
        params: &UniversalParams<G>
    ) -> Result<(DLogCommitterKey<G>, DLogVerifierKey<G>), PCError>
    {
        InnerProductArgPC::<G, D>::trim(
            params,
            self.segment_size - 1,
        )
    }
}

/// Trait for recursive circuit of a PCD scheme. Both witnesses and public inputs
/// are derived from previous proofs and some additional "payload".
pub trait PCDCircuit<G: AffineCurve>: ConstraintSynthesizer<G::ScalarField> {

    /// Any data that may be needed to bootstrap the circuit that is not covered by the other
    /// fields.
    type SetupData: Clone;

    /// Witnesses needed to enforce the business logic of the circuit,
    /// (e.g. all the things not related to recursion).
    type IncrementalData;

    /// Elements that are deferred in recursion. They should be derived from the previous proofs.
    type SystemInputs: ToConstraintField<G::ScalarField> + Debug + Clone;

    /// PCD type the circuit need to verify
    type PreviousPCD:  PCD;

    /// Initialize the circuit state without explicitly assigning inputs and witnesses.
    /// To be used to generate pk and vk.
    fn init(config: Self::SetupData) -> Self;

    /// Assign a concrete state to the circuit, using previous proofs and some "payload".
    /// As the circuit needs to verify previous proofs, it also needs the corresponding vks;
    fn init_state(
        config:               Self::SetupData,
        previous_proofs_data: Vec<Self::PreviousPCD>,
        previous_proofs_vks:  Vec<<Self::PreviousPCD as PCD>::PCDVerifierKey>,
        incremental_data:     Self::IncrementalData,
    ) -> Self;

    /// Extract the system inputs from a concrete instantiation of the circuit.
    /// Return Error if it's not possible to derive SystemInputs.
    fn get_sys_ins(&self) -> Result<&Self::SystemInputs, PCDError>;

    /// Extract the user inputs from a concrete instantiation of the circuit.
    /// Return Error if it's not possible to derive UserInputs.
    fn get_usr_ins(&self) -> Result<Vec<G::ScalarField>, PCDError>;
}

/// This trait expresses the functions for proof carrying data, in which the PCD is assumed
/// to be a set of data consisting of a statement, some deferred elements and a proof.
pub trait PCD: Sized + Send + Sync {
    type PCDAccumulator: ItemAccumulator;
    type PCDVerifierKey: AsRef<<Self::PCDAccumulator as ItemAccumulator>::AccumulatorVerifierKey>;

    /// Perform only cheap verification that tipically includes algebraic checks which
    /// are not MSM, e.g. verification of Marlin's sumcheck equations, Bullet reduction
    /// and so on. Return an accumulator for the proof if verification was succesfull,
    /// Error otherwise.
    fn succinct_verify(
        &self,
        vk:         &Self::PCDVerifierKey,
    ) -> Result<<Self::PCDAccumulator as ItemAccumulator>::Item, PCDError>;

    /// Check the current accumulator, tipically via one or more MSMs
    fn hard_verify<R: RngCore>(
        &self,
        acc:    <Self::PCDAccumulator as ItemAccumulator>::Item,
        vk:     &Self::PCDVerifierKey,
        rng:    &mut R,
    ) -> Result<bool, PCDError>
    { <Self::PCDAccumulator as ItemAccumulator>::check_items::<R>(vk.as_ref(), &[acc], rng)
        .map_err(|e| PCDError::FailedHardVerification(e.to_string()))
    }

    /// Perform full verification of `self`, i.e. both succinct and hard part.
    fn verify<R: RngCore>(
        &self,
        vk:         &Self::PCDVerifierKey,
        rng:        &mut R,
    ) -> Result<bool, PCDError>
    {
        let acc = self.succinct_verify(vk)?;
        self.hard_verify::<R>(acc, vk, rng)
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
/// Achieve polymorphism for PCD via an enumerable. This provides nice APIs for
/// the proof aggregation implementation and testing.
pub enum GeneralPCD<'a, G1: AffineCurve, G2: AffineCurve, D: Digest> {
    SimpleMarlin(SimpleMarlinPCD<'a, G1, D>),
    FinalDarlin(FinalDarlinPCD<'a, G1, G2, D>)
}

// Testing functions
impl<'a, G1, G2, D> GeneralPCD<'a, G1, G2, D>
    where
        G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
        G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
        D: Digest,
{
    pub fn randomize_usr_ins<R: RngCore>(
        &mut self,
        rng: &mut R
    )
    {
        match self {
            Self::SimpleMarlin(simple_marlin) => {
                // No sys ins (for now) for SimpleMarlin, so modify the usr_ins instead
                let ins_len = simple_marlin.usr_ins.len();
                simple_marlin.usr_ins = vec![G1::ScalarField::rand(rng); ins_len];
            },
            Self::FinalDarlin(final_darlin) => {
                let ins_len = final_darlin.usr_ins.len();
                final_darlin.usr_ins = vec![G1::ScalarField::rand(rng); ins_len];
            }
        }
    }

    pub fn randomize_sys_ins<R: RngCore>(
        &mut self,
        ck_g1: &DLogCommitterKey<G1>,
        ck_g2: &DLogCommitterKey<G2>,
        rng: &mut R
    )
    {
        match self {
            Self::SimpleMarlin(simple_marlin) => {
                // No sys ins (for now) for SimpleMarlin, so modify the usr_ins instead
                let ins_len = simple_marlin.usr_ins.len();
                simple_marlin.usr_ins = vec![G1::ScalarField::rand(rng); ins_len];
            },
            Self::FinalDarlin(final_darlin) => {
                final_darlin.final_darlin_proof.deferred = FinalDarlinDeferredData::<G1, G2>::generate_random::<R, D>(
                    rng, ck_g1, ck_g2
                );
            }
        }
    }
}

impl<'a, G1, G2, D> PCD for GeneralPCD<'a, G1, G2, D>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
    D: Digest + 'a,
{
    type PCDAccumulator = DualDLogItemAccumulator<'a, G1, G2, D>;
    type PCDVerifierKey = FinalDarlinPCDVerifierKey<'a, G1, G2, D>;

    fn succinct_verify(
        &self,
        vk: &Self::PCDVerifierKey,
    ) ->  Result<<Self::PCDAccumulator as ItemAccumulator>::Item, PCDError>
    {
        match self {
            Self::SimpleMarlin(simple_marlin) => {
                // Works because a FinalDarlinVk is a MarlinVk
                let simple_marlin_vk = SimpleMarlinPCDVerifierKey (vk.final_darlin_vk, vk.dlog_vks.0);
                let acc = simple_marlin.succinct_verify(&simple_marlin_vk)?;
                Ok(DualDLogItem (vec![acc], vec![]))
            },
            Self::FinalDarlin(final_darlin) => {
                final_darlin.succinct_verify(vk)
            }
        }
    }
}