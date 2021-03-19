use algebra::{AffineCurve, ToConstraintField, UniformRand};
use r1cs_core::ConstraintSynthesizer;
use marlin::{ProverKey as MarlinProverKey, VerifierKey as MarlinVerifierKey, Error as MarlinError, AHPForR1CS};
use poly_commit::{
    ipa_pc::{
        InnerProductArgPC, UniversalParams,
        CommitterKey as DLogCommitterKey, VerifierKey as DLogVerifierKey,
    },
    PolynomialCommitment, Error as PCError
};
use crate::darlin::accumulators::Accumulator;
use rand::RngCore;
use digest::Digest;
use crate::darlin::pcd::final_darlin::{FinalDarlinPCD, FinalDarlinPCDVerifierKey, FinalDarlinDeferredData};
use crate::darlin::pcd::simple_marlin::{SimpleMarlinPCD, SimpleMarlinPCDVerifierKey};
use crate::darlin::accumulators::dlog::DualDLogAccumulator;

pub mod simple_marlin;
pub mod final_darlin;

pub struct PCDParameters {
    pub segment_size: usize
}

impl PCDParameters {
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

    pub fn circuit_specific_setup<G: AffineCurve, C: ConstraintSynthesizer<G::ScalarField>, D: Digest>(
        &self,
        circuit: C,
        ck: &DLogCommitterKey<G>,
    ) -> Result<
        (
            MarlinProverKey<G::ScalarField, InnerProductArgPC<G, D>>,
            MarlinVerifierKey<G::ScalarField, InnerProductArgPC<G, D>>
        ), MarlinError<PCError>
    >
    {
        let index = AHPForR1CS::index(circuit)?;

        let (index_comms, index_comm_rands): (_, _) =
            InnerProductArgPC::<G, D>::commit(ck, index.iter(), None).map_err(MarlinError::from_pc_err)?;

        let index_comms = index_comms
            .into_iter()
            .map(|c| c.commitment().clone())
            .collect();

        let index_vk = MarlinVerifierKey {
            index_info: index.index_info,
            index_comms,
        };

        let index_pk = MarlinProverKey {
            index,
            index_comm_rands,
            index_vk: index_vk.clone(),
        };

        Ok((index_pk, index_vk))
    }
}

/// this trait expresses the functions for proof carrying data, in which the PCD is assumed
/// to be a set of data consisting of a statement, some deferred elements and a proof.
pub trait PCD<'a>: Sized + Send + Sync {
    type PCDAccumulator: Accumulator<'a>;
    type PCDVerifierKey: AsRef<<Self::PCDAccumulator as Accumulator<'a>>::AccumulatorVerifierKey>;

    fn succinct_verify(
        &self,
        vk:         &Self::PCDVerifierKey,
    ) -> Result<Self::PCDAccumulator, PCError>;

    fn hard_verify<R: RngCore, D: Digest>(
        &self,
        acc:    Self::PCDAccumulator,
        vk:     &Self::PCDVerifierKey,
        rng:    &mut R,
    ) -> Result<bool, PCError>
    { <Self::PCDAccumulator as Accumulator>::check_accumulators::<R, D>(vk.as_ref(), &[acc], rng) }

    fn verify<R: RngCore, D: Digest>(
        &self,
        vk:         &Self::PCDVerifierKey,
        rng:        &mut R,
    ) -> Result<bool, PCError>
    {
        let acc = self.succinct_verify(vk)?;
        self.hard_verify::<R, D>(acc, vk, rng)
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = ""))]
pub enum GeneralPCD<G1: AffineCurve, G2: AffineCurve, D: Digest> {
    SimpleMarlin(SimpleMarlinPCD<G1, D>),
    FinalDarlin(FinalDarlinPCD<G1, G2, D>)
}

// Testing functions
impl<G1, G2, D> GeneralPCD<G1, G2, D>
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
                final_darlin.deferred = FinalDarlinDeferredData::<G1, G2>::generate_random::<R, D>(
                    rng, ck_g1, ck_g2
                );
            }
        }
    }
}

impl<'a, G1, G2, D> PCD<'a> for GeneralPCD<G1, G2, D>
where
    G1: AffineCurve<BaseField = <G2 as AffineCurve>::ScalarField> + ToConstraintField<<G2 as AffineCurve>::ScalarField>,
    G2: AffineCurve<BaseField = <G1 as AffineCurve>::ScalarField> + ToConstraintField<<G1 as AffineCurve>::ScalarField>,
    D: Digest + 'a,
{
    type PCDAccumulator = DualDLogAccumulator<G1, G2>;
    type PCDVerifierKey = FinalDarlinPCDVerifierKey<'a, G1, G2, D>;

    fn succinct_verify(
        &self,
        vk: &Self::PCDVerifierKey,
    ) ->  Result<Self::PCDAccumulator, PCError>
    {
        match self {
            Self::SimpleMarlin(simple_marlin) => {
                let simple_marlin_vk = SimpleMarlinPCDVerifierKey (vk.marlin_vk, vk.dlog_vks.0);
                let acc = simple_marlin.succinct_verify(&simple_marlin_vk)?;
                Ok(DualDLogAccumulator (vec![acc], vec![]))
            },
            Self::FinalDarlin(final_darlin) => {
                final_darlin.succinct_verify(vk)
            }
        }
    }
}