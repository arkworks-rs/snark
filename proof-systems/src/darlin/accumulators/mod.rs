use algebra::{ToBytes, AffineCurve};
use rand::RngCore;
use poly_commit::{
    ipa_pc::Proof,
    Error
};
use digest::Digest;
use poly_commit::ipa_pc::Commitment;

pub mod dlog;

/// General struct of a proof for an amortization strategy not
/// based on a separate definition of an information-theoretic protocol.
#[derive(Default)]
pub struct AccumulationProof<G: AffineCurve> {
    /// Commitments to the polynomials produced by the prover.
    pub commitments: Vec<Vec<Commitment<G>>>,
    /// Evaluations of these polynomials.
    pub evaluations: Vec<G::ScalarField>,
    /// An evaluation proof from the polynomial commitment.
    pub pc_proof: Proof<G>,
}

/// This trait embraces our accumulator notion from the Darlin Proof Tree doc:
/// There, a (full) accumulator is a composite structure of dlog and inner sumcheck "single"
/// accumulators from the two groups of the cycle (the "current", and the "collected" ones).
/// It is considered as part of the PCDProof. Although within our PCD we do not separate the
/// accumulation strategy from the proving system, we nevertheless serve this functionality
/// for post processing outside the PCD.
pub trait Accumulator: Sized + ToBytes {
    type AccumulatorProverKey;
    type AccumulatorVerifierKey;
    type AccumulationProof;

    /// Decide whether the public accumulators are correct.
    /// Typically involves a non-succinct MSM.
    fn check_accumulators<R: RngCore, D: Digest>(
        vk: &Self::AccumulatorVerifierKey,
        accumulators: &[Self],
        rng: &mut R,
    ) -> Result<bool, Error>;

    /// Amortization strategy for Accumulator as a separate protocol.
    /// Return the new "updated" accumulator and a non-interactive
    /// proof of its correct derivation from the given accumulators
    /// to be aggregated.
    fn accumulate<R: RngCore, D: Digest>(
        ck: &Self::AccumulatorProverKey,
        accumulators: Vec<Self>,
        rng: &mut R,
    ) -> Result<(Self, Self::AccumulationProof), Error>;

    /// Fully verifies a proof produced by accumulate() given the accumulator.
    /// Depending on the PC it may involve a non-succinct MSM.
    fn verify_accumulate<R: RngCore, D: Digest>(
        &self,
        vk: &Self::AccumulatorVerifierKey,
        previous_accumulators: Vec<Self>,
        proof: &Self::AccumulationProof,
        rng: &mut R,
    ) -> Result<bool, Error>;
}