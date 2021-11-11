//! Trait for general (public, or "atomic") accumulation schemes [BCMS20](https://eprint.iacr.org/2020/499).
//! Comes with the aggregation/verification of "items", i.e. some data structure typically satisfying a
//! non-efficient predicate).  
//! The trait applies to mixed type accumulators as described in our Darlin Proof Tree document:
//! There, a (full) accumulator is a composite structure of dlog and inner sumcheck ("single") accumulators,
//! from both groups of the EC cycle (the "current", and the "collected" ones).
//! Although within recursion we do not separate accumulation strategy from the SNARK on protocol level,
//! we nevertheless serve this functionality for post processing outside the PCD.
use algebra::{serialize::*, AffineCurve};
use poly_commit::ipa_pc::Commitment;
use poly_commit::{ipa_pc::Proof, Error};
use rand::RngCore;

pub mod dlog;

/// General struct of an aggregation proof. Typically, such proof stems from an
/// interactive oracle protocol (IOP) and a polynomial commitment scheme.
#[derive(Clone, Default, Debug, Eq, PartialEq, CanonicalSerialize, CanonicalDeserialize)]
pub struct AccumulationProof<G: AffineCurve> {
    /// Commitments to the polynomials produced by the prover.
    pub commitments: Vec<Vec<Commitment<G>>>,
    /// Evaluations of these polynomials.
    pub evaluations: Vec<G::ScalarField>,
    /// An evaluation proof from the polynomial commitment.
    pub pc_proof: Proof<G>,
}

/// The `ItemAccumulator` trait comes with the essential functions for proving
/// and verifying aggregation, as well as checking ("deciding") if an item
/// satisfies the predicate.
/// It applies to mixed type accumulators as described in our [Darlin Proof Tree doc](TODO: add link):
/// There, a (full) accumulator is a composite structure of dlog and inner
/// sumcheck ("single") accumulators, from both groups of the EC cycle (the
/// "current", and the "collected" ones). Although within recursion we do
/// not separate accumulation strategy from the SNARK on protocol level,
/// we nevertheless serve this functionality for post processing outside the PCD.
pub trait ItemAccumulator {
    type AccumulatorProverKey;
    type AccumulatorVerifierKey;
    type AccumulationProof;
    type Item;

    /// Decide whether an/the public accumulator/s are correct,
    /// i.e. whether they satisfy the non-efficient predicate.
    /// Typically involves non-succinct MSMs.
    fn check_items<R: RngCore>(
        vk: &Self::AccumulatorVerifierKey,
        accumulators: &[Self::Item],
        rng: &mut R,
    ) -> Result<bool, Error>;

    /// Amortization strategy for items as a separate argument.  
    /// Returns the new/"updated" item and a non-interactive
    /// proof of its correct aggregation.
    fn accumulate_items(
        ck: &Self::AccumulatorProverKey,
        accumulators: Vec<Self::Item>,
    ) -> Result<(Self::Item, Self::AccumulationProof), Error>;

    /// Fully verifies a proof produced by accumulate_items() given the accumulators.
    /// Depending on the PC it may involve a non-succinct MSM.
    fn verify_accumulated_items<R: RngCore>(
        current_accumulator: &Self::Item,
        vk: &Self::AccumulatorVerifierKey,
        previous_accumulators: Vec<Self::Item>,
        proof: &Self::AccumulationProof,
        rng: &mut R,
    ) -> Result<bool, Error>;
}
