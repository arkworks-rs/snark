use crate::Error;
use rand::Rng;
use std::marker::PhantomData;

use super::{
    pedersen::{PedersenCommitment, PedersenParameters, PedersenRandomness, PedersenWindow},
    CommitmentScheme,
};
pub use crate::crh::injective_map::InjectiveMap;
use algebra::groups::Group;

#[cfg(feature = "r1cs")]
pub mod constraints;

pub struct PedersenCommCompressor<G: Group, I: InjectiveMap<G>, W: PedersenWindow> {
    _group:      PhantomData<G>,
    _compressor: PhantomData<I>,
    _comm:       PedersenCommitment<G, W>,
}

impl<G: Group, I: InjectiveMap<G>, W: PedersenWindow> CommitmentScheme
    for PedersenCommCompressor<G, I, W>
{
    type Output = I::Output;
    type Parameters = PedersenParameters<G>;
    type Randomness = PedersenRandomness<G>;

    fn setup<R: Rng>(rng: &mut R) -> Result<Self::Parameters, Error> {
        let time = start_timer!(|| format!("PedersenCompressor::Setup"));
        let params = PedersenCommitment::<G, W>::setup(rng);
        end_timer!(time);
        params
    }

    fn commit(
        parameters: &Self::Parameters,
        input: &[u8],
        randomness: &Self::Randomness,
    ) -> Result<Self::Output, Error> {
        let eval_time = start_timer!(|| "PedersenCompressor::Eval");
        let result = I::injective_map(&PedersenCommitment::<G, W>::commit(
            parameters, input, randomness,
        )?)?;
        end_timer!(eval_time);
        Ok(result)
    }
}
