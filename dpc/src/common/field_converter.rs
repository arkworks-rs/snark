use algebra::{utils::*, Group, PairingEngine};
use crate::Error;

use crate::{
    crypto_primitives::{
        commitment::pedersen::PedersenParameters as PCParameters,
        crh::{pedersen::PedersenParameters as PHParameters, FixedLengthCRH},
        signature::schnorr::SchnorrSigParameters,
    },
    ledger::ideal_ledger::Digest,
};
use digest::Digest as HashDigest;

impl<E: PairingEngine, G: Group + ToEngineFr<E>> ToEngineFr<E> for PCParameters<G> {
    #[inline]
    fn to_engine_fr(&self) -> Result<Vec<E::Fr>, Error> {
        Ok(Vec::new())
    }
}

impl<E: PairingEngine, G: Group + ToEngineFr<E>> ToEngineFr<E> for PHParameters<G> {
    #[inline]
    fn to_engine_fr(&self) -> Result<Vec<E::Fr>, Error> {
        Ok(Vec::new())
    }
}

impl<E: PairingEngine, G: Group + ToEngineFr<E>, D: HashDigest> ToEngineFr<E>
    for SchnorrSigParameters<G, D>
{
    #[inline]
    fn to_engine_fr(&self) -> Result<Vec<E::Fr>, Error> {
        self.generator.to_engine_fr()
    }
}

impl<E: PairingEngine, H: FixedLengthCRH> ToEngineFr<E> for Digest<H>
where
    H::Output: ToEngineFr<E>,
{
    #[inline]
    fn to_engine_fr(&self) -> Result<Vec<E::Fr>, Error> {
        self.0.to_engine_fr()
    }
}
