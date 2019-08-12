use algebra::{utils::*, Group, Field};
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

impl<ConstraintF: Field, G: Group + ToConstraintField<ConstraintF>> ToConstraintField<ConstraintF> for PCParameters<G> {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        Ok(Vec::new())
    }
}

impl<ConstraintF: Field, G: Group + ToConstraintField<ConstraintF>> ToConstraintField<ConstraintF> for PHParameters<G> {
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        Ok(Vec::new())
    }
}

impl<ConstraintF: Field, G: Group + ToConstraintField<ConstraintF>, D: HashDigest> ToConstraintField<ConstraintF>
    for SchnorrSigParameters<G, D>
{
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        self.generator.to_field_elements()
    }
}

impl<ConstraintF: Field, H: FixedLengthCRH> ToConstraintField<ConstraintF> for Digest<H>
where
    H::Output: ToConstraintField<ConstraintF>,
{
    #[inline]
    fn to_field_elements(&self) -> Result<Vec<ConstraintF>, Error> {
        self.0.to_field_elements()
    }
}
