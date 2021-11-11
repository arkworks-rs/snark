use crate::{
    bytes::ToBytes,
    curves::{
        bls12::Bls12Parameters,
        short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        AffineCurve,
    },
    FromBytes,
};
use serde::{Deserialize, Serialize};
use std::io;
use std::io::{Read, Result as IoResult, Write};

pub type G1Affine<P> = GroupAffine<<P as Bls12Parameters>::G1Parameters>;
pub type G1Projective<P> = GroupProjective<<P as Bls12Parameters>::G1Parameters>;

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: Bls12Parameters"),
    Debug(bound = "P: Bls12Parameters"),
    PartialEq(bound = "P: Bls12Parameters"),
    Eq(bound = "P: Bls12Parameters")
)]
#[derive(Serialize, Deserialize)]
#[serde(bound(serialize = "P: Bls12Parameters"))]
#[serde(bound(deserialize = "P: Bls12Parameters"))]
#[serde(transparent)]
pub struct G1Prepared<P: Bls12Parameters>(pub G1Affine<P>);

impl<P: Bls12Parameters> G1Prepared<P> {
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<P: Bls12Parameters> From<G1Affine<P>> for G1Prepared<P> {
    fn from(other: G1Affine<P>) -> Self {
        G1Prepared(other)
    }
}

impl<P: Bls12Parameters> Default for G1Prepared<P> {
    fn default() -> Self {
        G1Prepared(G1Affine::<P>::prime_subgroup_generator())
    }
}

impl<P: Bls12Parameters> ToBytes for G1Prepared<P> {
    fn write<W: Write>(&self, writer: W) -> IoResult<()> {
        self.0.write(writer)
    }
}

impl<P: Bls12Parameters> FromBytes for G1Prepared<P> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let g1a = G1Affine::<P>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Ok(G1Prepared(g1a))
    }
}
