use crate::{
    bytes::{FromBytes, ToBytes},
    curves::{
        bn::BnParameters,
        short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        AffineCurve,
    },
};

use serde::{Deserialize, Serialize};
use std::io::{Error, ErrorKind, Read, Result as IoResult, Write};

pub type G1Affine<P> = GroupAffine<<P as BnParameters>::G1Parameters>;
pub type G1Projective<P> = GroupProjective<<P as BnParameters>::G1Parameters>;

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: BnParameters"),
    Debug(bound = "P: BnParameters"),
    PartialEq(bound = "P: BnParameters"),
    Eq(bound = "P: BnParameters")
)]
#[derive(Serialize, Deserialize)]
#[serde(bound(serialize = "P: BnParameters"))]
#[serde(bound(deserialize = "P: BnParameters"))]
#[serde(transparent)]
pub struct G1Prepared<P: BnParameters>(pub G1Affine<P>);

impl<P: BnParameters> G1Prepared<P> {
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
}

impl<P: BnParameters> From<G1Affine<P>> for G1Prepared<P> {
    fn from(other: G1Affine<P>) -> Self {
        G1Prepared(other)
    }
}

impl<P: BnParameters> Default for G1Prepared<P> {
    fn default() -> Self {
        G1Prepared(G1Affine::<P>::prime_subgroup_generator())
    }
}

impl<P: BnParameters> ToBytes for G1Prepared<P> {
    fn write<W: Write>(&self, writer: W) -> IoResult<()> {
        self.0.write(writer)
    }
}

impl<P: BnParameters> FromBytes for G1Prepared<P> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let g1a =
            G1Affine::<P>::read(&mut reader).map_err(|e| Error::new(ErrorKind::InvalidData, e))?;
        Ok(G1Prepared(g1a))
    }
}
