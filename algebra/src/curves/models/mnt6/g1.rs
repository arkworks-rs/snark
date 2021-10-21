use crate::curves::models::mnt6::{MNT6Parameters, MNT6p};
use crate::curves::short_weierstrass_projective::{GroupAffine, GroupProjective};
use crate::{AffineCurve, Fp3, FromBytes, ToBytes};
use serde::{Deserialize, Serialize};
use std::io;
use std::io::{Read, Result as IoResult, Write};

pub type G1Affine<P> = GroupAffine<<P as MNT6Parameters>::G1Parameters>;
pub type G1Projective<P> = GroupProjective<<P as MNT6Parameters>::G1Parameters>;

#[derive(Derivative)]
#[derivative(
    Copy(bound = "P: MNT6Parameters"),
    Clone(bound = "P: MNT6Parameters"),
    Debug(bound = "P: MNT6Parameters"),
    PartialEq(bound = "P: MNT6Parameters"),
    Eq(bound = "P: MNT6Parameters")
)]
#[derive(Serialize, Deserialize)]
#[serde(bound(serialize = "P: MNT6Parameters"))]
#[serde(bound(deserialize = "P: MNT6Parameters"))]
pub struct G1Prepared<P: MNT6Parameters> {
    pub p: G1Affine<P>,
    pub py_twist_squared: Fp3<P::Fp3Params>,
}

impl<P: MNT6Parameters> ToBytes for G1Prepared<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.p.write(&mut writer)?;
        self.py_twist_squared.write(&mut writer)?;
        Ok(())
    }
}

impl<P: MNT6Parameters> FromBytes for G1Prepared<P> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let p = G1Affine::<P>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let py_twist_squared = Fp3::<P::Fp3Params>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Ok(G1Prepared {
            p,
            py_twist_squared,
        })
    }
}

impl<P: MNT6Parameters> From<G1Affine<P>> for G1Prepared<P> {
    fn from(other: G1Affine<P>) -> Self {
        MNT6p::<P>::ate_precompute_g1(&other)
    }
}

impl<P: MNT6Parameters> Default for G1Prepared<P> {
    fn default() -> Self {
        Self::from(G1Affine::<P>::prime_subgroup_generator())
    }
}
