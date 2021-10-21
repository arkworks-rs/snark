use crate::curves::models::mnt6::{MNT6Parameters, MNT6p};
use crate::curves::short_weierstrass_projective::{GroupAffine, GroupProjective};
use crate::{AffineCurve, Fp3, FromBytes, ToBytes};
use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use serde::{Deserialize, Serialize};
use std::io;
use std::io::{Read, Result as IoResult, Write};

pub type G2Affine<P> = GroupAffine<<P as MNT6Parameters>::G2Parameters>;
pub type G2Projective<P> = GroupProjective<<P as MNT6Parameters>::G2Parameters>;

#[derive(Derivative)]
#[derivative(
    Copy(bound = "P: MNT6Parameters"),
    Clone(bound = "P: MNT6Parameters"),
    Debug(bound = "P: MNT6Parameters"),
    PartialEq(bound = "P: MNT6Parameters"),
    Eq(bound = "P: MNT6Parameters")
)]
#[derive(Serialize, Deserialize)]
pub struct G2PreparedCoefficients<P: MNT6Parameters> {
    pub r_y: Fp3<P::Fp3Params>,
    pub gamma: Fp3<P::Fp3Params>,
    pub gamma_x: Fp3<P::Fp3Params>,
}

impl<P: MNT6Parameters> ToBytes for G2PreparedCoefficients<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.r_y.write(&mut writer)?;
        self.gamma.write(&mut writer)?;
        self.gamma_x.write(&mut writer)?;
        Ok(())
    }
}

impl<P: MNT6Parameters> FromBytes for G2PreparedCoefficients<P> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let r_y = Fp3::<P::Fp3Params>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let gamma = Fp3::<P::Fp3Params>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let gamma_x = Fp3::<P::Fp3Params>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Ok(G2PreparedCoefficients {
            r_y,
            gamma,
            gamma_x,
        })
    }
}

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: MNT6Parameters"),
    Debug(bound = "P: MNT6Parameters"),
    PartialEq(bound = "P: MNT6Parameters"),
    Eq(bound = "P: MNT6Parameters")
)]
#[derive(Serialize, Deserialize)]
#[serde(bound(serialize = "P: MNT6Parameters"))]
#[serde(bound(deserialize = "P: MNT6Parameters"))]
pub struct G2Prepared<P: MNT6Parameters> {
    pub q: G2Affine<P>,
    pub coeffs: Vec<G2PreparedCoefficients<P>>,
}

impl<P: MNT6Parameters> ToBytes for G2Prepared<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.q.write(&mut writer)?;
        writer.write_u32::<BigEndian>(self.coeffs.len() as u32)?;
        for c in &self.coeffs {
            c.write(&mut writer)?;
        }
        Ok(())
    }
}

impl<P: MNT6Parameters> FromBytes for G2Prepared<P> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let q = G2Affine::<P>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let coeffs_len = reader.read_u32::<BigEndian>()? as usize;
        let mut coeffs = vec![];

        for _ in 0..coeffs_len {
            let c = G2PreparedCoefficients::<P>::read(&mut reader)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            coeffs.push(c);
        }
        Ok(G2Prepared { q, coeffs })
    }
}

impl<P: MNT6Parameters> From<G2Affine<P>> for G2Prepared<P> {
    fn from(point: G2Affine<P>) -> Self {
        MNT6p::<P>::ate_precompute_g2(&point).unwrap()
    }
}

impl<P: MNT6Parameters> Default for G2Prepared<P> {
    fn default() -> Self {
        Self::from(G2Affine::<P>::prime_subgroup_generator())
    }
}
