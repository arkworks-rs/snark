use crate::curves::models::mnt4::{MNT4Parameters, MNT4p};
use crate::curves::short_weierstrass_projective::{GroupAffine, GroupProjective};
use crate::{AffineCurve, Fp2, FromBytes, ToBytes};
use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use serde::{Deserialize, Serialize};
use std::io;
use std::io::{Read, Result as IoResult, Write};

pub type G2Affine<P> = GroupAffine<<P as MNT4Parameters>::G2Parameters>;
pub type G2Projective<P> = GroupProjective<<P as MNT4Parameters>::G2Parameters>;

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: MNT4Parameters"),
    Debug(bound = "P: MNT4Parameters"),
    PartialEq(bound = "P: MNT4Parameters"),
    Eq(bound = "P: MNT4Parameters")
)]
#[derive(Serialize, Deserialize)]
pub struct G2PreparedCoefficients<P: MNT4Parameters> {
    pub r_y: Fp2<P::Fp2Params>,
    pub gamma: Fp2<P::Fp2Params>,
    pub gamma_x: Fp2<P::Fp2Params>,
}

impl<P: MNT4Parameters> ToBytes for G2PreparedCoefficients<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.r_y.write(&mut writer)?;
        self.gamma.write(&mut writer)?;
        self.gamma_x.write(&mut writer)?;
        Ok(())
    }
}

impl<P: MNT4Parameters> FromBytes for G2PreparedCoefficients<P> {
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let r_y = Fp2::<P::Fp2Params>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let gamma = Fp2::<P::Fp2Params>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let gamma_x = Fp2::<P::Fp2Params>::read(&mut reader)
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
    Clone(bound = "P: MNT4Parameters"),
    Debug(bound = "P: MNT4Parameters"),
    PartialEq(bound = "P: MNT4Parameters"),
    Eq(bound = "P: MNT4Parameters")
)]
#[derive(Serialize, Deserialize)]
#[serde(bound(serialize = "P: MNT4Parameters"))]
#[serde(bound(deserialize = "P: MNT4Parameters"))]
pub struct G2Prepared<P: MNT4Parameters> {
    pub q: G2Affine<P>,
    pub coeffs: Vec<G2PreparedCoefficients<P>>,
}

impl<P: MNT4Parameters> ToBytes for G2Prepared<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.q.write(&mut writer)?;
        writer.write_u32::<BigEndian>(self.coeffs.len() as u32)?;
        for c in &self.coeffs {
            c.write(&mut writer)?;
        }
        Ok(())
    }
}

impl<P: MNT4Parameters> FromBytes for G2Prepared<P> {
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

impl<P: MNT4Parameters> From<G2Affine<P>> for G2Prepared<P> {
    fn from(point: G2Affine<P>) -> Self {
        MNT4p::<P>::ate_precompute_g2(&point)
    }
}

impl<P: MNT4Parameters> Default for G2Prepared<P> {
    fn default() -> Self {
        Self::from(G2Affine::<P>::prime_subgroup_generator())
    }
}
