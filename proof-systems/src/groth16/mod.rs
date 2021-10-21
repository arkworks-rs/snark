//! An implementation of the [Groth][Groth16] zkSNARK.
//! [Groth16]: https://eprint.iacr.org/2016/260.pdf
use algebra::{bytes::{
    ToBytes, FromBytes,
}, PairingCurve, PairingEngine};
use r1cs_core::SynthesisError;
use std::io::{self, Read, Result as IoResult, Write};
use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};

/// Reduce an R1CS instance to a *Quadratic Arithmetic Program* instance.
pub mod r1cs_to_qap;

/// Generate public parameters for the Groth16 zkSNARK construction.
pub mod generator;

/// Create proofs for the Groth16 zkSNARK construction.
pub mod prover;

/// Verify proofs for the Groth16 zkSNARK construction.
pub mod verifier;

#[cfg(test)]
mod test;

pub use self::{generator::*, prover::*, verifier::*};

/// A proof in the Groth16 SNARK.
#[derive(Clone, Debug)]
pub struct Proof<E: PairingEngine> {
    pub a: E::G1Affine,
    pub b: E::G2Affine,
    pub c: E::G1Affine,
}

impl<E: PairingEngine> ToBytes for Proof<E> {
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> io::Result<()> {
        self.a.write(&mut writer)?;
        self.b.write(&mut writer)?;
        self.c.write(&mut writer)
    }
}

impl<E: PairingEngine> FromBytes for Proof<E> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let a = E::G1Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let b = E::G2Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let c = E::G1Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Ok(Proof{a, b, c})
    }
}

impl<E: PairingEngine> PartialEq for Proof<E> {
    fn eq(&self, other: &Self) -> bool {
        self.a == other.a && self.b == other.b && self.c == other.c
    }
}

impl<E: PairingEngine> Default for Proof<E> {
    fn default() -> Self {
        Self {
            a: E::G1Affine::default(),
            b: E::G2Affine::default(),
            c: E::G1Affine::default(),
        }
    }
}

use algebra::curves::AffineCurve;

fn read_affine_vec<G: AffineCurve, R: Read>(len: usize, check_for_zero: bool, mut reader: R) -> IoResult<Vec<G>> {
    let mut v = vec![];
    for _ in 0..len {
        let g = G::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            .and_then(|e| {
                if check_for_zero && e.is_zero() {
                    Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "point at infinity",
                    ))
                } else {
                    Ok(e)
                }
            })?;
        v.push(g);
    }
    Ok(v)
}

/// A verification key in the Groth16 SNARK.
#[derive(Clone, Debug)]
pub struct VerifyingKey<E: PairingEngine> {
    pub alpha_g1_beta_g2:   E::Fqk,
    pub gamma_g2:           E::G2Affine,
    pub delta_g2:           E::G2Affine,
    pub gamma_abc_g1:       Vec<E::G1Affine>,
}

impl<E: PairingEngine> ToBytes for VerifyingKey<E> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.alpha_g1_beta_g2.write(&mut writer)?;
        self.gamma_g2.write(&mut writer)?;
        self.delta_g2.write(&mut writer)?;
        writer.write_u32::<BigEndian>(self.gamma_abc_g1.len() as u32)?;
        for q in &self.gamma_abc_g1 {
            q.write(&mut writer)?;
        }
        Ok(())
    }
}

impl<E: PairingEngine> FromBytes for VerifyingKey<E> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let alpha_g1_beta_g2 = E::Fqk::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let gamma_g2 = E::G2Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let delta_g2 = E::G2Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let ic_len = reader.read_u32::<BigEndian>()? as usize;
        let gamma_abc_g1 = read_affine_vec::<E::G1Affine, _>(ic_len, true, &mut reader)?;

        Ok(VerifyingKey{alpha_g1_beta_g2, gamma_g2, delta_g2, gamma_abc_g1})
    }
}


impl<E: PairingEngine> Default for VerifyingKey<E> {
    fn default() -> Self {
        Self {
            alpha_g1_beta_g2:     E::Fqk::default(),
            gamma_g2:             E::G2Affine::default(),
            delta_g2:             E::G2Affine::default(),
            gamma_abc_g1:         Vec::new(),
        }
    }
}

impl<E: PairingEngine> PartialEq for VerifyingKey<E> {
    fn eq(&self, other: &Self) -> bool {
        self.alpha_g1_beta_g2 == other.alpha_g1_beta_g2
            && self.gamma_g2 == other.gamma_g2
            && self.delta_g2 == other.delta_g2
            && self.gamma_abc_g1 == other.gamma_abc_g1
    }
}


/// Full public (prover and verifier) parameters for the Groth16 zkSNARK.
#[derive(Clone, Debug)]
pub struct Parameters<E: PairingEngine> {
    pub vk:         VerifyingKey<E>,
    pub alpha_g1:   E::G1Affine,
    pub beta_g1:    E::G1Affine,
    pub beta_g2:    E::G2Affine,
    pub delta_g1:   E::G1Affine,
    pub delta_g2:   E::G2Affine,
    pub a_query:    Vec<E::G1Affine>,
    pub b_g1_query: Vec<E::G1Affine>,
    pub b_g2_query: Vec<E::G2Affine>,
    pub h_query:    Vec<E::G1Affine>,
    pub l_query:    Vec<E::G1Affine>,
}

impl<E: PairingEngine> PartialEq for Parameters<E> {
    fn eq(&self, other: &Self) -> bool {
        self.vk == other.vk
            && self.alpha_g1 == other.alpha_g1
            && self.beta_g1 == other.beta_g1
            && self.beta_g2 == other.beta_g2
            && self.delta_g1 == other.delta_g1
            && self.delta_g2 == other.delta_g2
            && self.a_query == other.a_query
            && self.b_g1_query == other.b_g1_query
            && self.b_g2_query == other.b_g2_query
            && self.h_query == other.h_query
            && self.l_query == other.l_query
    }
}

impl<E: PairingEngine> ToBytes for Parameters<E>{
    #[inline]
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.vk.write(&mut writer)?;
        self.alpha_g1.write(&mut writer)?;
        self.beta_g1.write(&mut writer)?;
        self.beta_g2.write(&mut writer)?;
        self.delta_g1.write(&mut writer)?;
        self.delta_g2.write(&mut writer)?;
        writer.write_u32::<BigEndian>(self.a_query.len() as u32)?;
        for a in self.a_query.clone() {a.write(&mut writer)?;}
        writer.write_u32::<BigEndian>(self.b_g1_query.len() as u32)?;
        for a in self.b_g1_query.clone() {a.write(&mut writer)?;}
        writer.write_u32::<BigEndian>(self.b_g2_query.len() as u32)?;
        for a in self.b_g2_query.clone() {a.write(&mut writer)?;}
        writer.write_u32::<BigEndian>(self.h_query.len() as u32)?;
        for a in self.h_query.clone() {a.write(&mut writer)?;}
        writer.write_u32::<BigEndian>(self.l_query.len() as u32)?;
        for a in self.l_query.clone() {a.write(&mut writer)?;}
        Ok(())
    }
}

impl<E: PairingEngine> FromBytes for Parameters<E>{
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let vk = VerifyingKey::<E>::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let alpha_g1 = E::G1Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let beta_g1 = E::G1Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let beta_g2 = E::G2Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let delta_g1 = E::G1Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let delta_g2 = E::G2Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let a_len = reader.read_u32::<BigEndian>()? as usize;
        let a_query = read_affine_vec::<E::G1Affine, _>(a_len, false, &mut reader)?;
        let b_g1_len = reader.read_u32::<BigEndian>()? as usize;
        let b_g1_query = read_affine_vec::<E::G1Affine, _>(b_g1_len, false, &mut reader)?;
        let b_g2_len = reader.read_u32::<BigEndian>()? as usize;
        let b_g2_query = read_affine_vec::<E::G2Affine, _>(b_g2_len, false, &mut reader)?;
        let h_len = reader.read_u32::<BigEndian>()? as usize;
        let h_query = read_affine_vec::<E::G1Affine, _>(h_len, false, &mut reader)?;
        let l_len = reader.read_u32::<BigEndian>()? as usize;
        let l_query = read_affine_vec::<E::G1Affine, _>(l_len, false, &mut reader)?;
        Ok(Parameters{vk, alpha_g1, beta_g1, beta_g2, delta_g1, delta_g2, a_query, b_g1_query, b_g2_query, h_query, l_query})

    }
}

/// Preprocessed verification key parameters that enable faster verification
/// at the expense of larger size in memory.
#[derive(Clone, Debug)]
pub struct PreparedVerifyingKey<E: PairingEngine> {
    pub alpha_g1_beta_g2: E::Fqk,
    pub gamma_g2_neg_pc:  <E::G2Affine as PairingCurve>::Prepared,
    pub delta_g2_neg_pc:  <E::G2Affine as PairingCurve>::Prepared,
    pub gamma_abc_g1:     Vec<E::G1Affine>,
}

impl<E: PairingEngine> From<VerifyingKey<E>> for PreparedVerifyingKey<E> {
    fn from(other: VerifyingKey<E>) -> Self {
        prepare_verifying_key(&other)
    }
}

impl<E: PairingEngine> Default for PreparedVerifyingKey<E> {
    fn default() -> Self {
        Self {
            alpha_g1_beta_g2: E::Fqk::default(),
            gamma_g2_neg_pc:  <E::G2Affine as PairingCurve>::Prepared::default(),
            delta_g2_neg_pc:  <E::G2Affine as PairingCurve>::Prepared::default(),
            gamma_abc_g1:     Vec::new(),
        }
    }
}

impl<E: PairingEngine> ToBytes for PreparedVerifyingKey<E> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.alpha_g1_beta_g2.write(&mut writer)?;
        self.gamma_g2_neg_pc.write(&mut writer)?;
        self.delta_g2_neg_pc.write(&mut writer)?;
        writer.write_u32::<BigEndian>(self.gamma_abc_g1.len() as u32)?;
        for q in &self.gamma_abc_g1 {
            q.write(&mut writer)?;
        }
        Ok(())
    }
}

impl<E: PairingEngine> FromBytes for PreparedVerifyingKey<E> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {

        let alpha_g1_beta_g2 = E::Fqk::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let gamma_g2_neg_pc = <E::G2Affine as PairingCurve>::Prepared::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let delta_g2_neg_pc = <E::G2Affine as PairingCurve>::Prepared::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let ic_len = reader.read_u32::<BigEndian>()? as usize;
        let gamma_abc_g1 = read_affine_vec::<E::G1Affine, _>(ic_len, true, &mut reader)?;

        Ok(PreparedVerifyingKey {
            alpha_g1_beta_g2,
            gamma_g2_neg_pc,
            delta_g2_neg_pc,
            gamma_abc_g1,
        })
    }
}

impl<E: PairingEngine> PartialEq for PreparedVerifyingKey<E> {
    fn eq(&self, other: &Self) -> bool {
        self.alpha_g1_beta_g2 == other.alpha_g1_beta_g2
            && self.gamma_g2_neg_pc == other.gamma_g2_neg_pc
            && self.delta_g2_neg_pc == other.delta_g2_neg_pc
            && self.gamma_abc_g1 == other.gamma_abc_g1
    }
}

impl<E: PairingEngine> Parameters<E> {
    pub fn get_vk(&self, _: usize) -> Result<VerifyingKey<E>, SynthesisError> {
        Ok(self.vk.clone())
    }

    pub fn get_a_query(
        &self,
        num_inputs: usize,
    ) -> Result<(&[E::G1Affine], &[E::G1Affine]), SynthesisError> {
        Ok((&self.a_query[1..num_inputs], &self.a_query[num_inputs..]))
    }

    pub fn get_b_g1_query(
        &self,
        num_inputs: usize,
    ) -> Result<(&[E::G1Affine], &[E::G1Affine]), SynthesisError> {
        Ok((
            &self.b_g1_query[1..num_inputs],
            &self.b_g1_query[num_inputs..],
        ))
    }

    pub fn get_b_g2_query(
        &self,
        num_inputs: usize,
    ) -> Result<(&[E::G2Affine], &[E::G2Affine]), SynthesisError> {
        Ok((
            &self.b_g2_query[1..num_inputs],
            &self.b_g2_query[num_inputs..],
        ))
    }

    pub fn get_h_query(
        &self,
        num_inputs: usize,
    ) -> Result<(&[E::G1Affine], &[E::G1Affine]), SynthesisError> {
        Ok((&self.h_query[0..num_inputs], &self.h_query[num_inputs..]))
    }

    pub fn get_a_query_full(&self) -> Result<&[E::G1Affine], SynthesisError> {
        Ok(&self.a_query)
    }

    pub fn get_b_g1_query_full(&self) -> Result<&[E::G1Affine], SynthesisError> {
        Ok(&self.b_g1_query)
    }

    pub fn get_b_g2_query_full(&self) -> Result<&[E::G2Affine], SynthesisError> {
        Ok(&self.b_g2_query)
    }

    pub fn get_h_query_full(&self) -> Result<&[E::G1Affine], SynthesisError> {
        Ok(&self.h_query)
    }

    pub fn get_l_query_full(&self) -> Result<&[E::G1Affine], SynthesisError> {
        Ok(&self.l_query)
    }
}