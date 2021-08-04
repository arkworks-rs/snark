//! An implementation of the [Groth][Groth16] zkSNARK.
//! [Groth16]: https://eprint.iacr.org/2016/260.pdf
use algebra::{Field, serialize::*, bytes::{
    ToBytes, FromBytes,
}, PairingEngine, FromBytesChecked, SemanticallyValid};
use r1cs_core::{SynthesisError, Index, LinearCombination};
use std::io::{self, Read, Result as IoResult, Write};
use byteorder::{BigEndian, ReadBytesExt, WriteBytesExt};
use serde::{Serialize, Deserialize};

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
#[derive(Clone, Debug, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
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
    /// Doesn't perform group membership check for deserialized points
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

impl<E: PairingEngine> SemanticallyValid for Proof<E> {
    #[inline]
    fn is_valid(&self) -> bool {
        self.a.is_valid() &&
        self.b.is_valid() &&
        self.c.is_valid()
    }
}

impl<E: PairingEngine> FromBytesChecked for Proof<E> {
    #[inline]
    fn read_checked<R: Read>(mut reader: R) -> IoResult<Self> {
        let a = E::G1Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point A: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point A: point at infinity")); }
                Ok(p)
            })?;

        let b = E::G2Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point B: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point B: point at infinity")); }
                Ok(p)
            })?;

        let c = E::G1Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point C: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point C: point at infinity")); }
                Ok(p)
            })?;

        Ok(Proof { a, b, c })
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

fn read_affine_vec_checked<G: AffineCurve, R: Read>(len: usize, zero_check: bool, mut reader: R) -> IoResult<Vec<G>> {
    let mut v = vec![];
    for i in 0..len {
        let g = G::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point {}: {}", i, e)))
            .and_then(|p| {
                if zero_check && p.is_zero()
                {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid point {}: point at infinity", i),
                    ));
                }
                Ok(p)
            })?;
        v.push(g);
    }
    Ok(v)
}

fn read_affine_vec<G: AffineCurve, R: Read>(len: usize, mut reader: R) -> IoResult<Vec<G>> {
    let mut v = vec![];
    for i in 0..len {
        let g = G::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point {}: {}", i, e)))?;
        v.push(g);
    }
    Ok(v)
}

/// A verification key in the Groth16 SNARK.
#[derive(Clone, Debug, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerifyingKey<E: PairingEngine> {
    pub alpha_g1_beta_g2:   E::Fqk,
    pub gamma_g2:           E::G2Affine,
    pub delta_g2:           E::G2Affine,
    pub gamma_abc_g1:       Vec<E::G1Affine>,
}

impl<E: PairingEngine> VerifyingKey<E> {
    fn check_gamma_abc_g1_points(gamma_abc_g1: &[E::G1Affine]) -> IoResult<()> {
        use std::ops::Neg;

        for (i, &p1) in gamma_abc_g1.iter().enumerate() {
            if p1.is_zero() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid gamma_abc_g1[{}]: point at infinity", i),
                ));
            }
            for (j, &p2) in gamma_abc_g1.iter().skip(i + 1).enumerate() {
                if p1 == p2 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("duplicate points: gamma_abc_g1[{}] = gamma_abc_g1[{}]", i, j),
                    ));
                }
                if p1 == p2.neg() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("inverse points: gamma_abc_g1[{}] = -gamma_abc_g1[{}]", i, j),
                    ));
                }
            }
        }
        Ok(())
    }
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

impl<E: PairingEngine> SemanticallyValid for VerifyingKey<E> {
    #[inline]
    fn is_valid(&self) -> bool {
        self.alpha_g1_beta_g2.is_valid() && !self.alpha_g1_beta_g2.is_zero() &&
        self.gamma_g2.is_valid() && !self.gamma_g2.is_zero() &&
        self.delta_g2.is_valid() && !self.delta_g2.is_zero() &&
        self.gamma_abc_g1.iter().filter(|&p| !p.is_valid() || p.is_zero())
            .collect::<Vec<_>>().is_empty() &&
        Self::check_gamma_abc_g1_points(self.gamma_abc_g1.as_slice()).is_ok()
    }
}

impl<E: PairingEngine> FromBytesChecked for VerifyingKey<E> {
    #[inline]
    fn read_checked<R: Read>(mut reader: R) -> IoResult<Self> {
        let alpha_g1_beta_g2 = E::Fqk::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid alpha_g1_beta_g2: {}", e)))
            .and_then(|f| {
                if f.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid alpha_g1_beta_g2: zero")); }
                Ok(f)
            })?;

        let gamma_g2 = E::G2Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point gamma_g2: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point gamma_g2: point at infinity")); }
                Ok(p)
            })?;

        let delta_g2 = E::G2Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point delta_g2: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point delta_g2: point at infinity")); }
                Ok(p)
            })?;

        let ic_len = reader.read_u32::<BigEndian>()? as usize;
        let gamma_abc_g1 = read_affine_vec_checked::<E::G1Affine, _>(ic_len, true, &mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid gamma_abc_g1: {}", e)))?;
        Self::check_gamma_abc_g1_points(gamma_abc_g1.as_slice())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid gamma_abc_g1: {}", e)))?;

        Ok(VerifyingKey { alpha_g1_beta_g2, gamma_g2, delta_g2, gamma_abc_g1 })
    }
}

impl<E: PairingEngine> FromBytes for VerifyingKey<E> {
    /// Doesn't perform group membership check for deserialized points
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {
        let alpha_g1_beta_g2 = E::Fqk::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let gamma_g2 = E::G2Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let delta_g2 = E::G2Affine::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let ic_len = reader.read_u32::<BigEndian>()? as usize;
        let gamma_abc_g1 = read_affine_vec::<E::G1Affine, _>(ic_len, &mut reader)?;

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

pub(crate) fn push_constraints<F: Field>(
    l: LinearCombination<F>,
    constraints: &mut [Vec<(F, Index)>],
    this_constraint: usize,
) {
    for (var, coeff) in l.as_ref() {
        match var.get_unchecked() {
            Index::Input(i) => constraints[this_constraint].push((*coeff, Index::Input(i))),
            Index::Aux(i) => constraints[this_constraint].push((*coeff, Index::Aux(i))),
        }
    }
}

/// Full public (prover and verifier) parameters for the Groth16 zkSNARK.
#[derive(Clone, Debug, Serialize, Deserialize, CanonicalSerialize, CanonicalDeserialize)]
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

impl<E: PairingEngine> SemanticallyValid for Parameters<E> {
    #[inline]
    fn is_valid(&self) -> bool {
        self.vk.is_valid() &&
        self.alpha_g1.is_valid() && !self.alpha_g1.is_zero() &&
        self.beta_g1.is_valid() && !self.beta_g1.is_zero() &&
        self.beta_g2.is_valid() && !self.beta_g2.is_zero() &&
        self.delta_g1.is_valid() && !self.delta_g1.is_zero() &&
        self.delta_g2.is_valid() && !self.delta_g2.is_zero() &&
        self.a_query.iter().filter(|&p| !p.is_valid())
            .collect::<Vec<_>>().is_empty() &&
        self.b_g1_query.iter().filter(|&p| !p.is_valid())
            .collect::<Vec<_>>().is_empty() &&
        self.b_g2_query.iter().filter(|&p| !p.is_valid())
            .collect::<Vec<_>>().is_empty() &&
        self.h_query.iter().filter(|&p| !p.is_valid() || p.is_zero())
            .collect::<Vec<_>>().is_empty() &&
        self.l_query.iter().filter(|&p| !p.is_valid() || p.is_zero())
            .collect::<Vec<_>>().is_empty()
    }
}

impl<E: PairingEngine> FromBytesChecked for Parameters<E> {
    #[inline]
    fn read_checked<R: Read>(mut reader: R) -> IoResult<Self> {
        let vk = VerifyingKey::<E>::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let alpha_g1 = E::G1Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point alpha_g1: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point alpha_g1: point at infinity")); }
                Ok(p)
            })?;

        let beta_g1 = E::G1Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point beta_g1: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point beta_g1: point at infinity")); }
                Ok(p)
            })?;

        let beta_g2 = E::G2Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point beta_g2: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point beta_g2: point at infinity")); }
                Ok(p)
            })?;

        let delta_g1 = E::G1Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point delta_g1: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point delta_g1: point at infinity")); }
                Ok(p)
            })?;

        let delta_g2 = E::G2Affine::read_checked(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid point delta_g2: {}", e)))
            .and_then(|p| {
                if p.is_zero() { return Err(io::Error::new(io::ErrorKind::InvalidData, "invalid point delta_g2: point at infinity")); }
                Ok(p)
            })?;

        // NOTE: a, b_g1 and b_g2 also contains polynomials that evaluate to zero, therefore the
        // zero check is disabled.
        // TODO: Exclude the points above from the generation procedure
        let a_len = reader.read_u32::<BigEndian>()? as usize;
        let a_query = read_affine_vec_checked::<E::G1Affine, _>(a_len, false, &mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid a_query: {}", e)))?;

        let b_g1_len = reader.read_u32::<BigEndian>()? as usize;
        let b_g1_query = read_affine_vec_checked::<E::G1Affine, _>(b_g1_len, false, &mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid b_g1_query: {}", e)))?;

        let b_g2_len = reader.read_u32::<BigEndian>()? as usize;
        let b_g2_query = read_affine_vec_checked::<E::G2Affine, _>(b_g2_len, false, &mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid b_g2_query: {}", e)))?;

        let h_len = reader.read_u32::<BigEndian>()? as usize;
        let h_query = read_affine_vec_checked::<E::G1Affine, _>(h_len, true, &mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid h_query: {}", e)))?;

        let l_len = reader.read_u32::<BigEndian>()? as usize;
        let l_query = read_affine_vec_checked::<E::G1Affine, _>(l_len, true, &mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("invalid l_query: {}", e)))?;

        Ok(Parameters { vk, alpha_g1, beta_g1, beta_g2, delta_g1, delta_g2, a_query, b_g1_query, b_g2_query, h_query, l_query })
    }
}

impl<E: PairingEngine> FromBytes for Parameters<E> {
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
        let a_query = read_affine_vec::<E::G1Affine, _>(a_len, &mut reader)?;
        let b_g1_len = reader.read_u32::<BigEndian>()? as usize;
        let b_g1_query = read_affine_vec::<E::G1Affine, _>(b_g1_len, &mut reader)?;
        let b_g2_len = reader.read_u32::<BigEndian>()? as usize;
        let b_g2_query = read_affine_vec::<E::G2Affine, _>(b_g2_len, &mut reader)?;
        let h_len = reader.read_u32::<BigEndian>()? as usize;
        let h_query = read_affine_vec::<E::G1Affine, _>(h_len, &mut reader)?;
        let l_len = reader.read_u32::<BigEndian>()? as usize;
        let l_query = read_affine_vec::<E::G1Affine, _>(l_len, &mut reader)?;
        Ok(Parameters{vk, alpha_g1, beta_g1, beta_g2, delta_g1, delta_g2, a_query, b_g1_query, b_g2_query, h_query, l_query})
    }
}

/// Preprocessed verification key parameters that enable faster verification
/// at the expense of larger size in memory.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PreparedVerifyingKey<E: PairingEngine> {
    pub alpha_g1_beta_g2: E::Fqk,
    pub gamma_g2_neg_pc:  E::G2Prepared,
    pub delta_g2_neg_pc:  E::G2Prepared,
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
            gamma_g2_neg_pc:  E::G2Prepared::default(),
            delta_g2_neg_pc:  E::G2Prepared::default(),
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

//TODO: Maybe implemet FromBytesChecked also for PreparedVk ?
impl<E: PairingEngine> FromBytes for PreparedVerifyingKey<E> {
    #[inline]
    fn read<R: Read>(mut reader: R) -> IoResult<Self> {

        let alpha_g1_beta_g2 = E::Fqk::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let gamma_g2_neg_pc = E::G2Prepared::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let delta_g2_neg_pc = E::G2Prepared::read(&mut reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

        let ic_len = reader.read_u32::<BigEndian>()? as usize;
        let gamma_abc_g1 = read_affine_vec::<E::G1Affine, _>(ic_len, &mut reader)?;

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