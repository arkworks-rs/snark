//! An implementation of the [`Groth16`] zkSNARK.
//!
//! [`Groth16`]: https://eprint.iacr.org/2016/260.pdf
#![cfg_attr(not(feature = "std"), no_std)]
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, private_in_public, variant_size_differences)]
#![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
#![deny(unused_attributes, unused_imports, unused_mut)]
#![deny(renamed_and_removed_lints, stable_features, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use, const_err)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate bench_utils;

#[cfg(not(feature = "std"))]
#[macro_use]
extern crate alloc;

#[cfg(not(feature = "std"))]
use alloc::{string::String, vec::Vec};

#[cfg(feature = "std")]
use std::{string::String, vec::Vec};

use algebra::{bytes::ToBytes, PairingCurve, PairingEngine};
use r1cs_core::SynthesisError;

use algebra::io::{self, Read, Result as IoResult, Write};

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
#[derive(Clone)]
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

impl<E: PairingEngine> Proof<E> {
    /// Serialize the proof into bytes, for storage on disk or transmission
    /// over the network.
    pub fn write<W: Write>(&self, mut _writer: W) -> io::Result<()> {
        // TODO: implement serialization
        unimplemented!()
    }

    /// Deserialize the proof from bytes.
    pub fn read<R: Read>(mut _reader: R) -> io::Result<Self> {
        // TODO: implement serialization
        unimplemented!()
    }
}

/// A verification key in the Groth16 SNARK.
#[derive(Clone)]
pub struct VerifyingKey<E: PairingEngine> {
    pub alpha_g1:     E::G1Affine,
    pub beta_g2:      E::G2Affine,
    pub gamma_g2:     E::G2Affine,
    pub delta_g2:     E::G2Affine,
    pub gamma_abc_g1: Vec<E::G1Affine>,
}

impl<E: PairingEngine> ToBytes for VerifyingKey<E> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.alpha_g1.write(&mut writer)?;
        self.beta_g2.write(&mut writer)?;
        self.gamma_g2.write(&mut writer)?;
        self.delta_g2.write(&mut writer)?;
        for q in &self.gamma_abc_g1 {
            q.write(&mut writer)?;
        }
        Ok(())
    }
}

impl<E: PairingEngine> Default for VerifyingKey<E> {
    fn default() -> Self {
        Self {
            alpha_g1:     E::G1Affine::default(),
            beta_g2:      E::G2Affine::default(),
            gamma_g2:     E::G2Affine::default(),
            delta_g2:     E::G2Affine::default(),
            gamma_abc_g1: Vec::new(),
        }
    }
}

impl<E: PairingEngine> PartialEq for VerifyingKey<E> {
    fn eq(&self, other: &Self) -> bool {
        self.alpha_g1 == other.alpha_g1
            && self.beta_g2 == other.beta_g2
            && self.gamma_g2 == other.gamma_g2
            && self.delta_g2 == other.delta_g2
            && self.gamma_abc_g1 == other.gamma_abc_g1
    }
}

impl<E: PairingEngine> VerifyingKey<E> {
    /// Serialize the verification key into bytes, for storage on disk
    /// or transmission over the network.
    pub fn write<W: Write>(&self, mut _writer: W) -> io::Result<()> {
        // TODO: implement serialization
        unimplemented!()
    }

    /// Deserialize the verification key from bytes.
    pub fn read<R: Read>(mut _reader: R) -> io::Result<Self> {
        // TODO: implement serialization
        unimplemented!()
    }
}

/// Full public (prover and verifier) parameters for the Groth16 zkSNARK.
#[derive(Clone)]
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

impl<E: PairingEngine> Parameters<E> {
    /// Serialize the parameters to bytes.
    pub fn write<W: Write>(&self, mut _writer: W) -> io::Result<()> {
        // TODO: implement serialization
        unimplemented!()
    }

    /// Deserialize the public parameters from bytes.
    pub fn read<R: Read>(mut _reader: R, _checked: bool) -> io::Result<Self> {
        // TODO: implement serialization
        unimplemented!()
    }
}

/// Preprocessed verification key parameters that enable faster verification
/// at the expense of larger size in memory.
#[derive(Clone)]
pub struct PreparedVerifyingKey<E: PairingEngine> {
    pub vk:               VerifyingKey<E>,
    pub alpha_g1_beta_g2: E::Fqk,
    pub gamma_g2_neg_pc:  <E::G2Affine as PairingCurve>::Prepared,
    pub delta_g2_neg_pc:  <E::G2Affine as PairingCurve>::Prepared,
    pub gamma_abc_g1:     Vec<E::G1Affine>,
}

impl<E: PairingEngine> From<PreparedVerifyingKey<E>> for VerifyingKey<E> {
    fn from(other: PreparedVerifyingKey<E>) -> Self {
        other.vk
    }
}

impl<E: PairingEngine> From<VerifyingKey<E>> for PreparedVerifyingKey<E> {
    fn from(other: VerifyingKey<E>) -> Self {
        prepare_verifying_key(&other)
    }
}

impl<E: PairingEngine> Default for PreparedVerifyingKey<E> {
    fn default() -> Self {
        Self {
            vk:               VerifyingKey::default(),
            alpha_g1_beta_g2: E::Fqk::default(),
            gamma_g2_neg_pc:  <E::G2Affine as PairingCurve>::Prepared::default(),
            delta_g2_neg_pc:  <E::G2Affine as PairingCurve>::Prepared::default(),
            gamma_abc_g1:     Vec::new(),
        }
    }
}

impl<E: PairingEngine> ToBytes for PreparedVerifyingKey<E> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.vk.write(&mut writer)?;
        self.alpha_g1_beta_g2.write(&mut writer)?;
        self.gamma_g2_neg_pc.write(&mut writer)?;
        self.delta_g2_neg_pc.write(&mut writer)?;
        for q in &self.gamma_abc_g1 {
            q.write(&mut writer)?;
        }
        Ok(())
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
