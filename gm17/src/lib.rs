//! An implementation of the [`Groth-Maller`] simulation extractable zkSNARK.
//!
//! [`Groth-Maller`]: https://eprint.iacr.org/2017/540
#![cfg_attr(not(feature = "std"), no_std)]
#![deny(unused_import_braces, unused_qualifications, trivial_casts)]
#![deny(trivial_numeric_casts, private_in_public, variant_size_differences)]
#![deny(stable_features, unreachable_pub, non_shorthand_field_patterns)]
#![deny(unused_attributes, unused_imports, unused_mut)]
#![deny(renamed_and_removed_lints, unused_allocation)]
#![deny(unused_comparisons, bare_trait_objects, unused_must_use, const_err)]
#![forbid(unsafe_code)]

#[macro_use]
extern crate bench_utils;

#[cfg(not(feature = "std"))]
#[macro_use]
extern crate alloc;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

#[cfg(feature = "std")]
use std::vec::Vec;

use algebra_core::{
    bytes::ToBytes,
    io::{self, Result as IoResult},
    serialize::*,
    PairingEngine,
};
use r1cs_core::SynthesisError;

/// Reduce an R1CS instance to a *Square Arithmetic Program* instance.
pub mod r1cs_to_sap;

/// Generate public parameters for the GM17 zkSNARK construction.
pub mod generator;

/// Create proofs for the GM17 zkSNARK construction.
pub mod prover;

/// Verify proofs for the GM17 zkSNARK construction.
pub mod verifier;

#[cfg(test)]
mod test;

pub use self::{generator::*, prover::*, verifier::*};

/// A proof in the GM17 SNARK.
#[derive(PartialEq, Eq, Clone, CanonicalSerialize, CanonicalDeserialize)]
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

impl<E: PairingEngine> Default for Proof<E> {
    fn default() -> Self {
        Self {
            a: E::G1Affine::default(),
            b: E::G2Affine::default(),
            c: E::G1Affine::default(),
        }
    }
}

/// A verification key in the GM17 SNARK.
#[derive(Eq, PartialEq, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct VerifyingKey<E: PairingEngine> {
    pub h_g2: E::G2Affine,
    pub g_alpha_g1: E::G1Affine,
    pub h_beta_g2: E::G2Affine,
    pub g_gamma_g1: E::G1Affine,
    pub h_gamma_g2: E::G2Affine,
    pub query: Vec<E::G1Affine>,
}

impl<E: PairingEngine> ToBytes for VerifyingKey<E> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.h_g2.write(&mut writer)?;
        self.g_alpha_g1.write(&mut writer)?;
        self.h_beta_g2.write(&mut writer)?;
        self.g_gamma_g1.write(&mut writer)?;
        self.h_gamma_g2.write(&mut writer)?;
        for q in &self.query {
            q.write(&mut writer)?;
        }
        Ok(())
    }
}

impl<E: PairingEngine> Default for VerifyingKey<E> {
    fn default() -> Self {
        Self {
            h_g2: E::G2Affine::default(),
            g_alpha_g1: E::G1Affine::default(),
            h_beta_g2: E::G2Affine::default(),
            g_gamma_g1: E::G1Affine::default(),
            h_gamma_g2: E::G2Affine::default(),
            query: Vec::new(),
        }
    }
}

/// Full public (prover and verifier) parameters for the GM17 zkSNARK.
#[derive(PartialEq, Eq, Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Parameters<E: PairingEngine> {
    pub vk: VerifyingKey<E>,
    pub a_query: Vec<E::G1Affine>,
    pub b_query: Vec<E::G2Affine>,
    pub c_query_1: Vec<E::G1Affine>,
    pub c_query_2: Vec<E::G1Affine>,
    pub g_gamma_z: E::G1Affine,
    pub h_gamma_z: E::G2Affine,
    pub g_ab_gamma_z: E::G1Affine,
    pub g_gamma2_z2: E::G1Affine,
    pub g_gamma2_z_t: Vec<E::G1Affine>,
}

/// Preprocessed verification key parameters that enable faster verification
/// at the expense of larger size in memory.
#[derive(PartialEq, Eq, Clone)]
pub struct PreparedVerifyingKey<E: PairingEngine> {
    pub vk: VerifyingKey<E>,
    pub g_alpha: E::G1Affine,
    pub h_beta: E::G2Affine,
    pub g_alpha_h_beta_ml: E::Fqk,
    pub g_gamma_pc: E::G1Prepared,
    pub h_gamma_pc: E::G2Prepared,
    pub h_pc: E::G2Prepared,
    pub query: Vec<E::G1Affine>,
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
            vk: VerifyingKey::default(),
            g_alpha: E::G1Affine::default(),
            h_beta: E::G2Affine::default(),
            g_alpha_h_beta_ml: E::Fqk::default(),
            g_gamma_pc: E::G1Prepared::default(),
            h_gamma_pc: E::G2Prepared::default(),
            h_pc: E::G2Prepared::default(),
            query: Vec::new(),
        }
    }
}

impl<E: PairingEngine> ToBytes for PreparedVerifyingKey<E> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.vk.write(&mut writer)?;
        self.g_alpha.write(&mut writer)?;
        self.h_beta.write(&mut writer)?;
        self.g_alpha_h_beta_ml.write(&mut writer)?;
        self.g_gamma_pc.write(&mut writer)?;
        self.h_gamma_pc.write(&mut writer)?;
        self.h_pc.write(&mut writer)?;
        for q in &self.query {
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

    pub fn get_b_query(
        &self,
        num_inputs: usize,
    ) -> Result<(&[E::G2Affine], &[E::G2Affine]), SynthesisError> {
        Ok((&self.b_query[1..num_inputs], &self.b_query[num_inputs..]))
    }

    pub fn get_c_query_1(
        &self,
        num_inputs: usize,
    ) -> Result<(&[E::G1Affine], &[E::G1Affine]), SynthesisError> {
        Ok((
            &self.c_query_1[0..num_inputs],
            &self.c_query_1[num_inputs..],
        ))
    }

    pub fn get_c_query_2(
        &self,
        num_inputs: usize,
    ) -> Result<(&[E::G1Affine], &[E::G1Affine]), SynthesisError> {
        Ok((
            &self.c_query_2[1..num_inputs],
            &self.c_query_2[num_inputs..],
        ))
    }

    pub fn get_g_gamma_z(&self) -> Result<E::G1Affine, SynthesisError> {
        Ok(self.g_gamma_z)
    }

    pub fn get_h_gamma_z(&self) -> Result<E::G2Affine, SynthesisError> {
        Ok(self.h_gamma_z)
    }

    pub fn get_g_ab_gamma_z(&self) -> Result<E::G1Affine, SynthesisError> {
        Ok(self.g_ab_gamma_z)
    }

    pub fn get_g_gamma2_z2(&self) -> Result<E::G1Affine, SynthesisError> {
        Ok(self.g_gamma2_z2)
    }

    pub fn get_g_gamma2_z_t(
        &self,
        num_inputs: usize,
    ) -> Result<(&[E::G1Affine], &[E::G1Affine]), SynthesisError> {
        Ok((
            &self.g_gamma2_z_t[0..num_inputs],
            &self.g_gamma2_z_t[num_inputs..],
        ))
    }

    pub fn get_a_query_full(&self) -> Result<&[E::G1Affine], SynthesisError> {
        Ok(&self.a_query)
    }

    pub fn get_b_query_full(&self) -> Result<&[E::G2Affine], SynthesisError> {
        Ok(&self.b_query)
    }

    pub fn get_c_query_1_full(&self) -> Result<&[E::G1Affine], SynthesisError> {
        Ok(&self.c_query_1)
    }

    pub fn get_c_query_2_full(&self) -> Result<&[E::G1Affine], SynthesisError> {
        Ok(&self.c_query_2)
    }

    pub fn get_g_gamma2_z_t_full(&self) -> Result<&[E::G1Affine], SynthesisError> {
        Ok(&self.g_gamma2_z_t)
    }
}
