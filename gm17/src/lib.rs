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
use alloc::{string::String, vec::Vec};

#[cfg(feature = "std")]
use std::{string::String, vec::Vec};

use algebra_core::{
    bytes::{FromBytes, ToBytes},
    io::{self, Result as IoResult},
    serialize::*,
    PairingEngine,
};
use r1cs_core::SynthesisError;

// use byteorder::{BigEndian, ReadBytesExt};

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
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
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

// impl<E: PairingEngine> Proof<E> {
//     /// Serialize the proof into bytes, for storage on disk or transmission
//     /// over the network.
//     pub fn write<W: Write>(&self, mut _writer: W) -> io::Result<()> {
//         // TODO: implement serialization
//         unimplemented!()
//     }
//
//     /// Deserialize the proof from bytes.
//     pub fn read<R: Read>(mut _reader: R) -> io::Result<Self> {
//         // TODO: implement serialization
//         unimplemented!()
//     }
// }

/// A verification key in the GM17 SNARK.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, Debug)]
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

        write_slice_with_len(&mut writer, &self.query)?;
        // for q in &self.query {
        //     q.write(&mut writer)?;
        // }
        Ok(())
    }
}

impl<E: PairingEngine> FromBytes for VerifyingKey<E> {
    fn read<R: Read>(reader: R) -> IoResult<Self> {
        let mut reader = reader;
        // use num_traits::Zero;

        let result_transform_g1 = |result: io::Result<E::G1Affine>| -> io::Result<E::G1Affine> {
            result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            // .and_then(|e| {
            //     if e.is_zero() {
            //         Err(io::Error::new(
            //             io::ErrorKind::InvalidData,
            //             "point at infinity",
            //         ))
            //     } else {
            //         Ok(e)
            //     }
            // })
        };
        let read_g1 = |reader: &mut R| -> io::Result<E::G1Affine> {
            result_transform_g1(<E::G1Affine>::read(reader))
        };

        let result_transform_g2 = |result: io::Result<E::G2Affine>| -> io::Result<E::G2Affine> {
            result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            // .and_then(|e| {
            //     if e.is_zero() {
            //         Err(io::Error::new(
            //             io::ErrorKind::InvalidData,
            //             "point at infinity",
            //         ))
            //     } else {
            //         Ok(e)
            //     }
            // })
        };
        let read_g2 = |reader: &mut R| -> io::Result<E::G2Affine> {
            result_transform_g2(<E::G2Affine>::read(reader))
        };

        let h_g2 = read_g2(&mut reader)?;
        let g_alpha_g1 = read_g1(&mut reader)?;
        let h_beta_g2 = read_g2(&mut reader)?;
        let g_gamma_g1 = read_g1(&mut reader)?;
        let h_gamma_g2 = read_g2(&mut reader)?;
        let query = read_vec_with_len(&mut reader, result_transform_g1)?;

        // {
        //     let len = reader.read_u32::<BigEndian>()? as usize;
        //     for _ in 0..len {
        //         query.push(read_g1(&mut reader)?);
        //     }
        // }

        Ok(VerifyingKey {
            h_g2,
            g_alpha_g1,
            h_beta_g2,
            g_gamma_g1,
            h_gamma_g2,
            query,
        })
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

impl<E: PairingEngine> PartialEq for VerifyingKey<E> {
    fn eq(&self, other: &Self) -> bool {
        self.h_g2 == other.h_g2
            && self.g_alpha_g1 == other.g_alpha_g1
            && self.h_beta_g2 == other.h_beta_g2
            && self.g_gamma_g1 == other.g_gamma_g1
            && self.h_gamma_g2 == other.h_gamma_g2
            && self.query == other.query
    }
}

// impl<E: PairingEngine> VerifyingKey<E> {
//     /// Serialize the verification key into bytes, for storage on disk
//     /// or transmission over the network.
//     pub fn write<W: Write>(&self, mut _writer: W) -> io::Result<()> {
//         // TODO: implement serialization
//         unimplemented!()
//     }
//
//     /// Deserialize the verification key from bytes.
//     pub fn read<R: Read>(mut _reader: R) -> io::Result<Self> {
//         // TODO: implement serialization
//         unimplemented!()
//     }
// }

/// Full public (prover and verifier) parameters for the GM17 zkSNARK.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize, Debug)]
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

fn write_slice_with_len<W: Write, T: ToBytes>(mut writer: W, data: &[T]) -> IoResult<()> {
    (data.len() as u32).write(&mut writer)?;
    for d in data {
        d.write(&mut writer)?;
    }
    Ok(())
}
fn read_vec_with_len<R: Read, T: FromBytes, F: FnMut(IoResult<T>) -> IoResult<T>>(
    mut reader: R,
    mut result_transform: F,
) -> IoResult<Vec<T>> {
    let len = u32::read(&mut reader)? as usize;
    (0..len)
        .map(|_| result_transform(T::read(&mut reader)))
        .collect()
}

impl<E: PairingEngine> ToBytes for Parameters<E> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        self.vk.write(&mut writer)?;
        write_slice_with_len(&mut writer, &self.a_query)?;
        write_slice_with_len(&mut writer, &self.b_query)?;
        write_slice_with_len(&mut writer, &self.c_query_1)?;
        write_slice_with_len(&mut writer, &self.c_query_2)?;
        // for a_q in &self.a_query {
        //     a_q.write(&mut writer)?;
        // }
        // for b_q in &self.b_query {
        //     b_q.write(&mut writer)?;
        // }
        // for c_q1 in &self.c_query_1 {
        //     c_q1.write(&mut writer)?;
        // }
        // for c_q2 in &self.c_query_2 {
        //     c_q2.write(&mut writer)?;
        // }
        self.g_gamma_z.write(&mut writer)?;
        self.h_gamma_z.write(&mut writer)?;
        self.g_ab_gamma_z.write(&mut writer)?;
        self.g_gamma2_z2.write(&mut writer)?;
        write_slice_with_len(&mut writer, &self.g_gamma2_z_t)?;
        // for g in &self.g_gamma2_z_t {
        //     g.write(&mut writer)?;
        // }
        Ok(())
    }
}

impl<E: PairingEngine> FromBytes for Parameters<E> {
    fn read<R: Read>(reader: R) -> IoResult<Self> {
        let mut reader = reader;
        // use num_traits::Zero;

        let result_transform_g1 = |result: io::Result<E::G1Affine>| -> io::Result<E::G1Affine> {
            result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            // .and_then(|e| {
            //     if e.is_zero() {
            //         Err(io::Error::new(
            //             io::ErrorKind::InvalidData,
            //             "point at infinity",
            //         ))
            //     } else {
            //         Ok(e)
            //     }
            // })
        };
        let read_g1 = |reader: &mut R| -> io::Result<E::G1Affine> {
            result_transform_g1(<E::G1Affine>::read(reader))
        };

        let result_transform_g2 = |result: io::Result<E::G2Affine>| -> io::Result<E::G2Affine> {
            result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            // .and_then(|e| {
            //     if e.is_zero() {
            //         Err(io::Error::new(
            //             io::ErrorKind::InvalidData,
            //             "point at infinity",
            //         ))
            //     } else {
            //         Ok(e)
            //     }
            // })
        };
        let read_g2 = |reader: &mut R| -> io::Result<E::G2Affine> {
            result_transform_g2(<E::G2Affine>::read(reader))
        };

        let vk = VerifyingKey::<E>::read(&mut reader)?;

        let a_query = read_vec_with_len(&mut reader, result_transform_g1)?;
        let b_query = read_vec_with_len(&mut reader, result_transform_g2)?;
        let c_query_1 = read_vec_with_len(&mut reader, result_transform_g1)?;
        let c_query_2 = read_vec_with_len(&mut reader, result_transform_g1)?;

        // {
        //     let len = reader.read_u32::<BigEndian>()? as usize;
        //     for _ in 0..len {
        //         a_query.push(read_g1(&mut reader)?);
        //     }
        // }
        //
        // {
        //     let len = reader.read_u32::<BigEndian>()? as usize;
        //     for _ in 0..len {
        //         b_query.push(read_g2(&mut reader)?);
        //     }
        // }
        //
        // {
        //     let len = reader.read_u32::<BigEndian>()? as usize;
        //     for _ in 0..len {
        //         c_query_1.push(read_g1(&mut reader)?);
        //     }
        // }
        //
        // {
        //     let len = reader.read_u32::<BigEndian>()? as usize;
        //     for _ in 0..len {
        //         c_query_2.push(read_g1(&mut reader)?);
        //     }
        // }

        let g_gamma_z = read_g1(&mut reader)?;
        let h_gamma_z = read_g2(&mut reader)?;
        let g_ab_gamma_z = read_g1(&mut reader)?;
        let g_gamma2_z2 = read_g1(&mut reader)?;

        let g_gamma2_z_t = read_vec_with_len(&mut reader, result_transform_g1)?;

        // {
        //     let len = reader.read_u32::<BigEndian>()? as usize;
        //     for _ in 0..len {
        //         g_gamma2_z_t.push(read_g1(&mut reader)?);
        //     }
        // }

        Ok(Parameters {
            vk,
            a_query,
            b_query,
            c_query_1,
            c_query_2,
            g_gamma_z,
            h_gamma_z,
            g_ab_gamma_z,
            g_gamma2_z2,
            g_gamma2_z_t,
        })
    }
}

impl<E: PairingEngine> PartialEq for Parameters<E> {
    fn eq(&self, other: &Self) -> bool {
        self.vk == other.vk
            && self.a_query == other.a_query
            && self.b_query == other.b_query
            && self.c_query_1 == other.c_query_1
            && self.c_query_2 == other.c_query_2
            && self.g_gamma_z == other.g_gamma_z
            && self.h_gamma_z == other.h_gamma_z
            && self.g_ab_gamma_z == other.g_ab_gamma_z
            && self.g_gamma2_z2 == other.g_gamma2_z2
            && self.g_gamma2_z_t == other.g_gamma2_z_t
    }
}

// impl<E: PairingEngine> Parameters<E> {
//     /// Serialize the parameters to bytes.
//     pub fn write<W: Write>(&self, mut _writer: W) -> io::Result<()> {
//         // TODO: implement serialization
//         unimplemented!()
//     }
//
//     /// Deserialize the public parameters from bytes.
//     pub fn read<R: Read>(mut _reader: R, _checked: bool) -> io::Result<Self>
// {         // TODO: implement serialization
//         unimplemented!()
//     }
// }

/// Preprocessed verification key parameters that enable faster verification
/// at the expense of larger size in memory.
#[derive(Clone)]
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
