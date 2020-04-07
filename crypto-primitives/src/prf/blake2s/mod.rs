use alloc::vec::Vec;
use blake2::{Blake2s as B2s, VarBlake2s};
use digest::Digest;

use super::PRF;
use crate::CryptoError;

#[cfg(feature = "r1cs")]
pub mod constraints;

#[derive(Clone)]
pub struct Blake2s;

impl PRF for Blake2s {
    type Input = [u8; 32];
    type Output = [u8; 32];
    type Seed = [u8; 32];

    fn evaluate(seed: &Self::Seed, input: &Self::Input) -> Result<Self::Output, CryptoError> {
        let eval_time = start_timer!(|| "Blake2s::Eval");
        let mut h = B2s::new();
        h.input(seed.as_ref());
        h.input(input.as_ref());
        let mut result = [0u8; 32];
        result.copy_from_slice(&h.result());
        end_timer!(eval_time);
        Ok(result)
    }
}

#[derive(Clone)]
pub struct Blake2sWithParameterBlock {
    pub digest_length: u8,
    pub key_length: u8,
    pub fan_out: u8,
    pub depth: u8,
    pub leaf_length: u32,
    pub node_offset: u32,
    pub xof_digest_length: u16,
    pub node_depth: u8,
    pub inner_length: u8,
    pub salt: [u8; 8],
    pub personalization: [u8; 8],
}

impl Blake2sWithParameterBlock {
    pub fn parameters(&self) -> [u32; 8] {
        let mut parameters = [0; 8];
        parameters[0] = u32::from_le_bytes([
            self.digest_length,
            self.key_length,
            self.fan_out,
            self.depth,
        ]);
        parameters[1] = self.leaf_length;
        parameters[2] = self.node_offset;
        parameters[3] = u32::from_le_bytes([
            self.xof_digest_length as u8,
            (self.xof_digest_length >> 8) as u8,
            self.node_depth,
            self.inner_length,
        ]);

        let mut salt_bytes_1 = [0; 4];
        let mut salt_bytes_2 = [0; 4];
        let mut personalization_bytes_1 = [0; 4];
        let mut personalization_bytes_2 = [0; 4];
        for i in 0..4 {
            salt_bytes_1[i] = self.salt[i];
            salt_bytes_2[i] = self.salt[4 + i];
            personalization_bytes_1[i] = self.personalization[i];
            personalization_bytes_2[i] = self.personalization[4 + i];
        }
        parameters[4] = u32::from_le_bytes(salt_bytes_1);
        parameters[5] = u32::from_le_bytes(salt_bytes_2);
        parameters[6] = u32::from_le_bytes(personalization_bytes_1);
        parameters[7] = u32::from_le_bytes(personalization_bytes_2);

        parameters
    }

    pub fn evaluate(&self, input: &[u8]) -> Vec<u8> {
        use digest::*;
        let eval_time = start_timer!(|| "Blake2sWithParameterBlock::Eval");
        let mut h = VarBlake2s::with_parameter_block(&self.parameters());
        h.input(input.as_ref());
        end_timer!(eval_time);
        let mut buf = Vec::with_capacity(h.output_size());
        h.variable_result(|res| buf.extend_from_slice(res));
        buf
    }
}
