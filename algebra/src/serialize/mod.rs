mod error;
pub use error::*;

pub trait CanonicalSerialize {
    fn serialize(&self, extra_info: &[bool], output_buf: &mut [u8]) -> Result<(), SerializationError>;
}

pub trait CanonicalDeserialize: CanonicalSerialize {
    fn deserialize(bytes: &[u8], extra_info_buf: &mut [bool]) -> Result<Self, SerializationError>
    where Self: Sized;
}

pub fn buffer_bit_byte_size(modulus_bits: usize) -> (usize, usize) {
    let byte_size = (modulus_bits + 63)/64;
    ((byte_size*64), byte_size*8)
}

macro_rules! impl_prime_field_serializer {
    ($field: ident, $params: ident, $byte_size: expr) => {
        impl<P: $params> CanonicalSerialize for $field<P> {
            fn serialize(&self, extra_info: &[bool], output_buf: &mut [u8]) -> Result<(), SerializationError> {
                const BYTE_SIZE: usize = $byte_size;
                const BIT_SIZE: usize = 8*BYTE_SIZE;

                let mut bytes = [0u8; BYTE_SIZE];
                self.write(&mut bytes[..])?;

                if BYTE_SIZE != output_buf.len() {
                    return Err(SerializationError::BufferWrongSize);
                }
                let extra_info_len = extra_info.len();
                if  extra_info_len > (BIT_SIZE - P::MODULUS_BITS as usize) {
                    return Err(SerializationError::NotEnoughSpace);
                }

                for i in 0..extra_info_len {
                    if extra_info[i] {
                        bytes[BYTE_SIZE - 1] |= 1 << (8 - extra_info_len + i);
                    }
                }
                output_buf.copy_from_slice(&bytes[..]);
                Ok(())
            }
        }

        impl<P: $params> CanonicalDeserialize for $field<P> {
            fn deserialize(bytes: &[u8], extra_info_buf: &mut [bool]) -> Result<Self, SerializationError>
                where Self: Sized {
                const BYTE_SIZE: usize = $byte_size;
                const BIT_SIZE: usize = 8*BYTE_SIZE;
                if BYTE_SIZE != bytes.len() {
                    return Err(SerializationError::BufferWrongSize);
                }
                let extra_info_len = extra_info_buf.len();
                if  extra_info_len > (BIT_SIZE - P::MODULUS_BITS as usize) {
                    return Err(SerializationError::NotEnoughSpace);
                }

                let mut masked_bytes = [0; BYTE_SIZE];
                masked_bytes[..BYTE_SIZE].copy_from_slice(bytes);
                for i in 0..extra_info_len {
                    let bitmask = 1 << (8 - extra_info_len + i);
                    extra_info_buf[i] = bytes[BYTE_SIZE - 1] & bitmask == bitmask;
                    masked_bytes[BYTE_SIZE - 1] &= 0xFF - bitmask;
                }
                Ok(Self::read(&masked_bytes[..])?)
            }
        }
    }
}

macro_rules! impl_curve_serializer {
    ($params: ident) => {
        impl<P: $params> CanonicalSerialize for GroupAffine<P> {
            fn serialize(&self, _: &[bool], output_buf: &mut [u8]) -> Result<(), SerializationError> {
                if self.is_zero() {
                    P::BaseField::zero().serialize(&[false, true], output_buf)
                } else {
                    if self.y > -self.y {
                        self.x.serialize(&[true, false], output_buf)
                    } else {
                        self.x.serialize(&[false, false], output_buf)
                    }
                }
            }
        }

        impl<P: $params> CanonicalDeserialize for GroupAffine<P> {
            fn deserialize(bytes: &[u8], _: &mut [bool]) -> Result<Self, SerializationError>
                where Self: Sized {
                let mut extra_info_buf = [false; 2];
                let x = P::BaseField::deserialize(bytes, &mut extra_info_buf)?;
                if extra_info_buf[1] {
                    return Ok(Self::zero())
                }
                GroupAffine::<P>::get_point_from_x(x, extra_info_buf[0])
                    .ok_or(SerializationError::InvalidData)
            }
        }
    }
}

#[cfg(test)]
mod test {
    //use crate::fields::{Fp256, Fp256Parameters};
    //use crate::serialize::{CanonicalSerialize, SerializationError};
    //use crate::curves::bls12_377::{G1Affine, g1::Bls12_377G1Parameters};
}
