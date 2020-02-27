mod error;
pub use error::*;

pub trait CanonicalSerialize {
    fn serialize(
        &self,
        extra_info: &[bool],
        output_buf: &mut [u8],
    ) -> Result<(), SerializationError>;
    fn buffer_size() -> usize;
}

pub trait CanonicalDeserialize: CanonicalSerialize {
    fn deserialize(bytes: &[u8], extra_info_buf: &mut [bool]) -> Result<Self, SerializationError>
    where
        Self: Sized;
}

pub fn buffer_bit_byte_size(modulus_bits: usize) -> (usize, usize) {
    let byte_size = (modulus_bits + 7) / 8;
    ((byte_size * 8), byte_size)
}

macro_rules! impl_prime_field_serializer {
    ($field: ident, $params: ident, $byte_size: expr) => {
        impl<P: $params> CanonicalSerialize for $field<P> {
            fn serialize(
                &self,
                extra_info: &[bool],
                output_buf: &mut [u8],
            ) -> Result<(), crate::serialize::SerializationError> {
                const BYTE_SIZE: usize = $byte_size;

                let mut bytes = [0u8; BYTE_SIZE];
                self.write(&mut bytes[..])?;

                let (output_bit_size, output_byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                if output_byte_size != output_buf.len() {
                    return Err(crate::serialize::SerializationError::BufferWrongSize);
                }
                let extra_info_len = extra_info.len();
                if extra_info_len > (output_bit_size - P::MODULUS_BITS as usize) {
                    return Err(crate::serialize::SerializationError::NotEnoughSpace);
                }

                for i in 0..extra_info_len {
                    if extra_info[i] {
                        bytes[output_byte_size - 1] |= 1 << (8 - extra_info_len + i);
                    }
                }

                output_buf.copy_from_slice(&bytes[..output_byte_size]);
                Ok(())
            }

            fn buffer_size() -> usize {
                let (_, byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                byte_size
            }
        }

        impl<P: $params> CanonicalDeserialize for $field<P> {
            fn deserialize(
                bytes: &[u8],
                extra_info_buf: &mut [bool],
            ) -> Result<Self, crate::serialize::SerializationError>
            where
                Self: Sized,
            {
                const BYTE_SIZE: usize = $byte_size;
                let (input_bit_size, input_byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                if input_byte_size != bytes.len() {
                    return Err(crate::serialize::SerializationError::BufferWrongSize);
                }
                let extra_info_len = extra_info_buf.len();
                if extra_info_len > (input_bit_size - P::MODULUS_BITS as usize) {
                    return Err(crate::serialize::SerializationError::NotEnoughSpace);
                }

                let mut masked_bytes = [0; BYTE_SIZE];
                masked_bytes[..input_byte_size].copy_from_slice(bytes);
                for i in 0..extra_info_len {
                    let bitmask = 1 << (8 - extra_info_len + i);
                    extra_info_buf[i] = bytes[input_byte_size - 1] & bitmask == bitmask;
                    masked_bytes[input_byte_size - 1] &= 0xFF - bitmask;
                }
                Ok(Self::read(&masked_bytes[..])?)
            }
        }
    };
}

macro_rules! impl_sw_curve_serializer {
    ($params: ident) => {
        impl<P: $params> CanonicalSerialize for GroupAffine<P> {
            fn serialize(
                &self,
                extra_info_buf: &[bool],
                output_buf: &mut [u8],
            ) -> Result<(), crate::serialize::SerializationError> {
                if extra_info_buf.len() != 0 {
                    return Err(crate::serialize::SerializationError::ExtraInfoWrongSize);
                }
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

            fn buffer_size() -> usize {
                <P::BaseField as CanonicalSerialize>::buffer_size()
            }
        }

        impl<P: $params> CanonicalDeserialize for GroupAffine<P> {
            fn deserialize(
                bytes: &[u8],
                extra_info_buf: &mut [bool],
            ) -> Result<Self, crate::serialize::SerializationError>
            where
                Self: Sized,
            {
                if extra_info_buf.len() != 0 {
                    return Err(crate::serialize::SerializationError::ExtraInfoWrongSize);
                }
                let mut extra_info = [false; 2];
                let x = P::BaseField::deserialize(bytes, &mut extra_info)?;
                if extra_info[1] {
                    return Ok(Self::zero());
                }
                let p = GroupAffine::<P>::get_point_from_x(x, extra_info[0])
                    .ok_or(crate::serialize::SerializationError::InvalidData)?;
                if !p.is_in_correct_subgroup_assuming_on_curve() {
                    return Err(crate::serialize::SerializationError::InvalidData);
                }
                Ok(p)
            }
        }
    };
}

macro_rules! impl_edwards_curve_serializer {
    ($params: ident) => {
        impl<P: $params> CanonicalSerialize for GroupAffine<P> {
            fn serialize(
                &self,
                extra_info_buf: &[bool],
                output_buf: &mut [u8],
            ) -> Result<(), crate::serialize::SerializationError> {
                if extra_info_buf.len() != 0 {
                    return Err(crate::serialize::SerializationError::ExtraInfoWrongSize);
                }
                if self.is_zero() {
                    P::BaseField::zero().serialize(&[false], output_buf)
                } else {
                    if self.y > -self.y {
                        self.x.serialize(&[true], output_buf)
                    } else {
                        self.x.serialize(&[false], output_buf)
                    }
                }
            }

            fn buffer_size() -> usize {
                <P::BaseField as CanonicalSerialize>::buffer_size()
            }
        }

        impl<P: $params> CanonicalDeserialize for GroupAffine<P> {
            fn deserialize(
                bytes: &[u8],
                extra_info_buf: &mut [bool],
            ) -> Result<Self, crate::serialize::SerializationError>
            where
                Self: Sized,
            {
                if extra_info_buf.len() != 0 {
                    return Err(crate::serialize::SerializationError::ExtraInfoWrongSize);
                }
                let mut extra_info = [false; 1];
                let x = P::BaseField::deserialize(bytes, &mut extra_info)?;
                if x == P::BaseField::zero() {
                    return Ok(Self::zero());
                }
                let p = GroupAffine::<P>::get_point_from_x(x, extra_info[0])
                    .ok_or(crate::serialize::SerializationError::InvalidData)?;
                if !p.is_in_correct_subgroup_assuming_on_curve() {
                    return Err(crate::serialize::SerializationError::InvalidData);
                }
                Ok(p)
            }
        }
    };
}

#[cfg(test)]
mod test {
    // use crate::fields::{Fp256, Fp256Parameters};
    // use crate::serialize::{CanonicalSerialize, SerializationError};
    // use crate::curves::bls12_377::{G1Affine, g1::Bls12_377G1Parameters};
}
