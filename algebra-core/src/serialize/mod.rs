mod error;
mod flags;
use crate::io::{Read, Write};
use core::cmp;
pub use error::*;
pub use flags::*;

#[cfg(not(feature = "std"))]
use alloc::vec::Vec;

#[cfg(feature = "std")]
use std::vec::Vec;

/// Serializer in little endian format allowing to encode flags.
pub trait CanonicalSerializeWithFlags: CanonicalSerialize {
    /// Serializes `self` and `flags` into `writer`.
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        writer: &mut W,
        flags: F,
    ) -> Result<(), SerializationError>;
}

/// Serializer in little endian format.
pub trait CanonicalSerialize {
    /// Serializes `self` into `writer`.
    fn serialize<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError>;
    fn serialized_size(&self) -> usize;

    /// Serializes `self` into `writer` without compression.
    fn serialize_uncompressed<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        self.serialize(writer)
    }
    fn uncompressed_size(&self) -> usize {
        self.serialized_size()
    }
}

/// Deserializer in little endian format allowing flags to be encoded.
pub trait CanonicalDeserializeWithFlags: Sized {
    /// Reads `Self` and `Flags` from `reader`.
    /// Returns empty flags by default.
    fn deserialize_with_flags<R: Read, F: Flags>(
        reader: &mut R,
    ) -> Result<(Self, F), SerializationError>;
}

/// Deserializer in little endian format.
pub trait CanonicalDeserialize: Sized {
    /// Reads `Self` from `reader`.
    fn deserialize<R: Read>(reader: &mut R) -> Result<Self, SerializationError>;

    /// Reads `Self` from `reader` without compression.
    fn deserialize_uncompressed<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        Self::deserialize(reader)
    }
}

impl CanonicalSerialize for u64 {
    fn serialize<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        Ok(writer.write_all(&self.to_le_bytes())?)
    }

    fn serialized_size(&self) -> usize {
        8
    }
}

impl CanonicalDeserialize for u64 {
    fn deserialize<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; 8];
        reader.read_exact(&mut bytes)?;
        Ok(u64::from_le_bytes(bytes))
    }
}

impl<T: CanonicalSerialize> CanonicalSerialize for Vec<T> {
    fn serialize<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize(writer)?;
        for item in self.iter() {
            item.serialize(writer)?;
        }
        Ok(())
    }

    fn serialized_size(&self) -> usize {
        8 + self
            .iter()
            .map(|item| item.serialized_size())
            .sum::<usize>()
    }

    fn serialize_uncompressed<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize(writer)?;
        for item in self.iter() {
            item.serialize_uncompressed(writer)?;
        }
        Ok(())
    }
    fn uncompressed_size(&self) -> usize {
        8 + self
            .iter()
            .map(|item| item.uncompressed_size())
            .sum::<usize>()
    }
}

impl<T: CanonicalDeserialize> CanonicalDeserialize for Vec<T> {
    fn deserialize<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        let len = u64::deserialize(reader)?;
        let mut values = vec![];
        for _ in 0..len {
            values.push(T::deserialize(reader)?);
        }
        Ok(values)
    }

    fn deserialize_uncompressed<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        let len = u64::deserialize(reader)?;
        let mut values = vec![];
        for _ in 0..len {
            values.push(T::deserialize_uncompressed(reader)?);
        }
        Ok(values)
    }
}

pub fn buffer_bit_byte_size(modulus_bits: usize) -> (usize, usize) {
    let byte_size = (modulus_bits + 7) / 8;
    ((byte_size * 8), byte_size)
}

pub(crate) fn serialize_num_limbs<W: Write, F: Flags>(
    writer: &mut W,
    limbs: &[u64],
    mut num_bytes: usize,
    flags: F,
) -> Result<(), SerializationError> {
    // Calculate number of limbs to encode to write `num_bytes`.
    let num_limbs = (num_bytes + 7) / 8;
    for (i, limb) in limbs.iter().take(num_limbs).enumerate() {
        // Convert each limb to little-endian.
        let mut bytes: [u8; 8] = limb.to_le_bytes();
        // Calculate the number of bytes to write.
        let bytes_to_write = cmp::min(8, num_bytes);

        // Encode flags into last limb, last byte.
        if i == num_limbs - 1 {
            bytes[bytes_to_write - 1] |= flags.u8_bitmask();
        }

        // Write bytes.
        writer.write_all(&bytes[..bytes_to_write])?;
        num_bytes -= bytes_to_write;
    }
    Ok(())
}

pub(crate) fn deserialize_num_limbs<R: Read, F: Flags>(
    reader: &mut R,
    limbs: &mut [u64],
    mut num_bytes: usize,
    remove_flags: bool,
) -> Result<F, SerializationError> {
    // Calculate number of limbs to encode to write `num_bytes`.
    let num_limbs = (num_bytes + 7) / 8;
    let mut flags = F::default();
    for (i, limb) in limbs.iter_mut().take(num_limbs).enumerate() {
        // Calculate the number of bytes to read.
        let bytes_to_read = cmp::min(8, num_bytes);

        // Read exactly bytes_to_read bytes.
        let mut bytes = [0u8; 8];
        reader.read_exact(&mut bytes[..bytes_to_read])?;

        // Read/Remove flags from last limb, last byte.
        if i == num_limbs - 1 {
            if remove_flags {
                flags = F::from_u8_remove_flags(&mut bytes[bytes_to_read - 1]);
            } else {
                flags = F::from_u8(bytes[bytes_to_read - 1]);
            }
        }

        // Convert from little-endian.
        *limb = u64::from_le_bytes(bytes);
        num_bytes -= bytes_to_read;
    }
    Ok(flags)
}

macro_rules! impl_prime_field_serializer {
    ($field: ident, $params: ident, $byte_size: expr) => {
        impl<P: $params> CanonicalSerializeWithFlags for $field<P> {
            #[allow(unused_qualifications)]
            fn serialize_with_flags<W: crate::io::Write, F: crate::serialize::Flags>(
                &self,
                writer: &mut W,
                flags: F,
            ) -> Result<(), crate::serialize::SerializationError> {
                let (output_bit_size, output_byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                if F::len() > (output_bit_size - P::MODULUS_BITS as usize) {
                    return Err(crate::serialize::SerializationError::NotEnoughSpace);
                }

                let repr = self.into_repr();
                crate::serialize::serialize_num_limbs(
                    writer,
                    repr.as_ref(),
                    output_byte_size,
                    flags,
                )
            }
        }

        impl<P: $params> CanonicalSerialize for $field<P> {
            #[allow(unused_qualifications)]
            fn serialize<W: crate::io::Write>(
                &self,
                writer: &mut W,
            ) -> Result<(), crate::serialize::SerializationError> {
                self.serialize_with_flags(writer, crate::serialize::EmptyFlags)
            }

            fn serialized_size(&self) -> usize {
                let (_, output_byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                output_byte_size
            }
        }

        impl<P: $params> CanonicalDeserializeWithFlags for $field<P> {
            #[allow(unused_qualifications)]
            fn deserialize_with_flags<R: crate::io::Read, F: crate::serialize::Flags>(
                reader: &mut R,
            ) -> Result<(Self, F), crate::serialize::SerializationError> {
                let (output_bit_size, output_byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                if F::len() > (output_bit_size - P::MODULUS_BITS as usize) {
                    return Err(crate::serialize::SerializationError::NotEnoughSpace);
                }

                let mut value = P::BigInt::default();
                let flags = crate::serialize::deserialize_num_limbs(
                    reader,
                    value.as_mut(),
                    output_byte_size,
                    true,
                )?;
                Ok((Self::from_repr(value), flags))
            }
        }

        impl<P: $params> CanonicalDeserialize for $field<P> {
            #[allow(unused_qualifications)]
            fn deserialize<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let (_, output_byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                let mut value = P::BigInt::default();
                crate::serialize::deserialize_num_limbs::<_, crate::serialize::EmptyFlags>(
                    reader,
                    value.as_mut(),
                    output_byte_size,
                    false,
                )?;
                Ok(Self::from_repr(value))
            }
        }
    };
}

macro_rules! impl_sw_curve_serializer {
    ($params: ident) => {
        impl<P: $params> CanonicalSerialize for GroupAffine<P> {
            #[allow(unused_qualifications)]
            fn serialize<W: crate::io::Write>(
                &self,
                writer: &mut W,
            ) -> Result<(), crate::serialize::SerializationError> {
                if self.is_zero() {
                    let flags = crate::serialize::SWFlags::infinity();
                    // Serialize 0.
                    P::BaseField::zero().serialize_with_flags(writer, flags)
                } else {
                    let flags = crate::serialize::SWFlags::from_y_sign(self.y > -self.y);
                    self.x.serialize_with_flags(writer, flags)
                }
            }

            fn serialized_size(&self) -> usize {
                self.x.serialized_size()
            }

            #[allow(unused_qualifications)]
            fn serialize_uncompressed<W: crate::io::Write>(
                &self,
                writer: &mut W,
            ) -> Result<(), crate::serialize::SerializationError> {
                let flags = if self.is_zero() {
                    crate::serialize::SWFlags::infinity()
                } else {
                    crate::serialize::SWFlags::default()
                };
                self.x.serialize(writer)?;
                self.y.serialize_with_flags(writer, flags)?;
                Ok(())
            }

            fn uncompressed_size(&self) -> usize {
                self.x.serialized_size() + self.y.serialized_size()
            }
        }

        impl<P: $params> CanonicalDeserialize for GroupAffine<P> {
            #[allow(unused_qualifications)]
            fn deserialize<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let (x, flags): (P::BaseField, crate::serialize::SWFlags) =
                    CanonicalDeserializeWithFlags::deserialize_with_flags(reader)?;
                if flags.is_infinity() {
                    Ok(Self::zero())
                } else {
                    let p = GroupAffine::<P>::get_point_from_x(x, flags.is_positive().unwrap())
                        .ok_or(crate::serialize::SerializationError::InvalidData)?;
                    if !p.is_in_correct_subgroup_assuming_on_curve() {
                        return Err(crate::serialize::SerializationError::InvalidData);
                    }
                    Ok(p)
                }
            }

            #[allow(unused_qualifications)]
            fn deserialize_uncompressed<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let x: P::BaseField = CanonicalDeserialize::deserialize(reader)?;
                let (y, flags): (P::BaseField, crate::serialize::SWFlags) =
                    CanonicalDeserializeWithFlags::deserialize_with_flags(reader)?;

                let p = GroupAffine::<P>::new(x, y, flags.is_infinity());
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
            #[allow(unused_qualifications)]
            fn serialize<W: crate::io::Write>(
                &self,
                writer: &mut W,
            ) -> Result<(), crate::serialize::SerializationError> {
                if self.is_zero() {
                    let flags = crate::serialize::EdwardsFlags::default();
                    // Serialize 0.
                    P::BaseField::zero().serialize_with_flags(writer, flags)
                } else {
                    let flags = crate::serialize::EdwardsFlags::from_y_sign(self.y > -self.y);
                    self.x.serialize_with_flags(writer, flags)
                }
            }

            fn serialized_size(&self) -> usize {
                CanonicalSerialize::serialized_size(&self.x)
            }

            #[allow(unused_qualifications)]
            fn serialize_uncompressed<W: crate::io::Write>(
                &self,
                writer: &mut W,
            ) -> Result<(), crate::serialize::SerializationError> {
                self.x.serialize_uncompressed(writer)?;
                self.y.serialize_uncompressed(writer)?;
                Ok(())
            }

            fn uncompressed_size(&self) -> usize {
                self.x.uncompressed_size() + self.y.uncompressed_size()
            }
        }

        impl<P: $params> CanonicalDeserialize for GroupAffine<P> {
            #[allow(unused_qualifications)]
            fn deserialize<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let (x, flags): (P::BaseField, crate::serialize::EdwardsFlags) =
                    CanonicalDeserializeWithFlags::deserialize_with_flags(reader)?;
                if x == P::BaseField::zero() {
                    Ok(Self::zero())
                } else {
                    let p = GroupAffine::<P>::get_point_from_x(x, flags.is_positive())
                        .ok_or(crate::serialize::SerializationError::InvalidData)?;
                    if !p.is_in_correct_subgroup_assuming_on_curve() {
                        return Err(crate::serialize::SerializationError::InvalidData);
                    }
                    Ok(p)
                }
            }

            #[allow(unused_qualifications)]
            fn deserialize_uncompressed<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let x: P::BaseField = CanonicalDeserialize::deserialize(reader)?;
                let y: P::BaseField = CanonicalDeserialize::deserialize(reader)?;

                let p = GroupAffine::<P>::new(x, y);
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
    use crate::{io::Cursor, CanonicalDeserialize, CanonicalSerialize};

    #[test]
    fn test_primitives() {
        let a = 192830918u64;
        let mut serialized = vec![0u8; a.serialized_size()];
        let mut cursor = Cursor::new(&mut serialized[..]);
        a.serialize(&mut cursor).unwrap();

        let mut cursor = Cursor::new(&serialized[..]);
        let b = u64::deserialize(&mut cursor).unwrap();
        assert_eq!(a, b);
    }
}
