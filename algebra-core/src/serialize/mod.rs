mod error;
mod flags;
use crate::io::{Read, Write};
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
    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        self.serialize(writer)
    }
    #[inline]
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
    #[inline]
    fn deserialize_uncompressed<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        Self::deserialize(reader)
    }
}

impl CanonicalSerialize for u64 {
    #[inline]
    fn serialize<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        Ok(writer.write_all(&self.to_le_bytes())?)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        8
    }
}

impl CanonicalDeserialize for u64 {
    #[inline]
    fn deserialize<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; 8];
        reader.read_exact(&mut bytes)?;
        Ok(u64::from_le_bytes(bytes))
    }
}

impl<T: CanonicalSerialize> CanonicalSerialize for Vec<T> {
    #[inline]
    fn serialize<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize(writer)?;
        for item in self.iter() {
            item.serialize(writer)?;
        }
        Ok(())
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        8 + self
            .iter()
            .map(|item| item.serialized_size())
            .sum::<usize>()
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize(writer)?;
        for item in self.iter() {
            item.serialize_uncompressed(writer)?;
        }
        Ok(())
    }

    #[inline]
    fn uncompressed_size(&self) -> usize {
        8 + self
            .iter()
            .map(|item| item.uncompressed_size())
            .sum::<usize>()
    }
}

impl<T: CanonicalDeserialize> CanonicalDeserialize for Vec<T> {
    #[inline]
    fn deserialize<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        let len = u64::deserialize(reader)?;
        let mut values = vec![];
        for _ in 0..len {
            values.push(T::deserialize(reader)?);
        }
        Ok(values)
    }

    #[inline]
    fn deserialize_uncompressed<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        let len = u64::deserialize(reader)?;
        let mut values = vec![];
        for _ in 0..len {
            values.push(T::deserialize_uncompressed(reader)?);
        }
        Ok(values)
    }
}

#[inline]
pub fn buffer_bit_byte_size(modulus_bits: usize) -> (usize, usize) {
    let byte_size = (modulus_bits + 7) / 8;
    ((byte_size * 8), byte_size)
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
                const BYTE_SIZE: usize = $byte_size;

                let (output_bit_size, output_byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                if F::len() > (output_bit_size - P::MODULUS_BITS as usize) {
                    return Err(crate::serialize::SerializationError::NotEnoughSpace);
                }

                let mut bytes = [0u8; BYTE_SIZE];
                self.write(&mut bytes[..])?;

                bytes[output_byte_size - 1] |= flags.u8_bitmask();

                writer.write_all(&bytes[..output_byte_size])?;
                Ok(())
            }
        }

        impl<P: $params> CanonicalSerialize for $field<P> {
            #[allow(unused_qualifications)]
            #[inline]
            fn serialize<W: crate::io::Write>(
                &self,
                writer: &mut W,
            ) -> Result<(), crate::serialize::SerializationError> {
                self.serialize_with_flags(writer, crate::serialize::EmptyFlags)
            }

            #[inline]
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
                const BYTE_SIZE: usize = $byte_size;

                let (output_bit_size, output_byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                if F::len() > (output_bit_size - P::MODULUS_BITS as usize) {
                    return Err(crate::serialize::SerializationError::NotEnoughSpace);
                }

                let mut masked_bytes = [0; BYTE_SIZE];
                reader.read_exact(&mut masked_bytes[..output_byte_size])?;

                let flags = F::from_u8_remove_flags(&mut masked_bytes[output_byte_size - 1]);

                Ok((Self::read(&masked_bytes[..])?, flags))
            }
        }

        impl<P: $params> CanonicalDeserialize for $field<P> {
            #[allow(unused_qualifications)]
            fn deserialize<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                const BYTE_SIZE: usize = $byte_size;

                let (_, output_byte_size) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());

                let mut masked_bytes = [0; BYTE_SIZE];
                reader.read_exact(&mut masked_bytes[..output_byte_size])?;

                Ok(Self::read(&masked_bytes[..])?)
            }
        }
    };
}

macro_rules! impl_sw_curve_serializer {
    ($params: ident) => {
        impl<P: $params> CanonicalSerialize for GroupAffine<P> {
            #[allow(unused_qualifications)]
            #[inline]
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

            #[inline]
            fn serialized_size(&self) -> usize {
                self.x.serialized_size()
            }

            #[allow(unused_qualifications)]
            #[inline]
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

            #[inline]
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
            #[inline]
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

            #[inline]
            fn serialized_size(&self) -> usize {
                CanonicalSerialize::serialized_size(&self.x)
            }

            #[allow(unused_qualifications)]
            #[inline]
            fn serialize_uncompressed<W: crate::io::Write>(
                &self,
                writer: &mut W,
            ) -> Result<(), crate::serialize::SerializationError> {
                self.x.serialize_uncompressed(writer)?;
                self.y.serialize_uncompressed(writer)?;
                Ok(())
            }

            #[inline]
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
