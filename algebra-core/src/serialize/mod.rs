mod error;
mod flags;
use crate::io::{Read, Write};
pub use error::*;
pub use flags::*;

/// Serializer in little endian format.
pub trait CanonicalSerialize {
    /// Serializes `self` into `writer`.
    fn serialize<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize_with_flags(self, writer, EmptyFlags::default())
    }
    /// Serializes `self` and `flags` into `writer`.
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        writer: &mut W,
        flags: F,
    ) -> Result<(), SerializationError>;
    fn serialized_size(&self) -> usize;

    /// Serializes `self` into `writer` without compression.
    fn serialize_uncompressed<W: Write>(&self, writer: &mut W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(self, writer)
    }
    fn uncompressed_size(&self) -> usize {
        CanonicalSerialize::serialized_size(self)
    }
}

/// Deserializer in little endian format.
pub trait CanonicalDeserialize: Sized {
    /// Reads `Self` from `reader`.
    fn deserialize<R: Read>(reader: &mut R) -> Result<Self, SerializationError>;
    /// Reads `Self` and `Flags` from `reader`.
    /// Returns empty flags by default.
    fn deserialize_with_flags<R: Read, F: Flags>(
        reader: &mut R,
    ) -> Result<(Self, F), SerializationError>;

    /// Reads `Self` from `reader` without compression.
    fn deserialize_uncompressed<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        CanonicalDeserialize::deserialize(reader)
    }
}

impl CanonicalSerialize for u64 {
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        writer: &mut W,
        flags: F,
    ) -> Result<(), SerializationError> {
        // Encode flags into very first u64 at most significant bits.
        let value = *self | flags.u64_bitmask();
        // Then write the u64 in little endian order.
        writer.write_all(&value.to_le_bytes())?;
        Ok(())
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

    fn deserialize_with_flags<R: Read, F: Flags>(
        reader: &mut R,
    ) -> Result<(Self, F), SerializationError> {
        let mut bytes = [0u8; 8];
        reader.read_exact(&mut bytes)?;
        let mut value = u64::from_le_bytes(bytes);
        let flags = Flags::from_u64_remove_flags(&mut value);
        Ok((value, flags))
    }
}

impl CanonicalSerialize for bool {
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        writer: &mut W,
        _flags: F,
    ) -> Result<(), SerializationError> {
        let value = if *self { 1 } else { 0 };
        writer.write_all(&[value])?;
        Ok(())
    }

    fn serialized_size(&self) -> usize {
        1
    }
}

impl CanonicalDeserialize for bool {
    fn deserialize<R: Read>(reader: &mut R) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; 1];
        reader.read_exact(&mut bytes)?;
        Ok(bytes[0] != 0)
    }

    fn deserialize_with_flags<R: Read, F: Flags>(
        reader: &mut R,
    ) -> Result<(Self, F), SerializationError> {
        Ok((CanonicalDeserialize::deserialize(reader)?, F::default()))
    }
}

pub fn buffer_bit_byte_size(modulus_bits: usize) -> (usize, usize) {
    let byte_size = (modulus_bits + 7) / 8;
    ((byte_size * 8), byte_size)
}

macro_rules! impl_prime_field_serializer {
    ($field: ident, $params: ident, $byte_size: expr) => {
        impl<P: $params> CanonicalSerialize for $field<P> {
            #[allow(unused_qualifications)]
            fn serialize_with_flags<W: crate::io::Write, F: crate::serialize::Flags>(
                &self,
                writer: &mut W,
                flags: F,
            ) -> Result<(), crate::serialize::SerializationError> {
                let (output_bit_size, _) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                if F::len() > (output_bit_size - P::MODULUS_BITS as usize) {
                    return Err(crate::serialize::SerializationError::NotEnoughSpace);
                }

                CanonicalSerialize::serialize_with_flags(&self.into_repr(), writer, flags)
            }

            fn serialized_size(&self) -> usize {
                CanonicalSerialize::serialized_size(&self.into_repr())
            }
        }

        impl<P: $params> CanonicalDeserialize for $field<P> {
            #[allow(unused_qualifications)]
            fn deserialize<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let value: P::BigInt = CanonicalDeserialize::deserialize(reader)?;
                Ok(Self::from_repr(value))
            }

            #[allow(unused_qualifications)]
            fn deserialize_with_flags<R: crate::io::Read, F: crate::serialize::Flags>(
                reader: &mut R,
            ) -> Result<(Self, F), crate::serialize::SerializationError> {
                let (output_bit_size, _) =
                    crate::serialize::buffer_bit_byte_size($field::<P>::size_in_bits());
                if F::len() > (output_bit_size - P::MODULUS_BITS as usize) {
                    return Err(crate::serialize::SerializationError::NotEnoughSpace);
                }

                let (value, flags): (P::BigInt, _) =
                    CanonicalDeserialize::deserialize_with_flags(reader)?;
                Ok((Self::from_repr(value), flags))
            }
        }
    };
}

macro_rules! impl_sw_curve_serializer {
    ($params: ident) => {
        impl<P: $params> CanonicalSerialize for GroupAffine<P> {
            #[allow(unused_qualifications)]
            fn serialize_with_flags<W: crate::io::Write, F: crate::serialize::Flags>(
                &self,
                writer: &mut W,
                _flags: F,
            ) -> Result<(), crate::serialize::SerializationError> {
                if F::len() != 0 {
                    return Err(crate::serialize::SerializationError::UnexpectedFlags);
                }

                // We always ignore flags here.
                if self.is_zero() {
                    let flags = crate::serialize::SWFlags::infinity();
                    // Serialize 0.
                    CanonicalSerialize::serialize_with_flags(&P::BaseField::zero(), writer, flags)
                } else {
                    let flags = crate::serialize::SWFlags::y_sign(self.y > -self.y);
                    CanonicalSerialize::serialize_with_flags(&self.x, writer, flags)
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
                CanonicalSerialize::serialize_uncompressed(&self.x, writer)?;
                CanonicalSerialize::serialize_uncompressed(&self.y, writer)?;
                CanonicalSerialize::serialize_uncompressed(&self.infinity, writer)?;
                Ok(())
            }

            fn uncompressed_size(&self) -> usize {
                CanonicalSerialize::uncompressed_size(&self.x)
                    + CanonicalSerialize::uncompressed_size(&self.y)
                    + CanonicalSerialize::uncompressed_size(&self.infinity)
            }
        }

        impl<P: $params> CanonicalDeserialize for GroupAffine<P> {
            #[allow(unused_qualifications)]
            fn deserialize<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let (point, _): (Self, crate::serialize::EmptyFlags) =
                    CanonicalDeserialize::deserialize_with_flags(reader)?;
                Ok(point)
            }

            #[allow(unused_qualifications)]
            fn deserialize_with_flags<R: crate::io::Read, F: crate::serialize::Flags>(
                reader: &mut R,
            ) -> Result<(Self, F), crate::serialize::SerializationError> {
                let (x, flags): (P::BaseField, crate::serialize::SWFlags) =
                    CanonicalDeserialize::deserialize_with_flags(reader)?;
                if flags.is_infinity {
                    Ok((Self::zero(), F::default()))
                } else {
                    let p = GroupAffine::<P>::get_point_from_x(x, flags.y_sign)
                        .ok_or(crate::serialize::SerializationError::InvalidData)?;
                    if !p.is_in_correct_subgroup_assuming_on_curve() {
                        return Err(crate::serialize::SerializationError::InvalidData);
                    }
                    Ok((p, F::default()))
                }
            }

            #[allow(unused_qualifications)]
            fn deserialize_uncompressed<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let x: P::BaseField = CanonicalDeserialize::deserialize(reader)?;
                let y: P::BaseField = CanonicalDeserialize::deserialize(reader)?;
                let infinity: bool = CanonicalDeserialize::deserialize(reader)?;

                let p = GroupAffine::<P>::new(x, y, infinity);
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
            fn serialize_with_flags<W: crate::io::Write, F: crate::serialize::Flags>(
                &self,
                writer: &mut W,
                _flags: F,
            ) -> Result<(), crate::serialize::SerializationError> {
                if F::len() != 0 {
                    return Err(crate::serialize::SerializationError::UnexpectedFlags);
                }

                // We always ignore flags here.
                if self.is_zero() {
                    let flags = crate::serialize::EdwardsFlags::default();
                    // Serialize 0.
                    CanonicalSerialize::serialize_with_flags(&P::BaseField::zero(), writer, flags)
                } else {
                    let flags = crate::serialize::EdwardsFlags::y_sign(self.y > -self.y);
                    CanonicalSerialize::serialize_with_flags(&self.x, writer, flags)
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
                CanonicalSerialize::serialize_uncompressed(&self.x, writer)?;
                CanonicalSerialize::serialize_uncompressed(&self.y, writer)?;
                Ok(())
            }

            fn uncompressed_size(&self) -> usize {
                CanonicalSerialize::uncompressed_size(&self.x)
                    + CanonicalSerialize::uncompressed_size(&self.y)
            }
        }

        impl<P: $params> CanonicalDeserialize for GroupAffine<P> {
            #[allow(unused_qualifications)]
            fn deserialize<R: crate::io::Read>(
                reader: &mut R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let (point, _): (Self, crate::serialize::EmptyFlags) =
                    CanonicalDeserialize::deserialize_with_flags(reader)?;
                Ok(point)
            }

            #[allow(unused_qualifications)]
            fn deserialize_with_flags<R: crate::io::Read, F: crate::serialize::Flags>(
                reader: &mut R,
            ) -> Result<(Self, F), crate::serialize::SerializationError> {
                let (x, flags): (P::BaseField, crate::serialize::EdwardsFlags) =
                    CanonicalDeserialize::deserialize_with_flags(reader)?;
                if x == P::BaseField::zero() {
                    Ok((Self::zero(), F::default()))
                } else {
                    let p = GroupAffine::<P>::get_point_from_x(x, flags.y_sign)
                        .ok_or(crate::serialize::SerializationError::InvalidData)?;
                    if !p.is_in_correct_subgroup_assuming_on_curve() {
                        return Err(crate::serialize::SerializationError::InvalidData);
                    }
                    Ok((p, F::default()))
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

        let a = true;
        let mut serialized = vec![0u8; a.serialized_size()];
        let mut cursor = Cursor::new(&mut serialized[..]);
        a.serialize(&mut cursor).unwrap();

        let mut cursor = Cursor::new(&serialized[..]);
        let b = bool::deserialize(&mut cursor).unwrap();
        assert_eq!(a, b);
    }
}
