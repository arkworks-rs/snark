mod error;
mod flags;
pub use crate::{
    bytes::{FromBytes, ToBytes},
    io::{Read, Write},
};
pub use error::*;
pub use flags::*;

#[cfg(feature = "derive")]
#[doc(hidden)]
pub use algebra_core_derive::*;

use crate::{BTreeMap, Cow, ToOwned, Vec};
use core::convert::TryFrom;

/// Serializer in little endian format allowing to encode flags.
pub trait CanonicalSerializeWithFlags: CanonicalSerialize {
    /// Serializes `self` and `flags` into `writer`.
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        writer: W,
        flags: F,
    ) -> Result<(), SerializationError>;
}

/// Helper trait to get serialized size for constant sized structs.
pub trait ConstantSerializedSize: CanonicalSerialize {
    const SERIALIZED_SIZE: usize;
    const UNCOMPRESSED_SIZE: usize;
}

/// Serializer in little endian format.
/// This trait can be derived if all fields of a struct implement
/// `CanonicalSerialize` and the `derive` feature is enabled.
///
/// # Example
/// ```
/// // The `derive` feature must be set for the derivation to work.
/// use algebra_core::serialize::*;
///
/// # #[cfg(feature = "derive")]
/// #[derive(CanonicalSerialize)]
/// struct TestStruct {
///     a: u64,
///     b: (u64, (u64, u64)),
/// }
/// ```
///
/// If your code depends on `algebra` instead, the example works analogously
/// when importing `algebra::serialize::*`.
pub trait CanonicalSerialize {
    /// Serializes `self` into `writer`.
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError>;

    fn serialized_size(&self) -> usize;

    /// Serializes `self` into `writer` without compression.
    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.serialize(writer)
    }

    /// Serializes `self` into `writer` without compression, and without
    /// performing validity checks. Should be used *only* when there is no
    /// danger of adversarial manipulation of the output.
    #[inline]
    fn serialize_unchecked<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.serialize_uncompressed(writer)
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
        reader: R,
    ) -> Result<(Self, F), SerializationError>;
}

/// Deserializer in little endian format.
/// This trait can be derived if all fields of a struct implement
/// `CanonicalDeserialize` and the `derive` feature is enabled.
///
/// # Example
/// ```
/// // The `derive` feature must be set for the derivation to work.
/// use algebra_core::serialize::*;
///
/// # #[cfg(feature = "derive")]
/// #[derive(CanonicalDeserialize)]
/// struct TestStruct {
///     a: u64,
///     b: (u64, (u64, u64)),
/// }
/// ```
///
/// If your code depends on `algebra` instead, the example works analogously
/// when importing `algebra::serialize::*`.
pub trait CanonicalDeserialize: Sized {
    /// Reads `Self` from `reader`.
    fn deserialize<R: Read>(reader: R) -> Result<Self, SerializationError>;

    /// Reads `Self` from `reader` without compression.
    #[inline]
    fn deserialize_uncompressed<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Self::deserialize(reader)
    }

    /// Reads `self` from `reader` without compression, and without performing
    /// validity checks. Should be used *only* when the input is trusted.
    #[inline]
    fn deserialize_unchecked<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Self::deserialize_uncompressed(reader)
    }
}

macro_rules! impl_uint {
    ($ty: ident) => {
        impl CanonicalSerialize for $ty {
            #[inline]
            fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
                Ok(writer.write_all(&self.to_le_bytes())?)
            }

            #[inline]
            fn serialized_size(&self) -> usize {
                Self::SERIALIZED_SIZE
            }
        }

        impl ConstantSerializedSize for $ty {
            const SERIALIZED_SIZE: usize = core::mem::size_of::<$ty>();
            const UNCOMPRESSED_SIZE: usize = Self::SERIALIZED_SIZE;
        }

        impl CanonicalDeserialize for $ty {
            #[inline]
            fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
                let mut bytes = [0u8; Self::SERIALIZED_SIZE];
                reader.read_exact(&mut bytes)?;
                Ok($ty::from_le_bytes(bytes))
            }
        }
    };
}

impl_uint!(u8);
impl_uint!(u16);
impl_uint!(u32);
impl_uint!(u64);

impl CanonicalSerialize for usize {
    #[inline]
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        Ok(writer.write_all(&(*self as u64).to_le_bytes())?)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        Self::SERIALIZED_SIZE
    }
}

impl ConstantSerializedSize for usize {
    const SERIALIZED_SIZE: usize = core::mem::size_of::<u64>();
    const UNCOMPRESSED_SIZE: usize = Self::SERIALIZED_SIZE;
}

impl CanonicalDeserialize for usize {
    #[inline]
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; Self::SERIALIZED_SIZE];
        reader.read_exact(&mut bytes)?;
        usize::try_from(u64::from_le_bytes(bytes)).map_err(|_| SerializationError::InvalidData)
    }
}

impl<T: CanonicalSerialize> CanonicalSerialize for [T] {
    #[inline]
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize(&mut writer)?;
        for item in self.iter() {
            item.serialize(&mut writer)?;
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
    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize(&mut writer)?;
        for item in self.iter() {
            item.serialize_uncompressed(&mut writer)?;
        }
        Ok(())
    }

    #[inline]
    fn serialize_unchecked<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize(&mut writer)?;
        for item in self.iter() {
            item.serialize_unchecked(&mut writer)?;
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

impl<T: CanonicalSerialize> CanonicalSerialize for Vec<T> {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.as_slice().serialize(writer)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        self.as_slice().serialized_size()
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.as_slice().serialize_uncompressed(writer)
    }

    #[inline]
    fn serialize_unchecked<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.as_slice().serialize_unchecked(writer)
    }

    #[inline]
    fn uncompressed_size(&self) -> usize {
        self.as_slice().uncompressed_size()
    }
}

impl<T: CanonicalDeserialize> CanonicalDeserialize for Vec<T> {
    #[inline]
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = u64::deserialize(&mut reader)?;
        let mut values = vec![];
        for _ in 0..len {
            values.push(T::deserialize(&mut reader)?);
        }
        Ok(values)
    }

    #[inline]
    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = u64::deserialize(&mut reader)?;
        let mut values = vec![];
        for _ in 0..len {
            values.push(T::deserialize_uncompressed(&mut reader)?);
        }
        Ok(values)
    }

    #[inline]
    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = u64::deserialize(&mut reader)?;
        let mut values = vec![];
        for _ in 0..len {
            values.push(T::deserialize_unchecked(&mut reader)?);
        }
        Ok(values)
    }
}

#[inline]
pub fn buffer_bit_byte_size(modulus_bits: usize) -> (usize, usize) {
    let byte_size = buffer_byte_size(modulus_bits);
    ((byte_size * 8), byte_size)
}

#[inline]
pub const fn buffer_byte_size(modulus_bits: usize) -> usize {
    (modulus_bits + 7) / 8
}

macro_rules! impl_prime_field_serializer {
    ($field: ident, $params: ident, $byte_size: expr) => {
        impl<P: $params> CanonicalSerializeWithFlags for $field<P> {
            #[allow(unused_qualifications)]
            fn serialize_with_flags<W: crate::io::Write, F: crate::serialize::Flags>(
                &self,
                mut writer: W,
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

        impl<P: $params> ConstantSerializedSize for $field<P> {
            const SERIALIZED_SIZE: usize = crate::serialize::buffer_byte_size(
                <$field<P> as crate::PrimeField>::Params::MODULUS_BITS as usize,
            );
            const UNCOMPRESSED_SIZE: usize = Self::SERIALIZED_SIZE;
        }

        impl<P: $params> CanonicalSerialize for $field<P> {
            #[allow(unused_qualifications)]
            #[inline]
            fn serialize<W: crate::io::Write>(
                &self,
                writer: W,
            ) -> Result<(), crate::serialize::SerializationError> {
                self.serialize_with_flags(writer, crate::serialize::EmptyFlags)
            }

            #[inline]
            fn serialized_size(&self) -> usize {
                Self::SERIALIZED_SIZE
            }
        }

        impl<P: $params> CanonicalDeserializeWithFlags for $field<P> {
            #[allow(unused_qualifications)]
            fn deserialize_with_flags<R: crate::io::Read, F: crate::serialize::Flags>(
                mut reader: R,
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
                mut reader: R,
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
                writer: W,
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
                Self::SERIALIZED_SIZE
            }

            #[allow(unused_qualifications)]
            #[inline]
            fn serialize_uncompressed<W: crate::io::Write>(
                &self,
                mut writer: W,
            ) -> Result<(), crate::serialize::SerializationError> {
                let flags = if self.is_zero() {
                    crate::serialize::SWFlags::infinity()
                } else {
                    crate::serialize::SWFlags::default()
                };
                self.x.serialize(&mut writer)?;
                self.y.serialize_with_flags(&mut writer, flags)?;
                Ok(())
            }

            #[inline]
            fn uncompressed_size(&self) -> usize {
                Self::UNCOMPRESSED_SIZE
            }
        }

        impl<P: $params> ConstantSerializedSize for GroupAffine<P> {
            const SERIALIZED_SIZE: usize =
                <P::BaseField as ConstantSerializedSize>::SERIALIZED_SIZE;
            const UNCOMPRESSED_SIZE: usize =
                2 * <P::BaseField as ConstantSerializedSize>::SERIALIZED_SIZE;
        }

        impl<P: $params> CanonicalDeserialize for GroupAffine<P> {
            #[allow(unused_qualifications)]
            fn deserialize<R: crate::io::Read>(
                reader: R,
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
                reader: R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let p = Self::deserialize_unchecked(reader)?;

                if !p.is_in_correct_subgroup_assuming_on_curve() {
                    return Err(crate::serialize::SerializationError::InvalidData);
                }
                Ok(p)
            }

            #[allow(unused_qualifications)]
            fn deserialize_unchecked<R: crate::io::Read>(
                mut reader: R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let x: P::BaseField = CanonicalDeserialize::deserialize(&mut reader)?;
                let (y, flags): (P::BaseField, crate::serialize::SWFlags) =
                    CanonicalDeserializeWithFlags::deserialize_with_flags(&mut reader)?;
                let p = GroupAffine::<P>::new(x, y, flags.is_infinity());
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
                writer: W,
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
                Self::SERIALIZED_SIZE
            }

            #[allow(unused_qualifications)]
            #[inline]
            fn serialize_uncompressed<W: crate::io::Write>(
                &self,
                mut writer: W,
            ) -> Result<(), crate::serialize::SerializationError> {
                self.x.serialize_uncompressed(&mut writer)?;
                self.y.serialize_uncompressed(&mut writer)?;
                Ok(())
            }

            #[inline]
            fn uncompressed_size(&self) -> usize {
                Self::UNCOMPRESSED_SIZE
            }
        }

        impl<P: $params> ConstantSerializedSize for GroupAffine<P> {
            const SERIALIZED_SIZE: usize =
                <P::BaseField as ConstantSerializedSize>::SERIALIZED_SIZE;
            const UNCOMPRESSED_SIZE: usize =
                2 * <P::BaseField as ConstantSerializedSize>::SERIALIZED_SIZE;
        }

        impl<P: $params> CanonicalDeserialize for GroupAffine<P> {
            #[allow(unused_qualifications)]
            fn deserialize<R: crate::io::Read>(
                mut reader: R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let (x, flags): (P::BaseField, crate::serialize::EdwardsFlags) =
                    CanonicalDeserializeWithFlags::deserialize_with_flags(&mut reader)?;
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
                reader: R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let p = Self::deserialize_unchecked(reader)?;

                if !p.is_in_correct_subgroup_assuming_on_curve() {
                    return Err(crate::serialize::SerializationError::InvalidData);
                }
                Ok(p)
            }

            #[allow(unused_qualifications)]
            fn deserialize_unchecked<R: crate::io::Read>(
                mut reader: R,
            ) -> Result<Self, crate::serialize::SerializationError> {
                let x: P::BaseField = CanonicalDeserialize::deserialize(&mut reader)?;
                let y: P::BaseField = CanonicalDeserialize::deserialize(&mut reader)?;

                let p = GroupAffine::<P>::new(x, y);
                Ok(p)
            }
        }
    };
}

// Implement Serialization for tuples
macro_rules! impl_tuple {
    ($( $ty: ident : $no: tt, )+) => {
        impl<$($ty, )+> CanonicalSerialize for ($($ty,)+) where
            $($ty: CanonicalSerialize,)+
        {
            #[inline]
            fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
                $(self.$no.serialize(&mut writer)?;)*
                Ok(())
            }

            #[inline]
            fn serialized_size(&self) -> usize {
                [$(
                    self.$no.serialized_size(),
                )*].iter().sum()
            }

            #[inline]
            fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
                $(self.$no.serialize_uncompressed(&mut writer)?;)*
                Ok(())
            }

            #[inline]
            fn serialize_unchecked<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
                $(self.$no.serialize_unchecked(&mut writer)?;)*
                Ok(())
            }

            #[inline]
            fn uncompressed_size(&self) -> usize {
                [$(
                    self.$no.uncompressed_size(),
                )*].iter().sum()
            }
        }

        impl<$($ty, )+> CanonicalDeserialize for ($($ty,)+) where
            $($ty: CanonicalDeserialize,)+
        {
            #[inline]
            fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
                Ok(($(
                    $ty::deserialize(&mut reader)?,
                )+))
            }

            #[inline]
            fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
                Ok(($(
                    $ty::deserialize_uncompressed(&mut reader)?,
                )+))
            }

            #[inline]
            fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
                Ok(($(
                    $ty::deserialize_unchecked(&mut reader)?,
                )+))
            }
        }
    }
}

impl_tuple!(A:0, B:1,);
impl_tuple!(A:0, B:1, C:2,);
impl_tuple!(A:0, B:1, C:2, D:3,);

// No-op
impl<T> CanonicalSerialize for core::marker::PhantomData<T> {
    #[inline]
    fn serialize<W: Write>(&self, _writer: W) -> Result<(), SerializationError> {
        Ok(())
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        0
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, _writer: W) -> Result<(), SerializationError> {
        Ok(())
    }
}

impl<T> CanonicalDeserialize for core::marker::PhantomData<T> {
    #[inline]
    fn deserialize<R: Read>(_reader: R) -> Result<Self, SerializationError> {
        Ok(core::marker::PhantomData)
    }

    #[inline]
    fn deserialize_uncompressed<R: Read>(_reader: R) -> Result<Self, SerializationError> {
        Ok(core::marker::PhantomData)
    }
}

impl<'a, T: CanonicalSerialize + ToOwned> CanonicalSerialize for Cow<'a, T> {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.as_ref().serialize(writer)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        self.as_ref().serialized_size()
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.as_ref().serialize_uncompressed(writer)
    }
}

impl<'a, T> CanonicalDeserialize for Cow<'a, T>
where
    T: ToOwned,
    <T as ToOwned>::Owned: CanonicalDeserialize,
{
    #[inline]
    fn deserialize<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Ok(Cow::Owned(<T as ToOwned>::Owned::deserialize(reader)?))
    }

    #[inline]
    fn deserialize_uncompressed<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Ok(Cow::Owned(<T as ToOwned>::Owned::deserialize_uncompressed(
            reader,
        )?))
    }
}

impl<T: CanonicalSerialize> CanonicalSerialize for Option<T> {
    #[inline]
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.is_some().serialize(&mut writer)?;
        if let Some(item) = self {
            item.serialize(&mut writer)?;
        }

        Ok(())
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        self.is_some().serialized_size()
            + if let Some(item) = self {
                item.serialized_size()
            } else {
                0
            }
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.is_some().serialize_uncompressed(&mut writer)?;
        if let Some(item) = self {
            item.serialize_uncompressed(&mut writer)?;
        }

        Ok(())
    }
}

impl<T: CanonicalDeserialize> CanonicalDeserialize for Option<T> {
    #[inline]
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let is_some = bool::deserialize(&mut reader)?;
        let data = if is_some {
            Some(T::deserialize(&mut reader)?)
        } else {
            None
        };

        Ok(data)
    }

    #[inline]
    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let is_some = bool::deserialize(&mut reader)?;
        let data = if is_some {
            Some(T::deserialize_uncompressed(&mut reader)?)
        } else {
            None
        };

        Ok(data)
    }
}

impl CanonicalSerialize for bool {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        Ok(self.write(writer)?)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        1
    }
}

impl CanonicalDeserialize for bool {
    #[inline]
    fn deserialize<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Ok(bool::read(reader)?)
    }
}

impl<K, V> CanonicalSerialize for BTreeMap<K, V>
where
    K: CanonicalSerialize,
    V: CanonicalSerialize,
{
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize(&mut writer)?;
        for (k, v) in self.iter() {
            k.serialize(&mut writer)?;
            v.serialize(&mut writer)?;
        }
        Ok(())
    }

    fn serialized_size(&self) -> usize {
        8 + self
            .iter()
            .map(|(k, v)| k.serialized_size() + v.serialized_size())
            .sum::<usize>()
    }
}

impl<K, V> CanonicalDeserialize for BTreeMap<K, V>
where
    K: Ord + CanonicalDeserialize,
    V: CanonicalDeserialize,
{
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = u64::deserialize(&mut reader)?;
        let mut map = BTreeMap::new();
        for _ in 0..len {
            map.insert(K::deserialize(&mut reader)?, V::deserialize(&mut reader)?);
        }
        Ok(map)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn test_serialize<
        T: PartialEq + core::fmt::Debug + CanonicalSerialize + CanonicalDeserialize,
    >(
        data: T,
    ) {
        let mut serialized = vec![0; data.serialized_size()];
        data.serialize(&mut serialized[..]).unwrap();
        let de = T::deserialize(&serialized[..]).unwrap();
        assert_eq!(data, de);

        let mut serialized = vec![0; data.uncompressed_size()];
        data.serialize_uncompressed(&mut serialized[..]).unwrap();
        let de = T::deserialize_uncompressed(&serialized[..]).unwrap();
        assert_eq!(data, de);

        let mut serialized = vec![0; data.uncompressed_size()];
        data.serialize_unchecked(&mut serialized[..]).unwrap();
        let de = T::deserialize_unchecked(&serialized[..]).unwrap();
        assert_eq!(data, de);
    }

    #[test]
    fn test_vec() {
        test_serialize(vec![1u64, 2, 3, 4, 5]);
        test_serialize(Vec::<u64>::new());
    }

    #[test]
    fn test_uint() {
        test_serialize(192830918usize);
        test_serialize(192830918u64);
        test_serialize(192830918u32);
        test_serialize(22313u16);
        test_serialize(123u8);
    }

    #[test]
    fn test_tuple() {
        test_serialize((123u64, 234u32, 999u16));
    }

    #[test]
    fn test_tuple_vec() {
        test_serialize(vec![
            (123u64, 234u32, 999u16),
            (123u64, 234u32, 999u16),
            (123u64, 234u32, 999u16),
        ]);
    }

    #[test]
    fn test_option() {
        test_serialize(Some(3u32));
        test_serialize(None::<u32>);
    }

    #[test]
    fn test_bool() {
        test_serialize(true);
        test_serialize(false);
    }

    #[test]
    fn test_btreemap() {
        let mut map = BTreeMap::new();
        map.insert(0u64, true);
        map.insert(5u64, false);
        test_serialize(map);
        let mut map = BTreeMap::new();
        map.insert(10u64, vec![1u8, 2u8, 3u8]);
        map.insert(50u64, vec![4u8, 5u8, 6u8]);
        test_serialize(map);
    }

    #[test]
    fn test_phantomdata() {
        test_serialize(core::marker::PhantomData::<u64>);
    }
}
