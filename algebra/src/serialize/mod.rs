mod error;
mod flags;

pub use error::*;
pub use flags::*;
pub use std::io::{Read, Write};
use std::{
    borrow::{Cow, ToOwned},
    collections::{BTreeMap, BTreeSet},
    convert::TryFrom,
    rc::Rc,
    string::String,
    vec::Vec,
};

#[cfg(feature = "derive")]
#[doc(hidden)]
pub use algebra_derive::*;

/// Serializer in little endian format allowing to encode flags.
pub trait CanonicalSerializeWithFlags: CanonicalSerialize {
    /// Serializes `self` and `flags` into `writer`.
    fn serialize_with_flags<W: Write, F: Flags>(
        &self,
        writer: W,
        flags: F,
    ) -> Result<(), SerializationError>;

    /// Get size of serialized self with flags.
    fn serialized_size_with_flags<F: Flags>(&self) -> usize;
}

/// Serializer in little endian format.
/// The serialization format must be 'length-extension' safe.
/// e.g. if T implements Canonical Serialize and Deserialize,
/// then for all strings `x, y`, if `a = T::deserialize(Reader(x))` and `a` is not an error,
/// then it must be the case that `a = T::deserialize(Reader(x || y))`,
/// and that both readers read the same number of bytes.
///
/// This trait can be derived if all fields of a struct implement
/// `CanonicalSerialize` and the `derive` feature is enabled.
pub trait CanonicalSerialize {
    /// Serializes `self` into `writer`.
    /// It is left up to a particular type for how it strikes the
    /// serialization efficiency vs compression tradeoff.
    /// For standard types (e.g. `bool`, lengths, etc.) typically an uncompressed
    /// form is used, whereas for algebraic types compressed forms are used.
    ///
    /// Particular examples of interest:
    /// `bool` - 1 byte encoding
    /// uints - Direct encoding
    /// Length prefixing (for any container implemented by default) - 8 byte encoding
    /// Elliptic curves - compressed point encoding
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError>;

    fn serialized_size(&self) -> usize;

    /// Like `serialize_uncompressed()`, but doesn't write (if present) any additional information
    /// required to reconstruct `self` (e.g. the length of a container type).
    /// For this reason, there isn't any deserialization counterpart function in
    /// CanonicalDeserialize trait.
    fn serialize_without_metadata<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize_uncompressed(self, writer)
    }

    /// Serializes `self` into `writer` without compression.
    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(self, writer)
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
pub trait CanonicalDeserialize: Sized {
    /// Reads `Self` from `reader`.
    fn deserialize<R: Read>(reader: R) -> Result<Self, SerializationError>;

    /// Reads `Self` from `reader` without performing validity checks.
    /// Should be used *only* when the input is trusted.
    fn deserialize_unchecked<R: Read>(reader: R) -> Result<Self, SerializationError> {
        CanonicalDeserialize::deserialize(reader)
    }

    /// Reads `Self` from `reader` without compression.
    #[inline]
    fn deserialize_uncompressed<R: Read>(reader: R) -> Result<Self, SerializationError> {
        CanonicalDeserialize::deserialize(reader)
    }

    /// Reads `self` from `reader` without compression, and without performing
    /// validity checks. Should be used *only* when the input is trusted.
    #[inline]
    fn deserialize_uncompressed_unchecked<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Self::deserialize_uncompressed(reader)
    }
}

// Macro for implementing serialize for u8, u16, u32, u64
macro_rules! impl_uint {
    ($ty: ident) => {
        impl CanonicalSerialize for $ty {
            #[inline]
            fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
                Ok(writer.write_all(&self.to_le_bytes())?)
            }

            #[inline]
            fn serialized_size(&self) -> usize {
                std::mem::size_of::<$ty>()
            }
        }

        impl CanonicalDeserialize for $ty {
            #[inline]
            fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
                let mut bytes = [0u8; std::mem::size_of::<$ty>()];
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
impl_uint!(u128);

// Serialize usize with 8 bytes
impl CanonicalSerialize for usize {
    #[inline]
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        Ok(writer.write_all(&(*self as u64).to_le_bytes())?)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        std::mem::size_of::<u64>()
    }
}

impl CanonicalDeserialize for usize {
    #[inline]
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let mut bytes = [0u8; std::mem::size_of::<u64>()];
        reader.read_exact(&mut bytes)?;
        usize::try_from(u64::from_le_bytes(bytes)).map_err(|_| SerializationError::InvalidData)
    }
}

impl<'a, T: 'a + CanonicalSerialize> CanonicalSerialize for &'a T {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(*self, writer)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        (*self).serialized_size()
    }

    #[inline]
    fn serialize_without_metadata<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize_without_metadata(*self, writer)
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize_uncompressed(*self, writer)
    }

    #[inline]
    fn uncompressed_size(&self) -> usize {
        (*self).uncompressed_size()
    }
}

// Implement Serialization for `String`
// It is serialized by obtaining its byte representation as a Vec<u8> and
// serializing that. This yields an end serialization of
// `string.len() || string_bytes`.
impl CanonicalSerialize for String {
    #[inline]
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(&self.clone().into_bytes(), &mut writer)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        self.clone().into_bytes().serialized_size()
    }
}

impl CanonicalDeserialize for String {
    #[inline]
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        String::from_utf8(<Vec<u8> as CanonicalDeserialize>::deserialize(&mut reader)?)
            .map_err(|_| SerializationError::InvalidData)
    }
}

impl<T: CanonicalSerialize> CanonicalSerialize for [T] {
    #[inline]
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        CanonicalSerialize::serialize(&len, &mut writer)?;
        for item in self.iter() {
            CanonicalSerialize::serialize(item, &mut writer)?;
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
    fn serialize_without_metadata<W: Write>(
        &self,
        mut writer: W,
    ) -> Result<(), SerializationError> {
        for item in self.iter() {
            item.serialize_without_metadata(&mut writer)?;
        }
        Ok(())
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        CanonicalSerialize::serialize(&len, &mut writer)?;
        for item in self.iter() {
            item.serialize_uncompressed(&mut writer)?;
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

impl<'a, T: 'a + CanonicalSerialize> CanonicalSerialize for &'a [T] {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(*self, writer)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        (*self).serialized_size()
    }

    #[inline]
    fn serialize_without_metadata<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize_without_metadata(*self, writer)
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize_uncompressed(*self, writer)
    }

    #[inline]
    fn uncompressed_size(&self) -> usize {
        (*self).uncompressed_size()
    }
}

impl<T: CanonicalSerialize> CanonicalSerialize for Vec<T> {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(self.as_slice(), writer)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        self.as_slice().serialized_size()
    }

    #[inline]
    fn serialize_without_metadata<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.as_slice().serialize_without_metadata(writer)
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.as_slice().serialize_uncompressed(writer)
    }

    #[inline]
    fn uncompressed_size(&self) -> usize {
        self.as_slice().uncompressed_size()
    }
}

impl<T: CanonicalDeserialize> CanonicalDeserialize for Vec<T> {
    #[inline]
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = <u64 as CanonicalDeserialize>::deserialize(&mut reader)?;
        let mut values = Vec::new();
        for _ in 0..len {
            values.push(T::deserialize(&mut reader)?);
        }
        Ok(values)
    }

    #[inline]
    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = <u64 as CanonicalDeserialize>::deserialize(&mut reader)?;
        let mut values = Vec::new();
        for _ in 0..len {
            values.push(T::deserialize_unchecked(&mut reader)?);
        }
        Ok(values)
    }

    #[inline]
    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = <u64 as CanonicalDeserialize>::deserialize(&mut reader)?;
        let mut values = Vec::new();
        for _ in 0..len {
            values.push(T::deserialize_uncompressed(&mut reader)?);
        }
        Ok(values)
    }

    #[inline]
    fn deserialize_uncompressed_unchecked<R: Read>(
        mut reader: R,
    ) -> Result<Self, SerializationError> {
        let len = <u64 as CanonicalDeserialize>::deserialize(&mut reader)?;
        let mut values = Vec::new();
        for _ in 0..len {
            values.push(T::deserialize_uncompressed_unchecked(&mut reader)?);
        }
        Ok(values)
    }
}

#[inline]
pub fn buffer_bit_byte_size(modulus_bits: usize) -> (usize, usize) {
    let byte_size = buffer_byte_size(modulus_bits);
    ((byte_size * 8), byte_size)
}

/// Converts the number of bits required to represent a number
/// into the number of bytes required to represent it.
#[inline]
pub const fn buffer_byte_size(modulus_bits: usize) -> usize {
    (modulus_bits + 7) / 8
}

// Implement Serialization for tuples
macro_rules! impl_tuple {
    ($( $ty: ident : $no: tt, )*) => {
        impl<$($ty, )*> CanonicalSerialize for ($($ty,)*) where
            $($ty: CanonicalSerialize,)*
        {
            #[inline]
            fn serialize<W: Write>(&self, mut _writer: W) -> Result<(), SerializationError> {
                $(CanonicalSerialize::serialize(&self.$no, &mut _writer)?;)*
                Ok(())
            }

            #[inline]
            fn serialized_size(&self) -> usize {
                [$(
                    self.$no.serialized_size(),
                )*].iter().sum()
            }

            #[inline]
            fn serialize_uncompressed<W: Write>(&self, mut _writer: W) -> Result<(), SerializationError> {
                $(self.$no.serialize_uncompressed(&mut _writer)?;)*
                Ok(())
            }

            #[inline]
            fn uncompressed_size(&self) -> usize {
                [$(
                    self.$no.uncompressed_size(),
                )*].iter().sum()
            }
        }

        impl<$($ty, )*> CanonicalDeserialize for ($($ty,)*) where
            $($ty: CanonicalDeserialize,)*
        {
            #[inline]
            fn deserialize<R: Read>(mut _reader: R) -> Result<Self, SerializationError> {
                Ok(($(
                    $ty::deserialize(&mut _reader)?,
                )*))
            }

            #[inline]
            fn deserialize_unchecked<R: Read>(mut _reader: R) -> Result<Self, SerializationError> {
                Ok(($(
                    $ty::deserialize_unchecked(&mut _reader)?,
                )*))
            }

            #[inline]
            fn deserialize_uncompressed<R: Read>(mut _reader: R) -> Result<Self, SerializationError> {
                Ok(($(
                    $ty::deserialize_uncompressed(&mut _reader)?,
                )*))
            }

            #[inline]
            fn deserialize_uncompressed_unchecked<R: Read>(mut _reader: R) -> Result<Self, SerializationError> {
                Ok(($(
                    $ty::deserialize_uncompressed_unchecked(&mut _reader)?,
                )*))
            }
        }
    }
}

impl_tuple!();
impl_tuple!(A:0, B:1,);
impl_tuple!(A:0, B:1, C:2,);
impl_tuple!(A:0, B:1, C:2, D:3,);

// No-op
impl<T> CanonicalSerialize for std::marker::PhantomData<T> {
    #[inline]
    fn serialize<W: Write>(&self, _writer: W) -> Result<(), SerializationError> {
        Ok(())
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        0
    }
}

impl<T> CanonicalDeserialize for std::marker::PhantomData<T> {
    #[inline]
    fn deserialize<R: Read>(_reader: R) -> Result<Self, SerializationError> {
        Ok(std::marker::PhantomData)
    }
}

// Serialize cow objects by serializing the underlying object.
impl<'a, T: CanonicalSerialize + ToOwned> CanonicalSerialize for Cow<'a, T> {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(self.as_ref(), writer)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        self.as_ref().serialized_size()
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        self.as_ref().serialize_uncompressed(writer)
    }

    fn uncompressed_size(&self) -> usize {
        self.as_ref().uncompressed_size()
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
    fn deserialize_unchecked<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Ok(Cow::Owned(<T as ToOwned>::Owned::deserialize_unchecked(
            reader,
        )?))
    }

    #[inline]
    fn deserialize_uncompressed<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Ok(Cow::Owned(<T as ToOwned>::Owned::deserialize_uncompressed(
            reader,
        )?))
    }

    #[inline]
    fn deserialize_uncompressed_unchecked<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Ok(Cow::Owned(
            <T as ToOwned>::Owned::deserialize_uncompressed_unchecked(reader)?,
        ))
    }
}

// If Option<T> is None, serialize as serialize(False).
// If its Some, serialize as serialize(True) || serialize(T)
impl<T: CanonicalSerialize> CanonicalSerialize for Option<T> {
    #[inline]
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(&self.is_some(), &mut writer)?;
        if let Some(item) = self {
            CanonicalSerialize::serialize(item, &mut writer)?;
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
    fn serialize_without_metadata<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        if self.is_some() {
            self.serialize_without_metadata(writer)?;
        }
        Ok(())
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.is_some().serialize_uncompressed(&mut writer)?;
        if let Some(item) = self {
            item.serialize_uncompressed(&mut writer)?;
        }

        Ok(())
    }

    #[inline]
    fn uncompressed_size(&self) -> usize {
        self.is_some().uncompressed_size()
            + if let Some(item) = self {
                item.uncompressed_size()
            } else {
                0
            }
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
    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let is_some = bool::deserialize_unchecked(&mut reader)?;
        let data = if is_some {
            Some(T::deserialize_unchecked(&mut reader)?)
        } else {
            None
        };

        Ok(data)
    }

    #[inline]
    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let is_some = bool::deserialize_uncompressed(&mut reader)?;
        let data = if is_some {
            Some(T::deserialize_uncompressed(&mut reader)?)
        } else {
            None
        };

        Ok(data)
    }

    #[inline]
    fn deserialize_uncompressed_unchecked<R: Read>(
        mut reader: R,
    ) -> Result<Self, SerializationError> {
        let is_some = bool::deserialize_uncompressed_unchecked(&mut reader)?;
        let data = if is_some {
            Some(T::deserialize_uncompressed_unchecked(&mut reader)?)
        } else {
            None
        };

        Ok(data)
    }
}

// Implement Serialization for `Rc<T>`
impl<T: CanonicalSerialize> CanonicalSerialize for Rc<T> {
    #[inline]
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        CanonicalSerialize::serialize(self.as_ref(), &mut writer)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        self.as_ref().serialized_size()
    }

    #[inline]
    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        self.as_ref().serialize_uncompressed(&mut writer)
    }

    #[inline]
    fn uncompressed_size(&self) -> usize {
        self.as_ref().uncompressed_size()
    }
}

impl<T: CanonicalDeserialize> CanonicalDeserialize for Rc<T> {
    #[inline]
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        Ok(Rc::new(T::deserialize(&mut reader)?))
    }

    #[inline]
    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        Ok(Rc::new(T::deserialize_unchecked(&mut reader)?))
    }

    #[inline]
    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        Ok(Rc::new(T::deserialize_uncompressed(&mut reader)?))
    }

    #[inline]
    fn deserialize_uncompressed_unchecked<R: Read>(
        mut reader: R,
    ) -> Result<Self, SerializationError> {
        Ok(Rc::new(T::deserialize_uncompressed_unchecked(&mut reader)?))
    }
}

// Serialize boolean with a full byte
impl CanonicalSerialize for bool {
    #[inline]
    fn serialize<W: Write>(&self, writer: W) -> Result<(), SerializationError> {
        Ok(CanonicalSerialize::serialize(&(*self as u8), writer)?)
    }

    #[inline]
    fn serialized_size(&self) -> usize {
        1
    }
}

impl CanonicalDeserialize for bool {
    #[inline]
    fn deserialize<R: Read>(reader: R) -> Result<Self, SerializationError> {
        let val = u8::deserialize(reader)?;
        if val == 0 {
            return Ok(false);
        } else if val == 1 {
            return Ok(true);
        }

        Err(SerializationError::InvalidData)
    }

    #[inline]
    fn deserialize_unchecked<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Ok(u8::deserialize(reader)? == 1)
    }

    #[inline]
    fn deserialize_uncompressed_unchecked<R: Read>(reader: R) -> Result<Self, SerializationError> {
        Self::deserialize_unchecked(reader)
    }
}

// Serialize BTreeMap as `len(map) || key 1 || value 1 || ... || key n || value n`
impl<K, V> CanonicalSerialize for BTreeMap<K, V>
where
    K: CanonicalSerialize,
    V: CanonicalSerialize,
{
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        CanonicalSerialize::serialize(&len, &mut writer)?;
        for (k, v) in self.iter() {
            CanonicalSerialize::serialize(k, &mut writer)?;
            CanonicalSerialize::serialize(v, &mut writer)?;
        }
        Ok(())
    }

    fn serialized_size(&self) -> usize {
        8 + self
            .iter()
            .map(|(k, v)| k.serialized_size() + v.serialized_size())
            .sum::<usize>()
    }

    #[inline]
    fn serialize_without_metadata<W: Write>(
        &self,
        mut writer: W,
    ) -> Result<(), SerializationError> {
        for (k, v) in self.iter() {
            k.serialize_without_metadata(&mut writer)?;
            v.serialize_without_metadata(&mut writer)?;
        }
        Ok(())
    }

    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize_uncompressed(&mut writer)?;
        for (k, v) in self.iter() {
            k.serialize_uncompressed(&mut writer)?;
            v.serialize_uncompressed(&mut writer)?;
        }
        Ok(())
    }

    fn uncompressed_size(&self) -> usize {
        8 + self
            .iter()
            .map(|(k, v)| k.uncompressed_size() + v.uncompressed_size())
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

    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = u64::deserialize_unchecked(&mut reader)?;
        let mut map = BTreeMap::new();
        for _ in 0..len {
            map.insert(
                K::deserialize_unchecked(&mut reader)?,
                V::deserialize_unchecked(&mut reader)?,
            );
        }
        Ok(map)
    }

    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = u64::deserialize_uncompressed(&mut reader)?;
        let mut map = BTreeMap::new();
        for _ in 0..len {
            map.insert(
                K::deserialize_uncompressed(&mut reader)?,
                V::deserialize_uncompressed(&mut reader)?,
            );
        }
        Ok(map)
    }

    fn deserialize_uncompressed_unchecked<R: Read>(
        mut reader: R,
    ) -> Result<Self, SerializationError> {
        let len = u64::deserialize_uncompressed_unchecked(&mut reader)?;
        let mut map = BTreeMap::new();
        for _ in 0..len {
            map.insert(
                K::deserialize_uncompressed_unchecked(&mut reader)?,
                V::deserialize_uncompressed_unchecked(&mut reader)?,
            );
        }
        Ok(map)
    }
}

// Serialize BTreeSet as `len(set) || value_1 || ... || value_n`.
impl<T: CanonicalSerialize> CanonicalSerialize for BTreeSet<T> {
    fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        CanonicalSerialize::serialize(&len, &mut writer)?;
        for elem in self.iter() {
            CanonicalSerialize::serialize(elem, &mut writer)?;
        }
        Ok(())
    }

    fn serialized_size(&self) -> usize {
        8 + self
            .iter()
            .map(|elem| elem.serialized_size())
            .sum::<usize>()
    }

    #[inline]
    fn serialize_without_metadata<W: Write>(
        &self,
        mut writer: W,
    ) -> Result<(), SerializationError> {
        for elem in self.iter() {
            elem.serialize_without_metadata(&mut writer)?;
        }
        Ok(())
    }

    fn serialize_uncompressed<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
        let len = self.len() as u64;
        len.serialize_uncompressed(&mut writer)?;
        for elem in self.iter() {
            elem.serialize_uncompressed(&mut writer)?;
        }
        Ok(())
    }

    fn uncompressed_size(&self) -> usize {
        8 + self
            .iter()
            .map(|elem| elem.uncompressed_size())
            .sum::<usize>()
    }
}

impl<T: CanonicalDeserialize + Ord> CanonicalDeserialize for BTreeSet<T> {
    fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = u64::deserialize(&mut reader)?;
        let mut set = BTreeSet::new();
        for _ in 0..len {
            set.insert(T::deserialize(&mut reader)?);
        }
        Ok(set)
    }

    fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = u64::deserialize_unchecked(&mut reader)?;
        let mut set = BTreeSet::new();
        for _ in 0..len {
            set.insert(T::deserialize_unchecked(&mut reader)?);
        }
        Ok(set)
    }

    fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
        let len = u64::deserialize_uncompressed(&mut reader)?;
        let mut set = BTreeSet::new();
        for _ in 0..len {
            set.insert(T::deserialize_uncompressed(&mut reader)?);
        }
        Ok(set)
    }

    fn deserialize_uncompressed_unchecked<R: Read>(
        mut reader: R,
    ) -> Result<Self, SerializationError> {
        let len = u64::deserialize_uncompressed_unchecked(&mut reader)?;
        let mut set = BTreeSet::new();
        for _ in 0..len {
            set.insert(T::deserialize_uncompressed_unchecked(&mut reader)?);
        }
        Ok(set)
    }
}

/// Positive test: Performs a serialize/deserialize (using all available variants), and checks:
///                1) Serialized size is equal to the one returned by the `serialized_size()` function;
///                2) Deserialized elem is equal to the original one;
/// Negative test: Performs a serialize/deserialize (using all available variants), and checks:
///                1) Serialization/deserialization to/from a buffer of insufficient length will result in an error;
///                2) Modifying the serialized data and then deserialize will lead to a deserialization error or
///                   to a elem different from the original one.
#[allow(dead_code)]
pub fn test_canonical_serialize_deserialize<
    T: PartialEq + std::fmt::Debug + CanonicalSerialize + CanonicalDeserialize,
>(
    negative_test: bool,
    data: &T,
) {
    // serialize/deserialize
    {
        let buf_size = data.serialized_size();

        let mut serialized = Vec::with_capacity(buf_size);
        CanonicalSerialize::serialize(data, &mut serialized).unwrap();
        assert_eq!(serialized.len(), buf_size);
        let de = T::deserialize(&serialized[..]).unwrap();
        assert_eq!(data, &de);

        if negative_test {
            let wrong_buf_size = buf_size - 1;
            T::deserialize(&serialized[..wrong_buf_size]).unwrap_err();
            CanonicalSerialize::serialize(data, &mut serialized[..wrong_buf_size]).unwrap_err();

            let wrong_ser_data = serialized.into_iter().map(|b| !b).collect::<Vec<_>>();
            let deser_result = T::deserialize(&wrong_ser_data[..]);
            assert!(deser_result.is_err() || &deser_result.unwrap() != data);
        }
    }

    // serialize/deserialize_unchecked
    {
        let buf_size = data.serialized_size();
        let mut serialized = Vec::with_capacity(buf_size);
        CanonicalSerialize::serialize(data, &mut serialized).unwrap();
        assert_eq!(serialized.len(), buf_size);
        let de = T::deserialize_unchecked(&serialized[..]).unwrap();
        assert_eq!(data, &de);

        if negative_test {
            let wrong_buf_size = buf_size - 1;
            T::deserialize_unchecked(&serialized[..wrong_buf_size]).unwrap_err();
            CanonicalSerialize::serialize(data, &mut serialized[..wrong_buf_size]).unwrap_err();
            let wrong_ser_data = serialized.into_iter().map(|b| !b).collect::<Vec<_>>();
            let deser_result = T::deserialize_unchecked(&wrong_ser_data[..]);
            assert!(deser_result.is_err() || &deser_result.unwrap() != data);
        }
    }

    // serialize_uncompressed/deserialize_uncompressed
    {
        let buf_size = data.uncompressed_size();
        let mut serialized = Vec::with_capacity(buf_size);
        CanonicalSerialize::serialize_uncompressed(data, &mut serialized).unwrap();
        assert_eq!(serialized.len(), buf_size);
        let de = T::deserialize_uncompressed(&serialized[..]).unwrap();
        assert_eq!(data, &de);

        if negative_test {
            let wrong_buf_size = buf_size - 1;
            T::deserialize_uncompressed(&serialized[..wrong_buf_size]).unwrap_err();
            CanonicalSerialize::serialize_uncompressed(data, &mut serialized[..wrong_buf_size])
                .unwrap_err();

            let wrong_ser_data = serialized.into_iter().map(|b| !b).collect::<Vec<_>>();
            let deser_result = T::deserialize_uncompressed(&wrong_ser_data[..]);
            assert!(deser_result.is_err() || &deser_result.unwrap() != data);
        }
    }

    // serialize_uncompressed/deserialize_uncompressed_unchecked
    {
        let buf_size = data.uncompressed_size();
        let mut serialized = Vec::with_capacity(buf_size);
        CanonicalSerialize::serialize_uncompressed(data, &mut serialized).unwrap();
        assert_eq!(serialized.len(), buf_size);
        let de = T::deserialize_uncompressed_unchecked(&serialized[..]).unwrap();
        assert_eq!(data, &de);

        if negative_test {
            let wrong_buf_size = buf_size - 1;
            T::deserialize_uncompressed_unchecked(&serialized[..wrong_buf_size]).unwrap_err();
            CanonicalSerialize::serialize_uncompressed(data, &mut serialized[..wrong_buf_size])
                .unwrap_err();

            let wrong_ser_data = serialized.into_iter().map(|b| !b).collect::<Vec<_>>();
            let deser_result = T::deserialize_uncompressed_unchecked(&wrong_ser_data[..]);
            assert!(deser_result.is_err() || &deser_result.unwrap() != data);
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::RngCore;
    use std::vec;

    // Serialize T, randomly mutate the data, and deserialize it.
    // Ensure it fails.
    // Up to the caller to provide a valid mutation criterion
    // to ensure that this test always fails.
    // This method requires a concrete instance of the data to be provided,
    // to get the serialized size.
    pub fn ensure_non_malleable_encoding<
        T: PartialEq + std::fmt::Debug + CanonicalSerialize + CanonicalDeserialize,
    >(
        data: T,
        valid_mutation: fn(&[u8]) -> bool,
    ) {
        let r = &mut rand::thread_rng();
        let mut serialized = vec![0; data.serialized_size()];
        r.fill_bytes(&mut serialized);
        while !valid_mutation(&serialized) {
            r.fill_bytes(&mut serialized);
        }
        let de = T::deserialize(&serialized[..]);
        assert!(de.is_err());

        let mut serialized = vec![0; data.uncompressed_size()];
        r.fill_bytes(&mut serialized);
        while !valid_mutation(&serialized) {
            r.fill_bytes(&mut serialized);
        }
        let de = T::deserialize_uncompressed(&serialized[..]);
        assert!(de.is_err());
    }

    #[derive(Copy, Clone, Ord, PartialOrd, Eq, PartialEq, Debug)]
    struct Dummy;

    impl CanonicalSerialize for Dummy {
        #[inline]
        fn serialize<W: Write>(&self, mut writer: W) -> Result<(), SerializationError> {
            CanonicalSerialize::serialize(&100u8, &mut writer)
        }

        #[inline]
        fn serialized_size(&self) -> usize {
            100u8.serialized_size()
        }

        #[inline]
        fn serialize_uncompressed<W: Write>(
            &self,
            mut writer: W,
        ) -> Result<(), SerializationError> {
            (&[100u8, 200u8]).serialize_uncompressed(&mut writer)
        }

        #[inline]
        fn uncompressed_size(&self) -> usize {
            (&[100u8, 200u8]).uncompressed_size()
        }
    }

    impl CanonicalDeserialize for Dummy {
        #[inline]
        fn deserialize<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
            let result = u8::deserialize(&mut reader)?;
            assert_eq!(result, 100u8);
            Ok(Dummy)
        }

        #[inline]
        fn deserialize_unchecked<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
            let result = u8::deserialize_unchecked(&mut reader)?;
            assert_eq!(result, 100u8);
            Ok(Dummy)
        }

        #[inline]
        fn deserialize_uncompressed<R: Read>(mut reader: R) -> Result<Self, SerializationError> {
            let result = Vec::<u8>::deserialize_uncompressed(&mut reader)?;
            assert_eq!(result.as_slice(), &[100u8, 200u8]);

            Ok(Dummy)
        }

        #[inline]
        fn deserialize_uncompressed_unchecked<R: Read>(
            mut reader: R,
        ) -> Result<Self, SerializationError> {
            let result = Vec::<u8>::deserialize_uncompressed_unchecked(&mut reader)?;
            assert_eq!(result.as_slice(), &[100u8, 200u8]);

            Ok(Dummy)
        }
    }

    #[test]
    fn test_vec() {
        test_canonical_serialize_deserialize(true, &vec![1u64, 2, 3, 4, 5]);
        test_canonical_serialize_deserialize(true, &Vec::<u64>::new());
    }

    #[test]
    fn test_uint() {
        test_canonical_serialize_deserialize(true, &192830918usize);
        test_canonical_serialize_deserialize(true, &192830918u64);
        test_canonical_serialize_deserialize(true, &192830918u32);
        test_canonical_serialize_deserialize(true, &22313u16);
        test_canonical_serialize_deserialize(true, &123u8);
    }

    #[test]
    fn test_string() {
        test_canonical_serialize_deserialize(true, &String::from("arkworks"));
    }

    #[test]
    fn test_tuple() {
        test_canonical_serialize_deserialize(false, &());
        test_canonical_serialize_deserialize(false, &(123u64, Dummy));
        test_canonical_serialize_deserialize(false, &(123u64, 234u32, Dummy));
    }

    #[test]
    fn test_tuple_vec() {
        test_canonical_serialize_deserialize(
            false,
            &vec![
                (Dummy, Dummy, Dummy),
                (Dummy, Dummy, Dummy),
                (Dummy, Dummy, Dummy),
            ],
        );
        test_canonical_serialize_deserialize(
            false,
            &vec![
                (86u8, 98u64, Dummy),
                (86u8, 98u64, Dummy),
                (86u8, 98u64, Dummy),
            ],
        );
    }

    #[test]
    fn test_option() {
        test_canonical_serialize_deserialize(false, &Some(Dummy));
        test_canonical_serialize_deserialize(false, &None::<Dummy>);

        test_canonical_serialize_deserialize(false, &Some(10u64));
        test_canonical_serialize_deserialize(false, &None::<u64>);
    }

    #[test]
    fn test_rc() {
        test_canonical_serialize_deserialize(false, &Rc::new(Dummy));
    }

    #[test]
    fn test_bool() {
        test_canonical_serialize_deserialize(true, &true);
        test_canonical_serialize_deserialize(false, &false);

        let valid_mutation = |data: &[u8]| -> bool {
            return data.len() == 1 && data[0] > 1;
        };
        for _ in 0..10 {
            ensure_non_malleable_encoding(true, valid_mutation);
            ensure_non_malleable_encoding(false, valid_mutation);
        }
    }

    #[test]
    fn test_btreemap() {
        let mut map = BTreeMap::new();
        map.insert(0u64, Dummy);
        map.insert(5u64, Dummy);
        test_canonical_serialize_deserialize(false, &map);
        let mut map = BTreeMap::new();
        map.insert(10u64, vec![1u8, 2u8, 3u8]);
        map.insert(50u64, vec![4u8, 5u8, 6u8]);
        test_canonical_serialize_deserialize(true, &map);
    }

    #[test]
    fn test_btreeset() {
        let mut set = BTreeSet::new();
        set.insert(Dummy);
        set.insert(Dummy);
        test_canonical_serialize_deserialize(false, &set);
        let mut set = BTreeSet::new();
        set.insert(vec![1u8, 2u8, 3u8]);
        set.insert(vec![4u8, 5u8, 6u8]);
        test_canonical_serialize_deserialize(true, &set);
    }

    #[test]
    fn test_phantomdata() {
        test_canonical_serialize_deserialize(false, &std::marker::PhantomData::<Dummy>);
    }
}
