use crate::{
    bits::{boolean::Boolean, uint8::UInt8},
    Vec,
};
use algebra::Field;
use r1cs_core::SynthesisError;

pub mod boolean;
pub mod uint8;
#[macro_use]
pub mod uint;

make_uint!(UInt16, 16, u16, uint16);
make_uint!(UInt32, 32, u32, uint32);
make_uint!(UInt64, 64, u64, uint64);

pub trait ToBitsGadget<F: Field> {
    /// Outputs the canonical bit-wise representation of `self`.
    ///
    /// This is the correct default for 99% of use cases.
    fn to_bits(&self) -> Result<Vec<Boolean<F>>, SynthesisError>;

    /// Outputs a possibly non-unique bit-wise representation of `self`.
    ///
    /// If you're not absolutely certain that your usecase can get away with a
    /// non-canonical representation, please use `self.to_bits(cs)` instead.
    fn to_non_unique_bits(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        self.to_bits()
    }
}

impl<F: Field> ToBitsGadget<F> for Boolean<F> {
    fn to_bits(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        Ok(vec![self.clone()])
    }
}

impl<F: Field> ToBitsGadget<F> for [Boolean<F>] {
    fn to_bits(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        Ok(self.to_vec())
    }
}

impl<F: Field> ToBitsGadget<F> for Vec<Boolean<F>> {
    fn to_bits(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        Ok(self.clone())
    }
}

impl<F: Field> ToBitsGadget<F> for UInt8<F> {
    fn to_bits(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        Ok(self.into_bits_le())
    }
}

impl<F: Field> ToBitsGadget<F> for [UInt8<F>] {
    fn to_bits(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        let mut result = Vec::with_capacity(&self.len() * 8);
        for byte in self {
            result.extend_from_slice(&byte.into_bits_le());
        }
        Ok(result)
    }
}

impl<F: Field> ToBitsGadget<F> for Vec<UInt8<F>> {
    fn to_bits(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        let mut result = Vec::with_capacity(&self.len() * 8);
        for byte in self {
            result.extend_from_slice(&byte.into_bits_le());
        }
        Ok(result)
    }
}

pub trait ToBytesGadget<F: Field> {
    /// Outputs a canonical byte-wise representation of `self`.
    ///
    /// This is the correct default for 99% of use cases.
    fn to_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError>;

    /// Outputs a possibly non-unique byte decomposition of `self`.
    ///
    /// If you're not absolutely certain that your usecase can get away with a
    /// non-canonical representation, please use `self.to_bytes(cs)` instead.
    fn to_non_unique_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        self.to_bytes()
    }
}

impl<F: Field> ToBytesGadget<F> for [UInt8<F>] {
    fn to_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        Ok(self.to_vec())
    }
}

impl<'a, F: Field, T: 'a + ToBytesGadget<F>> ToBytesGadget<F> for &'a T {
    fn to_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        (*self).to_bytes()
    }
}

impl<'a, F: Field> ToBytesGadget<F> for &'a [UInt8<F>] {
    fn to_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        Ok(self.to_vec())
    }
}
