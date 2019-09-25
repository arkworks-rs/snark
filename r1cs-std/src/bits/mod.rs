use crate::bits::{boolean::Boolean, uint8::UInt8};
use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};

pub mod boolean;
pub mod uint32;
pub mod uint8;

pub trait ToBitsGadget<ConstraintF: Field> {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Vec<Boolean>, SynthesisError>;

    /// Additionally checks if the produced list of booleans is 'valid'.
    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError>;
}

impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for Boolean {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, _: CS) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(vec![self.clone()])
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(vec![self.clone()])
    }
}

impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for [Boolean] {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, _cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.to_vec())
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.to_vec())
    }
}
impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for Vec<Boolean> {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, _cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.clone())
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        Ok(self.clone())
    }
}

impl<ConstraintF: Field> ToBitsGadget<ConstraintF> for [UInt8] {
    fn to_bits<CS: ConstraintSystem<ConstraintF>>(&self, _cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        let mut result = Vec::with_capacity(&self.len() * 8);
        for byte in self {
            result.extend_from_slice(&byte.into_bits_le());
        }
        Ok(result)
    }

    fn to_bits_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        self.to_bits(cs)
    }
}

pub trait ToBytesGadget<ConstraintF: Field> {
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Vec<UInt8>, SynthesisError>;

    /// Additionally checks if the produced list of booleans is 'valid'.
    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError>;
}

impl<ConstraintF: Field> ToBytesGadget<ConstraintF> for [UInt8] {
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(&self, _cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        Ok(self.to_vec())
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}

impl<'a, ConstraintF: Field, T: 'a + ToBytesGadget<ConstraintF>> ToBytesGadget<ConstraintF> for &'a T {
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        (*self).to_bytes(cs)
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}

impl<'a, ConstraintF: Field> ToBytesGadget<ConstraintF> for &'a [UInt8] {
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(&self, _cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        Ok(self.to_vec())
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}
