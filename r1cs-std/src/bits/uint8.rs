use algebra::Field;
use algebra::{FpParameters, PrimeField, ToConstraintField};

use r1cs_core::{ConstraintSystemRef, Namespace, SynthesisError};

use crate::{fields::fp::AllocatedFp, prelude::*, Assignment, Vec};
use core::borrow::Borrow;

/// Represents an interpretation of 8 `Boolean` objects as an
/// unsigned integer.
#[derive(Clone, Debug)]
pub struct UInt8<F: Field> {
    /// Little-endian representation: least significant bit first
    pub(crate) bits: Vec<Boolean<F>>,
    /// Little-endian representation: least significant bit first
    pub(crate) value: Option<u8>,
}

impl<F: Field> R1CSVar<F> for UInt8<F> {
    type Value = u8;

    fn cs(&self) -> Option<ConstraintSystemRef<F>> {
        self.bits.as_slice().cs()
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        let mut value = None;
        for (i, bit) in self.bits.iter().enumerate() {
            let b = u8::from(bit.value()?);
            value = match value {
                Some(value) => Some(value + (b << i)),
                None => Some(b << i),
            };
        }
        debug_assert_eq!(self.value, value);
        value.get()
    }
}

impl<F: Field> UInt8<F> {
    /// Construct a constant vector of `UInt8` from a vector of `u8`
    pub fn constant_vec(values: &[u8]) -> Vec<Self> {
        let mut result = Vec::new();
        for value in values {
            result.push(UInt8::constant(*value));
        }
        result
    }

    /// Construct a constant `UInt8` from a `u8`
    pub fn constant(value: u8) -> Self {
        let mut bits = Vec::with_capacity(8);

        let mut tmp = value;
        for _ in 0..8 {
            // If last bit is one, push one.
            bits.push(Boolean::constant(tmp & 1 == 1));
            tmp >>= 1;
        }

        Self {
            bits,
            value: Some(value),
        }
    }

    pub fn new_witness_vec(
        cs: impl Into<Namespace<F>>,
        values: &[impl Into<Option<u8>> + Copy],
    ) -> Result<Vec<Self>, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        let mut output_vec = Vec::with_capacity(values.len());
        for value in values {
            let byte: Option<u8> = Into::into(*value);
            output_vec.push(Self::new_witness(cs.clone(), || byte.get())?);
        }
        Ok(output_vec)
    }

    /// Allocates a vector of `u8`'s by first converting (chunks of) them to
    /// `ConstraintF` elements, (thus reducing the number of input allocations),
    /// and then converts this list of `ConstraintF` gadgets back into
    /// bytes.
    pub fn new_input_vec(
        cs: impl Into<Namespace<F>>,
        values: &[u8],
    ) -> Result<Vec<Self>, SynthesisError>
    where
        F: PrimeField,
    {
        let ns = cs.into();
        let cs = ns.cs();
        let values_len = values.len();
        let field_elements: Vec<F> = ToConstraintField::<F>::to_field_elements(values).unwrap();

        let max_size = 8 * (F::Params::CAPACITY / 8) as usize;
        let mut allocated_bits = Vec::new();
        for field_element in field_elements.into_iter() {
            let fe = AllocatedFp::new_input(cs.clone(), || Ok(field_element))?;
            let mut fe_bits = fe.to_bits()?;
            // FpGadget::to_bits outputs a big-endian binary representation of
            // fe_gadget's value, so we have to reverse it to get the little-endian
            // form.
            fe_bits.reverse();

            // Remove the most significant bit, because we know it should be zero
            // because `values.to_field_elements()` only
            // packs field elements up to the penultimate bit.
            // That is, the most significant bit (`ConstraintF::NUM_BITS`-th bit) is
            // unset, so we can just pop it off.
            allocated_bits.extend_from_slice(&fe_bits[0..max_size]);
        }

        // Chunk up slices of 8 bit into bytes.
        Ok(allocated_bits[0..8 * values_len]
            .chunks(8)
            .map(Self::from_bits_le)
            .collect())
    }

    /// Turns this `UInt8` into its little-endian byte order representation.
    /// LSB-first means that we can easily get the corresponding field element
    /// via double and add.
    pub fn into_bits_le(&self) -> Vec<Boolean<F>> {
        self.bits.to_vec()
    }

    /// Converts a little-endian byte order representation of bits into a
    /// `UInt8`.
    pub fn from_bits_le(bits: &[Boolean<F>]) -> Self {
        assert_eq!(bits.len(), 8);

        let bits = bits.to_vec();

        let mut value = Some(0u8);
        for b in bits.iter().rev() {
            value.as_mut().map(|v| *v <<= 1);

            match *b {
                Boolean::Constant(b) => {
                    value.as_mut().map(|v| *v |= u8::from(b));
                }
                Boolean::Is(ref b) => match b.value() {
                    Ok(b) => {
                        value.as_mut().map(|v| *v |= u8::from(b));
                    }
                    Err(_) => value = None,
                },
                Boolean::Not(ref b) => match b.value() {
                    Ok(b) => {
                        value.as_mut().map(|v| *v |= u8::from(!b));
                    }
                    Err(_) => value = None,
                },
            }
        }

        Self { value, bits }
    }

    /// XOR this `UInt8` with another `UInt8`
    pub fn xor(&self, other: &Self) -> Result<Self, SynthesisError> {
        let new_value = match (self.value, other.value) {
            (Some(a), Some(b)) => Some(a ^ b),
            _ => None,
        };

        let bits = self
            .bits
            .iter()
            .zip(other.bits.iter())
            .map(|(a, b)| a.xor(b))
            .collect::<Result<_, _>>()?;

        Ok(Self {
            bits,
            value: new_value,
        })
    }
}

impl<ConstraintF: Field> EqGadget<ConstraintF> for UInt8<ConstraintF> {
    fn is_eq(&self, other: &Self) -> Result<Boolean<ConstraintF>, SynthesisError> {
        self.bits.as_slice().is_eq(&other.bits)
    }

    fn conditional_enforce_equal(
        &self,
        other: &Self,
        condition: &Boolean<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        self.bits.conditional_enforce_equal(&other.bits, condition)
    }

    fn conditional_enforce_not_equal(
        &self,
        other: &Self,
        condition: &Boolean<ConstraintF>,
    ) -> Result<(), SynthesisError> {
        self.bits
            .conditional_enforce_not_equal(&other.bits, condition)
    }
}

impl<ConstraintF: Field> AllocVar<u8, ConstraintF> for UInt8<ConstraintF> {
    fn new_variable<T: Borrow<u8>>(
        cs: impl Into<Namespace<ConstraintF>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        let value = f().map(|f| *f.borrow());
        let values = match value {
            Ok(val) => (0..8).map(|i| Some((val >> i) & 1 == 1)).collect(),
            _ => vec![None; 8],
        };
        let bits = values
            .into_iter()
            .map(|v| Boolean::new_variable(cs.clone(), || v.get(), mode))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            bits,
            value: value.ok(),
        })
    }
}

#[cfg(test)]
mod test {
    use super::UInt8;
    use crate::{prelude::*, Vec};
    use algebra::bls12_381::Fr;
    use r1cs_core::{ConstraintSystem, SynthesisError};
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_uint8_from_bits_to_bits() -> Result<(), SynthesisError> {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let byte_val = 0b01110001;
        let byte = UInt8::new_witness(cs.ns("alloc value"), || Ok(byte_val)).unwrap();
        let bits = byte.into_bits_le();
        for (i, bit) in bits.iter().enumerate() {
            assert_eq!(bit.value()?, (byte_val >> i) & 1 == 1)
        }
        Ok(())
    }

    #[test]
    fn test_uint8_new_input_vec() -> Result<(), SynthesisError> {
        let cs = ConstraintSystem::<Fr>::new_ref();
        let byte_vals = (64u8..128u8).collect::<Vec<_>>();
        let bytes = UInt8::new_input_vec(cs.ns("alloc value"), &byte_vals).unwrap();
        for (native, variable) in byte_vals.into_iter().zip(bytes) {
            let bits = variable.into_bits_le();
            for (i, bit) in bits.iter().enumerate() {
                assert_eq!(bit.value()?, (native >> i) & 1 == 1)
            }
        }
        Ok(())
    }

    #[test]
    fn test_uint8_from_bits() -> Result<(), SynthesisError> {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..1000 {
            let v = (0..8)
                .map(|_| Boolean::<Fr>::Constant(rng.gen()))
                .collect::<Vec<_>>();

            let val = UInt8::from_bits_le(&v);

            for (i, bit) in val.bits.iter().enumerate() {
                match bit {
                    Boolean::Constant(b) => assert!(*b == ((val.value()? >> i) & 1 == 1)),
                    _ => unreachable!(),
                }
            }

            let expected_to_be_same = val.into_bits_le();

            for x in v.iter().zip(expected_to_be_same.iter()) {
                match x {
                    (&Boolean::Constant(true), &Boolean::Constant(true)) => {}
                    (&Boolean::Constant(false), &Boolean::Constant(false)) => {}
                    _ => unreachable!(),
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_uint8_xor() -> Result<(), SynthesisError> {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..1000 {
            let cs = ConstraintSystem::<Fr>::new_ref();

            let a: u8 = rng.gen();
            let b: u8 = rng.gen();
            let c: u8 = rng.gen();

            let mut expected = a ^ b ^ c;

            let a_bit = UInt8::new_witness(cs.ns("a_bit"), || Ok(a)).unwrap();
            let b_bit = UInt8::constant(b);
            let c_bit = UInt8::new_witness(cs.ns("c_bit"), || Ok(c)).unwrap();

            let r = a_bit.xor(&b_bit).unwrap();
            let r = r.xor(&c_bit).unwrap();

            assert!(cs.is_satisfied().unwrap());

            assert!(r.value == Some(expected));

            for b in r.bits.iter() {
                match b {
                    Boolean::Is(b) => assert!(b.value()? == (expected & 1 == 1)),
                    Boolean::Not(b) => assert!(!b.value()? == (expected & 1 == 1)),
                    Boolean::Constant(b) => assert!(*b == (expected & 1 == 1)),
                }

                expected >>= 1;
            }
        }
        Ok(())
    }
}
