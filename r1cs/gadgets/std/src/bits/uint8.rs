use algebra::{Field, FpParameters, PrimeField, ToConstraintField};

use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{boolean::AllocatedBit, fields::fp::FpGadget, prelude::*, Assignment};
use std::borrow::Borrow;

/// Represents an interpretation of 8 `Boolean` objects as an
/// unsigned integer.
#[derive(Clone, Debug)]
pub struct UInt8 {
    // Least significant bit_gadget first
    pub(crate) bits: Vec<Boolean>,
    pub(crate) value: Option<u8>,
}

impl UInt8 {
    pub fn get_value(&self) -> Option<u8> {
        self.value
    }

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
            if tmp & 1 == 1 {
                bits.push(Boolean::constant(true))
            } else {
                bits.push(Boolean::constant(false))
            }

            tmp >>= 1;
        }

        Self {
            bits,
            value: Some(value),
        }
    }

    pub fn alloc_vec<ConstraintF, CS, T>(
        mut cs: CS,
        values: &[T],
    ) -> Result<Vec<Self>, SynthesisError>
    where
        ConstraintF: Field,
        CS: ConstraintSystem<ConstraintF>,
        T: Into<Option<u8>> + Copy,
    {
        let mut output_vec = Vec::with_capacity(values.len());
        for (i, value) in values.into_iter().enumerate() {
            let byte: Option<u8> = Into::into(*value);
            let alloc_byte = Self::alloc(&mut cs.ns(|| format!("byte_{}", i)), || byte.get())?;
            output_vec.push(alloc_byte);
        }
        Ok(output_vec)
    }

    /// Allocates a vector of `u8`'s by first converting (chunks of) them to
    /// `ConstraintF` elements, (thus reducing the number of input allocations),
    /// and then converts this list of `ConstraintF` gadgets back into
    /// bytes.
    pub fn alloc_input_vec<ConstraintF, CS>(
        mut cs: CS,
        values: &[u8],
    ) -> Result<Vec<Self>, SynthesisError>
    where
        ConstraintF: PrimeField,
        CS: ConstraintSystem<ConstraintF>,
    {
        let values_len = values.len();
        let field_elements: Vec<ConstraintF> =
            ToConstraintField::<ConstraintF>::to_field_elements(values).unwrap();

        let max_size = (<ConstraintF as PrimeField>::Params::CAPACITY / 8) as usize;

        let mut allocated_bits = Vec::new();
        for (i, (field_element, byte_chunk)) in field_elements
            .into_iter()
            .zip(values.chunks(max_size))
            .enumerate()
        {
            let fe = FpGadget::alloc_input(&mut cs.ns(|| format!("Field element {}", i)), || {
                Ok(field_element)
            })?;

            // Let's use the length-restricted variant of the ToBitsGadget to remove the
            // padding: the padding bits are not constrained to be zero, so any field element
            // passed as input (as long as it has the last bits set to the proper value) can
            // satisfy the constraints. This kind of freedom might not be desiderable in
            // recursive SNARK circuits, where the public inputs of the inner circuit are
            // usually involved in other kind of constraints inside the wrap circuit.
            let to_skip: usize =
                <ConstraintF as PrimeField>::Params::MODULUS_BITS as usize - (byte_chunk.len() * 8);
            let mut fe_bits = fe.to_bits_with_length_restriction(
                cs.ns(|| format!("Convert fe to bits {}", i)),
                to_skip,
            )?;

            // FpGadget::to_bits outputs a big-endian binary representation of
            // fe_gadget's value, so we have to reverse it to get the little-endian
            // form.
            fe_bits.reverse();

            allocated_bits.extend_from_slice(fe_bits.as_slice());
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
    pub fn into_bits_le(&self) -> Vec<Boolean> {
        self.bits.iter().cloned().collect()
    }

    /// Converts a little-endian byte order representation of bits into a
    /// `UInt8`.
    pub fn from_bits_le(bits: &[Boolean]) -> Self {
        assert_eq!(bits.len(), 8);

        let bits = bits.to_vec();

        let mut value = Some(0u8);
        for b in bits.iter().rev() {
            value.as_mut().map(|v| *v <<= 1);

            match *b {
                Boolean::Constant(b) => {
                    if b {
                        value.as_mut().map(|v| *v |= 1);
                    }
                }
                Boolean::Is(ref b) => match b.get_value() {
                    Some(true) => {
                        value.as_mut().map(|v| *v |= 1);
                    }
                    Some(false) => {}
                    None => value = None,
                },
                Boolean::Not(ref b) => match b.get_value() {
                    Some(false) => {
                        value.as_mut().map(|v| *v |= 1);
                    }
                    Some(true) => {}
                    None => value = None,
                },
            }
        }

        Self { value, bits }
    }

    /// XOR this `UInt8` with another `UInt8`
    pub fn xor<ConstraintF, CS>(&self, mut cs: CS, other: &Self) -> Result<Self, SynthesisError>
    where
        ConstraintF: Field,
        CS: ConstraintSystem<ConstraintF>,
    {
        let new_value = match (self.value, other.value) {
            (Some(a), Some(b)) => Some(a ^ b),
            _ => None,
        };

        let bits = self
            .bits
            .iter()
            .zip(other.bits.iter())
            .enumerate()
            .map(|(i, (a, b))| Boolean::xor(cs.ns(|| format!("xor of bit_gadget {}", i)), a, b))
            .collect::<Result<_, _>>()?;

        Ok(Self {
            bits,
            value: new_value,
        })
    }

    /// OR this `UInt8` with another `UInt8`
    pub fn or<ConstraintF, CS>(&self, mut cs: CS, other: &Self) -> Result<Self, SynthesisError>
    where
        ConstraintF: Field,
        CS: ConstraintSystem<ConstraintF>,
    {
        let new_value = match (self.value, other.value) {
            (Some(a), Some(b)) => Some(a | b),
            _ => None,
        };

        let bits = self
            .bits
            .iter()
            .zip(other.bits.iter())
            .enumerate()
            .map(|(i, (a, b))| Boolean::or(cs.ns(|| format!("or of bit_gadget {}", i)), a, b))
            .collect::<Result<_, _>>()?;

        Ok(Self {
            bits,
            value: new_value,
        })
    }
}

impl PartialEq for UInt8 {
    fn eq(&self, other: &Self) -> bool {
        !self.value.is_none() && !other.value.is_none() && self.value == other.value
    }
}

impl Eq for UInt8 {}

impl<ConstraintF: Field> EqGadget<ConstraintF> for UInt8 {
    fn is_eq<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<Boolean, SynthesisError> {
        self.bits.as_slice().is_eq(cs, &other.bits)
    }

    fn conditional_enforce_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.bits
            .conditional_enforce_equal(cs, &other.bits, condition)
    }

    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.bits
            .conditional_enforce_not_equal(cs, &other.bits, condition)
    }
}

impl<ConstraintF: Field> AllocGadget<u8, ConstraintF> for UInt8 {
    fn alloc<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<u8>,
    {
        let value = value_gen().map(|val| *val.borrow());
        let values = match value {
            Ok(mut val) => {
                let mut v = Vec::with_capacity(8);

                for _ in 0..8 {
                    v.push(Some(val & 1 == 1));
                    val >>= 1;
                }

                v
            }
            _ => vec![None; 8],
        };

        let bits = values
            .into_iter()
            .enumerate()
            .map(|(i, v)| {
                Ok(Boolean::from(AllocatedBit::alloc(
                    &mut cs.ns(|| format!("allocated bit_gadget {}", i)),
                    || v.ok_or(SynthesisError::AssignmentMissing),
                )?))
            })
            .collect::<Result<Vec<_>, SynthesisError>>()?;

        Ok(Self {
            bits,
            value: value.ok(),
        })
    }

    fn alloc_input<F, T, CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<u8>,
    {
        let value = value_gen().map(|val| *val.borrow());
        let values = match value {
            Ok(mut val) => {
                let mut v = Vec::with_capacity(8);
                for _ in 0..8 {
                    v.push(Some(val & 1 == 1));
                    val >>= 1;
                }

                v
            }
            _ => vec![None; 8],
        };

        let bits = values
            .into_iter()
            .enumerate()
            .map(|(i, v)| {
                Ok(Boolean::from(AllocatedBit::alloc_input(
                    &mut cs.ns(|| format!("allocated bit_gadget {}", i)),
                    || v.ok_or(SynthesisError::AssignmentMissing),
                )?))
            })
            .collect::<Result<Vec<_>, SynthesisError>>()?;

        Ok(Self {
            bits,
            value: value.ok(),
        })
    }
}

#[cfg(test)]
mod test {
    use super::UInt8;
    use crate::{prelude::*, test_constraint_system::TestConstraintSystem};
    use algebra::fields::bls12_381::Fr;
    use r1cs_core::ConstraintSystem;
    use rand::{Rng, RngCore, SeedableRng};
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_uint8_from_bits_to_bits() {
        let mut cs = TestConstraintSystem::<Fr>::new();
        let byte_val = 0b01110001;
        let byte = UInt8::alloc(cs.ns(|| "alloc value"), || Ok(byte_val)).unwrap();
        let bits = byte.into_bits_le();
        for (i, bit) in bits.iter().enumerate() {
            assert_eq!(bit.get_value().unwrap(), (byte_val >> i) & 1 == 1)
        }
    }

    #[test]
    fn test_uint8_alloc_input_vec() {
        use algebra::{to_bytes, Field, FpParameters, PrimeField, ToBytes, UniformRand};
        use rand::thread_rng;

        let mut cs = TestConstraintSystem::<Fr>::new();
        let rng = &mut thread_rng();

        //Random test
        let samples = 100;
        for i in 0..samples {
            // Test with random field
            let byte_vals = to_bytes!(Fr::rand(rng)).unwrap();
            let bytes =
                UInt8::alloc_input_vec(cs.ns(|| format!("alloc value {}", i)), &byte_vals).unwrap();
            assert_eq!(byte_vals.len(), bytes.len());
            for (native_byte, gadget_byte) in byte_vals.into_iter().zip(bytes) {
                assert_eq!(gadget_byte.get_value().unwrap(), native_byte);
            }

            // Test with random bytes
            let mut byte_vals = vec![0u8; rng.gen_range(1..200)];
            rng.fill_bytes(byte_vals.as_mut_slice());
            let bytes = UInt8::alloc_input_vec(cs.ns(|| format!("alloc random {}", i)), &byte_vals)
                .unwrap();
            assert_eq!(byte_vals.len(), bytes.len());
            for (native_byte, gadget_byte) in byte_vals.into_iter().zip(bytes) {
                assert_eq!(gadget_byte.get_value().unwrap(), native_byte);
            }
        }

        //Test one
        let byte_vals = to_bytes!(Fr::one()).unwrap();
        let bytes = UInt8::alloc_input_vec(cs.ns(|| "alloc one bytes"), &byte_vals).unwrap();
        assert_eq!(byte_vals.len(), bytes.len());
        for (native_byte, gadget_byte) in byte_vals.into_iter().zip(bytes) {
            assert_eq!(gadget_byte.get_value().unwrap(), native_byte);
        }

        //Test zero
        let byte_vals = to_bytes!(Fr::zero()).unwrap();
        let bytes = UInt8::alloc_input_vec(cs.ns(|| "alloc zero bytes"), &byte_vals).unwrap();
        assert_eq!(byte_vals.len(), bytes.len());
        for (native_byte, gadget_byte) in byte_vals.into_iter().zip(bytes) {
            assert_eq!(gadget_byte.get_value().unwrap(), native_byte);
        }

        //Test over the modulus byte vec
        let byte_vals = vec![
            std::u8::MAX;
            ((<Fr as PrimeField>::Params::MODULUS_BITS
                + <Fr as PrimeField>::Params::REPR_SHAVE_BITS)
                / 8) as usize
        ];
        let bytes = UInt8::alloc_input_vec(cs.ns(|| "alloc all 1s byte vec"), &byte_vals).unwrap();
        assert_eq!(byte_vals.len(), bytes.len());
        for (native_byte, gadget_byte) in byte_vals.into_iter().zip(bytes) {
            assert_eq!(gadget_byte.get_value().unwrap(), native_byte);
        }
    }

    #[test]
    fn test_uint8_from_bits() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..1000 {
            let v = (0..8)
                .map(|_| Boolean::constant(rng.gen()))
                .collect::<Vec<_>>();

            let b = UInt8::from_bits_le(&v);

            for (i, bit_gadget) in b.bits.iter().enumerate() {
                match bit_gadget {
                    &Boolean::Constant(bit_gadget) => {
                        assert!(bit_gadget == ((b.value.unwrap() >> i) & 1 == 1));
                    }
                    _ => unreachable!(),
                }
            }

            let expected_to_be_same = b.into_bits_le();

            for x in v.iter().zip(expected_to_be_same.iter()) {
                match x {
                    (&Boolean::Constant(true), &Boolean::Constant(true)) => {}
                    (&Boolean::Constant(false), &Boolean::Constant(false)) => {}
                    _ => unreachable!(),
                }
            }
        }
    }

    #[test]
    fn test_uint8_xor() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..1000 {
            let mut cs = TestConstraintSystem::<Fr>::new();

            let a: u8 = rng.gen();
            let b: u8 = rng.gen();
            let c: u8 = rng.gen();

            let mut expected = a ^ b ^ c;

            let a_bit = UInt8::alloc(cs.ns(|| "a_bit"), || Ok(a)).unwrap();
            let b_bit = UInt8::constant(b);
            let c_bit = UInt8::alloc(cs.ns(|| "c_bit"), || Ok(c)).unwrap();

            let r = a_bit.xor(cs.ns(|| "first xor"), &b_bit).unwrap();
            let r = r.xor(cs.ns(|| "second xor"), &c_bit).unwrap();

            assert!(cs.is_satisfied());

            assert!(r.value == Some(expected));

            for b in r.bits.iter() {
                match b {
                    &Boolean::Is(ref b) => {
                        assert!(b.get_value().unwrap() == (expected & 1 == 1));
                    }
                    &Boolean::Not(ref b) => {
                        assert!(!b.get_value().unwrap() == (expected & 1 == 1));
                    }
                    &Boolean::Constant(b) => {
                        assert!(b == (expected & 1 == 1));
                    }
                }

                expected >>= 1;
            }
        }
    }
}
