use algebra::{Field, FpParameters, PrimeField};

use r1cs_core::{ConstraintSystem, LinearCombination, SynthesisError};

use crate::{
    boolean::{AllocatedBit, Boolean},
    prelude::*,
    Assignment,
};

/// Represents an interpretation of 64 `Boolean` objects as an
/// unsigned integer.
#[derive(Clone, Debug)]
pub struct UInt64 {
    // Least significant bit_gadget first
    bits: Vec<Boolean>,
    value: Option<u64>,
}

impl UInt64 {
    pub fn get_value(&self) -> Option<u64> {
        self.value
    }

    /// Construct a constant `UInt64` from a `u64`
    pub fn constant(value: u64) -> Self {
        let mut bits = Vec::with_capacity(64);

        let mut tmp = value;
        for _ in 0..64 {
            if tmp & 1 == 1 {
                bits.push(Boolean::constant(true))
            } else {
                bits.push(Boolean::constant(false))
            }

            tmp >>= 1;
        }

        UInt64 {
            bits,
            value: Some(value),
        }
    }

    /// Allocate a `UInt64` in the constraint system
    pub fn alloc<ConstraintF, CS>(mut cs: CS, value: Option<u64>) -> Result<Self, SynthesisError>
    where
        ConstraintF: Field,
        CS: ConstraintSystem<ConstraintF>,
    {
        let values = match value {
            Some(mut val) => {
                let mut v = Vec::with_capacity(64);

                for _ in 0..64 {
                    v.push(Some(val & 1 == 1));
                    val >>= 1;
                }

                v
            }
            None => vec![None; 64],
        };

        let bits = values
            .into_iter()
            .enumerate()
            .map(|(i, v)| {
                Ok(Boolean::from(AllocatedBit::alloc(
                    cs.ns(|| format!("allocated bit_gadget {}", i)),
                    || v.get(),
                )?))
            })
            .collect::<Result<Vec<_>, SynthesisError>>()?;

        Ok(UInt64 { bits, value })
    }

    /// Turns this `UInt64` into its little-endian byte order representation.
    pub fn to_bits_le(&self) -> Vec<Boolean> {
        self.bits.clone()
    }

    /// Converts a little-endian byte order representation of bits into a
    /// `UInt64`.
    pub fn from_bits_le(bits: &[Boolean]) -> Self {
        assert_eq!(bits.len(), 64);

        let bits = bits.to_vec();

        let mut value = Some(0u64);
        for b in bits.iter().rev() {
            value.as_mut().map(|v| *v <<= 1);

            match b {
                &Boolean::Constant(b) => {
                    if b {
                        value.as_mut().map(|v| *v |= 1);
                    }
                }
                &Boolean::Is(ref b) => match b.get_value() {
                    Some(true) => {
                        value.as_mut().map(|v| *v |= 1);
                    }
                    Some(false) => {}
                    None => value = None,
                },
                &Boolean::Not(ref b) => match b.get_value() {
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

    pub fn rotr(&self, by: usize) -> Self {
        let by = by % 64;

        let new_bits = self
            .bits
            .iter()
            .skip(by)
            .chain(self.bits.iter())
            .take(64)
            .cloned()
            .collect();

        UInt64 {
            bits: new_bits,
            value: self.value.map(|v| v.rotate_right(by as u32)),
        }
    }

    /// XOR this `UInt64` with another `UInt64`
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

        Ok(UInt64 {
            bits,
            value: new_value,
        })
    }

    /// Perform modular addition of several `UInt64` objects.
    pub fn addmany<ConstraintF, CS>(mut cs: CS, operands: &[Self]) -> Result<Self, SynthesisError>
    where
        ConstraintF: PrimeField,
        CS: ConstraintSystem<ConstraintF>,
    {
        // Make some arbitrary bounds for ourselves to avoid overflows
        // in the scalar field
        assert!(ConstraintF::Params::MODULUS_BITS >= 128);

        assert!(operands.len() >= 1);
        assert!(operands.len() <= 10);

        if operands.len() == 1 {
            return Ok(operands[0].clone());
        }

        // Compute the maximum value of the sum so we allocate enough bits for
        // the result
        let mut max_value = (operands.len() as u128) * u128::from(u64::max_value());

        // Keep track of the resulting value
        let mut result_value = Some(0u64 as u128);

        // This is a linear combination that we will enforce to be "zero"
        let mut lc = LinearCombination::zero();

        let mut all_constants = true;

        // Iterate over the operands
        for op in operands {
            // Accumulate the value
            match op.value {
                Some(val) => {
                    result_value.as_mut().map(|v| *v += u128::from(val));
                }
                None => {
                    // If any of our operands have unknown value, we won't
                    // know the value of the result
                    result_value = None;
                }
            }

            // Iterate over each bit_gadget of the operand and add the operand to
            // the linear combination
            let mut coeff = ConstraintF::one();
            for bit in &op.bits {
                match *bit {
                    Boolean::Is(ref bit) => {
                        all_constants = false;

                        // Add coeff * bit_gadget
                        lc += (coeff, bit.get_variable());
                    }
                    Boolean::Not(ref bit) => {
                        all_constants = false;

                        // Add coeff * (1 - bit_gadget) = coeff * ONE - coeff * bit_gadget
                        lc = lc + (coeff, CS::one()) - (coeff, bit.get_variable());
                    }
                    Boolean::Constant(bit) => {
                        if bit {
                            lc += (coeff, CS::one());
                        }
                    }
                }

                coeff.double_in_place();
            }
        }

        // The value of the actual result is modulo 2^64
        let modular_value = result_value.map(|v| v as u64);

        if all_constants && modular_value.is_some() {
            // We can just return a constant, rather than
            // unpacking the result into allocated bits.

            return Ok(UInt64::constant(modular_value.unwrap()));
        }

        // Storage area for the resulting bits
        let mut result_bits = vec![];

        // Allocate each bit_gadget of the result
        let mut coeff = ConstraintF::one();
        let mut i = 0;
        while max_value != 0 {
            // Allocate the bit_gadget
            let b = AllocatedBit::alloc(cs.ns(|| format!("result bit_gadget {}", i)), || {
                result_value.map(|v| (v >> i) & 1 == 1).get()
            })?;

            // Subtract this bit_gadget from the linear combination to ensure the sums
            // balance out
            lc = lc - (coeff, b.get_variable());

            result_bits.push(b.into());

            max_value >>= 1;
            i += 1;
            coeff.double_in_place();
        }

        // Enforce that the linear combination equals zero
        cs.enforce(|| "modular addition", |lc| lc, |lc| lc, |_| lc);

        // Discard carry bits that we don't care about
        result_bits.truncate(64);

        Ok(UInt64 {
            bits: result_bits,
            value: modular_value,
        })
    }
}

impl<ConstraintF: Field> ToBytesGadget<ConstraintF> for UInt64 {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let value_chunks = match self.value.map(|val| {
            use algebra::bytes::ToBytes;
            let mut bytes = [0u8; 8];
            val.write(bytes.as_mut()).unwrap();
            bytes
        }) {
            Some(chunks) => [
                Some(chunks[0]),
                Some(chunks[1]),
                Some(chunks[2]),
                Some(chunks[3]),
                Some(chunks[4]),
                Some(chunks[5]),
                Some(chunks[6]),
                Some(chunks[7]),
            ],
            None => [None, None, None, None, None, None, None, None],
        };
        let mut bytes = Vec::new();
        for (i, chunk8) in self.to_bits_le().chunks(8).enumerate() {
            let byte = UInt8 {
                bits: chunk8.to_vec(),
                value: value_chunks[i],
            };
            bytes.push(byte);
        }

        Ok(bytes)
    }

    fn to_bytes_strict<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.to_bytes(cs)
    }
}

impl PartialEq for UInt64 {
    fn eq(&self, other: &Self) -> bool {
        self.value.is_some() && other.value.is_some() && self.value == other.value
    }
}

impl Eq for UInt64 {}

impl<ConstraintF: Field> EqGadget<ConstraintF> for UInt64 {
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
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.bits
            .conditional_enforce_equal(cs, &other.bits, should_enforce)
    }

    fn conditional_enforce_not_equal<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
        should_enforce: &Boolean,
    ) -> Result<(), SynthesisError> {
        self.bits
            .conditional_enforce_not_equal(cs, &other.bits, should_enforce)
    }
}

#[cfg(test)]
mod test {
    use super::UInt64;
    use crate::{bits::boolean::Boolean, test_constraint_system::TestConstraintSystem};
    use algebra::fields::{bls12_381::Fr, Field};
    use r1cs_core::ConstraintSystem;
    use rand::{Rng, SeedableRng};
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_uint64_from_bits() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..1000 {
            let v = (0..64)
                .map(|_| Boolean::constant(rng.gen()))
                .collect::<Vec<_>>();

            let b = UInt64::from_bits_le(&v);

            for (i, bit_gadget) in b.bits.iter().enumerate() {
                match bit_gadget {
                    &Boolean::Constant(bit_gadget) => {
                        assert!(bit_gadget == ((b.value.unwrap() >> i) & 1 == 1));
                    }
                    _ => unreachable!(),
                }
            }

            let expected_to_be_same = b.to_bits_le();

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
    fn test_uint64_xor() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..1000 {
            let mut cs = TestConstraintSystem::<Fr>::new();

            let a: u64 = rng.gen();
            let b: u64 = rng.gen();
            let c: u64 = rng.gen();

            let mut expected = a ^ b ^ c;

            let a_bit = UInt64::alloc(cs.ns(|| "a_bit"), Some(a)).unwrap();
            let b_bit = UInt64::constant(b);
            let c_bit = UInt64::alloc(cs.ns(|| "c_bit"), Some(c)).unwrap();

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

    #[test]
    fn test_uint64_addmany_constants() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..1000 {
            let mut cs = TestConstraintSystem::<Fr>::new();

            let a: u64 = rng.gen();
            let b: u64 = rng.gen();
            let c: u64 = rng.gen();

            let a_bit = UInt64::constant(a);
            let b_bit = UInt64::constant(b);
            let c_bit = UInt64::constant(c);

            let mut expected = a.wrapping_add(b).wrapping_add(c);

            let r = UInt64::addmany(cs.ns(|| "addition"), &[a_bit, b_bit, c_bit]).unwrap();

            assert!(r.value == Some(expected));

            for b in r.bits.iter() {
                match b {
                    &Boolean::Is(_) => panic!(),
                    &Boolean::Not(_) => panic!(),
                    &Boolean::Constant(b) => {
                        assert!(b == (expected & 1 == 1));
                    }
                }

                expected >>= 1;
            }
        }
    }

    #[test]
    fn test_uint64_addmany() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..1000 {
            let mut cs = TestConstraintSystem::<Fr>::new();

            let a: u64 = rng.gen();
            let b: u64 = rng.gen();
            let c: u64 = rng.gen();
            let d: u64 = rng.gen();

            let mut expected = (a ^ b).wrapping_add(c).wrapping_add(d);

            let a_bit = UInt64::alloc(cs.ns(|| "a_bit"), Some(a)).unwrap();
            let b_bit = UInt64::constant(b);
            let c_bit = UInt64::constant(c);
            let d_bit = UInt64::alloc(cs.ns(|| "d_bit"), Some(d)).unwrap();

            let r = a_bit.xor(cs.ns(|| "xor"), &b_bit).unwrap();
            let r = UInt64::addmany(cs.ns(|| "addition"), &[r, c_bit, d_bit]).unwrap();

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
                    &Boolean::Constant(_) => unreachable!(),
                }

                expected >>= 1;
            }

            // Flip a bit_gadget and see if the addition constraint still works
            if cs.get("addition/result bit_gadget 0/boolean").is_zero() {
                cs.set("addition/result bit_gadget 0/boolean", Fr::one());
            } else {
                cs.set("addition/result bit_gadget 0/boolean", Fr::zero());
            }

            assert!(!cs.is_satisfied());
        }
    }

    #[test]
    fn test_uint64_rotr() {
        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        let mut num = rng.gen();

        let a = UInt64::constant(num);

        for i in 0..64 {
            let b = a.rotr(i);

            assert!(b.value.unwrap() == num);

            let mut tmp = num;
            for b in &b.bits {
                match b {
                    &Boolean::Constant(b) => {
                        assert_eq!(b, tmp & 1 == 1);
                    }
                    _ => unreachable!(),
                }

                tmp >>= 1;
            }

            num = num.rotate_right(1);
        }
    }
}
