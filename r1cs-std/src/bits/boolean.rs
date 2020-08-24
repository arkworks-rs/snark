use algebra::{BitIterator, Field};

use crate::{prelude::*, Assignment, Vec};
use core::borrow::Borrow;
use r1cs_core::{lc, ConstraintSystemRef, LinearCombination, Namespace, SynthesisError, Variable};

/// Represents a variable in the constraint system which is guaranteed
/// to be either zero or one.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AllocatedBit<F: Field> {
    variable: Variable,
    cs: ConstraintSystemRef<F>,
}

pub(crate) fn bool_to_field<F: Field>(val: impl Borrow<bool>) -> F {
    if *val.borrow() {
        F::one()
    } else {
        F::zero()
    }
}

impl<F: Field> AllocatedBit<F> {
    /// Get the assigned value for `self`.
    pub fn value(&self) -> Result<bool, SynthesisError> {
        let value = self.cs.assigned_value(self.variable).get()?;
        if value.is_zero() {
            Ok(false)
        } else if value.is_one() {
            Ok(true)
        } else {
            unreachable!("Incorrect value assigned: {:?}", value);
        }
    }

    /// Get the R1CS variable for `self`.
    pub fn variable(&self) -> Variable {
        self.variable
    }

    /// Allocate a witness variable without a booleanity check.
    fn new_witness_without_booleanity_check<T: Borrow<bool>>(
        cs: ConstraintSystemRef<F>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
    ) -> Result<Self, SynthesisError> {
        let variable = cs.new_witness_variable(|| f().map(bool_to_field))?;
        Ok(Self { variable, cs })
    }

    /// Performs an XOR operation over the two operands, returning
    /// an `AllocatedBit`.
    pub fn xor(&self, b: &Self) -> Result<Self, SynthesisError> {
        let result = Self::new_witness_without_booleanity_check(self.cs.clone(), || {
            Ok(self.value()? ^ b.value()?)
        })?;

        // Constrain (a + a) * (b) = (a + b - c)
        // Given that a and b are boolean constrained, if they
        // are equal, the only solution for c is 0, and if they
        // are different, the only solution for c is 1.
        //
        // ¬(a ∧ b) ∧ ¬(¬a ∧ ¬b) = c
        // (1 - (a * b)) * (1 - ((1 - a) * (1 - b))) = c
        // (1 - ab) * (1 - (1 - a - b + ab)) = c
        // (1 - ab) * (a + b - ab) = c
        // a + b - ab - (a^2)b - (b^2)a + (a^2)(b^2) = c
        // a + b - ab - ab - ab + ab = c
        // a + b - 2ab = c
        // -2a * b = c - a - b
        // 2a * b = a + b - c
        // (a + a) * b = a + b - c
        self.cs.enforce_constraint(
            lc!() + self.variable + self.variable,
            lc!() + b.variable,
            lc!() + self.variable + b.variable - result.variable,
        )?;

        Ok(result)
    }

    /// Performs an AND operation over the two operands, returning
    /// an `AllocatedBit`.
    pub fn and(&self, b: &Self) -> Result<Self, SynthesisError> {
        let result = Self::new_witness_without_booleanity_check(self.cs.clone(), || {
            Ok(self.value()? & b.value()?)
        })?;

        // Constrain (a) * (b) = (c), ensuring c is 1 iff
        // a AND b are both 1.
        self.cs.enforce_constraint(
            lc!() + self.variable,
            lc!() + b.variable,
            lc!() + result.variable,
        )?;

        Ok(result)
    }

    /// Performs an OR operation over the two operands, returning
    /// an `AllocatedBit`.
    pub fn or(&self, b: &Self) -> Result<Self, SynthesisError> {
        let result = Self::new_witness_without_booleanity_check(self.cs.clone(), || {
            Ok(self.value()? | b.value()?)
        })?;

        // Constrain (1 - a) * (1 - b) = (c), ensuring c is 1 iff
        // a and b are both false, and otherwise c is 0.
        self.cs.enforce_constraint(
            lc!() + Variable::One - self.variable,
            lc!() + Variable::One - b.variable,
            lc!() + Variable::One - result.variable,
        )?;

        Ok(result)
    }

    /// Calculates `a AND (NOT b)`.
    pub fn and_not(&self, b: &Self) -> Result<Self, SynthesisError> {
        let result = Self::new_witness_without_booleanity_check(self.cs.clone(), || {
            Ok(self.value()? & !b.value()?)
        })?;

        // Constrain (a) * (1 - b) = (c), ensuring c is 1 iff
        // a is true and b is false, and otherwise c is 0.
        self.cs.enforce_constraint(
            lc!() + self.variable,
            lc!() + Variable::One - b.variable,
            lc!() + result.variable,
        )?;

        Ok(result)
    }

    /// Calculates `(NOT a) AND (NOT b)`.
    pub fn nor(&self, b: &Self) -> Result<Self, SynthesisError> {
        let result = Self::new_witness_without_booleanity_check(self.cs.clone(), || {
            Ok(!(self.value()? | b.value()?))
        })?;

        // Constrain (1 - a) * (1 - b) = (c), ensuring c is 1 iff
        // a and b are both false, and otherwise c is 0.
        self.cs.enforce_constraint(
            lc!() + Variable::One - self.variable,
            lc!() + Variable::One - b.variable,
            lc!() + result.variable,
        )?;

        Ok(result)
    }
}

impl<F: Field> AllocVar<bool, F> for AllocatedBit<F> {
    /// If `self.mode` == `AllocationMode::Constant`, this method simply outputs
    /// a `Boolean::Constant`.
    ///
    /// Otherwise, it produces a new variable of the appropriate type
    /// (instance or witness), with a booleanity check.
    ///
    /// N.B.: we could omit the booleanity check when allocating `self`
    /// as a new public input, but that places an additional burden on
    /// protocol designers. Better safe than sorry!
    fn new_variable<T: Borrow<bool>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        if mode == AllocationMode::Constant {
            let variable = if *f()?.borrow() {
                Variable::One
            } else {
                Variable::Zero
            };
            Ok(Self { variable, cs })
        } else {
            let variable = if mode == AllocationMode::Input {
                cs.new_input_variable(|| f().map(bool_to_field))?
            } else {
                cs.new_witness_variable(|| f().map(bool_to_field))?
            };

            // Constrain: (1 - a) * a = 0
            // This constrains a to be either 0 or 1.

            cs.enforce_named_constraint(
                "Booleanity check",
                lc!() + Variable::One - variable,
                lc!() + variable,
                lc!(),
            )?;

            Ok(Self { variable, cs })
        }
    }
}

impl<F: Field> CondSelectGadget<F> for AllocatedBit<F> {
    fn conditionally_select(
        cond: &Boolean<F>,
        true_val: &Self,
        false_val: &Self,
    ) -> Result<Self, SynthesisError> {
        let res = Boolean::conditionally_select(
            cond,
            &true_val.clone().into(),
            &false_val.clone().into(),
        )?;
        match res {
            Boolean::Is(a) => Ok(a),
            _ => unreachable!("Impossible"),
        }
    }
}

/// This is a boolean value which may be either a constant or
/// an interpretation of an `AllocatedBit`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Boolean<F: Field> {
    /// Existential view of the boolean variable
    Is(AllocatedBit<F>),
    /// Negated view of the boolean variable
    Not(AllocatedBit<F>),
    /// Constant (not an allocated variable)
    Constant(bool),
}

impl<F: Field> R1CSVar<F> for Boolean<F> {
    type Value = bool;

    fn cs(&self) -> Option<ConstraintSystemRef<F>> {
        match self {
            Self::Is(a) | Self::Not(a) => Some(a.cs.clone()),
            _ => None,
        }
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        match self {
            Boolean::Constant(c) => Ok(*c),
            Boolean::Is(ref v) => v.value(),
            Boolean::Not(ref v) => v.value().map(|b| !b),
        }
    }
}

impl<F: Field> Boolean<F> {
    pub fn lc(&self) -> LinearCombination<F> {
        match self {
            Boolean::Constant(false) => lc!(),
            Boolean::Constant(true) => lc!() + Variable::One,
            Boolean::Is(v) => v.variable().into(),
            Boolean::Not(v) => lc!() + Variable::One - v.variable(),
        }
    }

    /// Construct a boolean vector from a vector of u8
    pub fn constant_vec_from_bytes(values: &[u8]) -> Vec<Self> {
        let mut input_bits = vec![];
        for input_byte in values {
            for bit_i in (0..8).rev() {
                input_bits.push(Self::Constant(((input_byte >> bit_i) & 1u8) == 1u8));
            }
        }
        input_bits
    }

    /// Construct a boolean from a known constant
    pub fn constant(b: bool) -> Self {
        Boolean::Constant(b)
    }

    /// Return a negated interpretation of this boolean.
    pub fn not(&self) -> Self {
        match *self {
            Boolean::Constant(c) => Boolean::Constant(!c),
            Boolean::Is(ref v) => Boolean::Not(v.clone()),
            Boolean::Not(ref v) => Boolean::Is(v.clone()),
        }
    }
}
impl<F: Field> Boolean<F> {
    /// Perform XOR over two boolean operands
    pub fn xor<'a>(&'a self, b: &'a Self) -> Result<Self, SynthesisError> {
        use Boolean::*;
        match (self, b) {
            (&Constant(false), x) | (x, &Constant(false)) => Ok(x.clone()),
            (&Constant(true), x) | (x, &Constant(true)) => Ok(x.not()),
            // a XOR (NOT b) = NOT(a XOR b)
            (is @ &Is(_), not @ &Not(_)) | (not @ &Not(_), is @ &Is(_)) => {
                Ok(is.xor(&not.not())?.not())
            }
            // a XOR b = (NOT a) XOR (NOT b)
            (&Is(ref a), &Is(ref b)) | (&Not(ref a), &Not(ref b)) => Ok(Is(a.xor(b)?)),
        }
    }

    /// Perform OR over two boolean operands
    pub fn or<'a>(&'a self, b: &'a Self) -> Result<Self, SynthesisError> {
        use Boolean::*;
        match (self, b) {
            (&Constant(false), x) | (x, &Constant(false)) => Ok(x.clone()),
            (&Constant(true), _) | (_, &Constant(true)) => Ok(Constant(true)),
            // a OR b = NOT ((NOT a) AND b)
            (a @ &Is(_), b @ &Not(_)) | (b @ &Not(_), a @ &Is(_)) | (b @ &Not(_), a @ &Not(_)) => {
                Ok(a.not().and(&b.not())?.not())
            }
            (&Is(ref a), &Is(ref b)) => a.or(b).map(From::from),
        }
    }

    /// Perform AND over two boolean operands
    pub fn and<'a>(&'a self, b: &'a Self) -> Result<Self, SynthesisError> {
        use Boolean::*;
        match (self, b) {
            // false AND x is always false
            (&Constant(false), _) | (_, &Constant(false)) => Ok(Constant(false)),
            // true AND x is always x
            (&Constant(true), x) | (x, &Constant(true)) => Ok(x.clone()),
            // a AND (NOT b)
            (&Is(ref is), &Not(ref not)) | (&Not(ref not), &Is(ref is)) => Ok(Is(is.and_not(not)?)),
            // (NOT a) AND (NOT b) = a NOR b
            (&Not(ref a), &Not(ref b)) => Ok(Is(a.nor(b)?)),
            // a AND b
            (&Is(ref a), &Is(ref b)) => Ok(Is(a.and(b)?)),
        }
    }

    pub fn kary_and(bits: &[Self]) -> Result<Self, SynthesisError> {
        assert!(!bits.is_empty());
        let mut cur: Option<Self> = None;
        for next in bits {
            cur = if let Some(b) = cur {
                Some(b.and(next)?)
            } else {
                Some(next.clone())
            };
        }

        Ok(cur.expect("should not be 0"))
    }

    pub fn kary_or(bits: &[Self]) -> Result<Self, SynthesisError> {
        assert!(!bits.is_empty());
        let mut cur: Option<Self> = None;
        for next in bits {
            cur = if let Some(b) = cur {
                Some(b.or(next)?)
            } else {
                Some(next.clone())
            };
        }

        Ok(cur.expect("should not be 0"))
    }

    pub fn kary_nand(bits: &[Self]) -> Result<Self, SynthesisError> {
        Ok(Self::kary_and(bits)?.not())
    }

    /// Assert that at least one input is false.
    fn enforce_kary_nand(bits: &[Self]) -> Result<(), SynthesisError> {
        use Boolean::*;
        let r = Self::kary_nand(bits)?;
        match r {
            Constant(true) => Ok(()),
            Constant(false) => Err(SynthesisError::AssignmentMissing),
            Is(_) | Not(_) => r.cs().unwrap().enforce_constraint(
                r.lc(),
                lc!() + Variable::One,
                lc!() + Variable::One,
            ),
        }
    }

    /// Asserts that this bit_gadget representation is "in
    /// the field" when interpreted in big endian.
    pub fn enforce_in_field(bits: &[Self]) -> Result<(), SynthesisError> {
        // b = char() - 1
        let mut b = F::characteristic().to_vec();
        assert_eq!(b[0] % 2, 1);
        b[0] -= 1;
        let run = Self::enforce_smaller_or_equal_than(bits, b)?;

        // We should always end in a "run" of zeros, because
        // the characteristic is an odd prime. So, this should
        // be empty.
        assert!(run.is_empty());

        Ok(())
    }

    /// Asserts that this bit_gadget representation is smaller
    /// or equal than the provided element
    pub fn enforce_smaller_or_equal_than(
        bits: &[Self],
        element: impl AsRef<[u64]>,
    ) -> Result<Vec<Self>, SynthesisError> {
        let mut bits_iter = bits.iter();
        let b: &[u64] = element.as_ref();

        // Runs of ones in r
        let mut last_run = Boolean::constant(true);
        let mut current_run = vec![];

        let mut found_one = false;

        let char_num_bits = {
            let mut leading_zeros = 0;
            let mut total_bits = 0;
            let mut found_one = false;
            for b in BitIterator::new(b.clone()) {
                total_bits += 1;
                if !b && !found_one {
                    leading_zeros += 1
                }
                if b {
                    found_one = true;
                }
            }

            total_bits - leading_zeros
        };

        if bits.len() > char_num_bits {
            let num_extra_bits = bits.len() - char_num_bits;
            let mut or_result = Boolean::constant(false);
            for should_be_zero in &bits[0..num_extra_bits] {
                or_result = or_result.or(should_be_zero)?;
                let _ = bits_iter.next().unwrap();
            }
            or_result.enforce_equal(&Boolean::constant(false))?;
        }

        for b in BitIterator::new(b) {
            // Skip over unset bits at the beginning
            found_one |= b;
            if !found_one {
                continue;
            }

            let a = bits_iter.next().unwrap();

            if b {
                // This is part of a run of ones.
                current_run.push(a.clone());
            } else {
                if !current_run.is_empty() {
                    // This is the start of a run of zeros, but we need
                    // to k-ary AND against `last_run` first.

                    current_run.push(last_run.clone());
                    last_run = Self::kary_and(&current_run)?;
                    current_run.truncate(0);
                }

                // If `last_run` is true, `a` must be false, or it would
                // not be in the field.
                //
                // If `last_run` is false, `a` can be true or false.
                //
                // Ergo, at least one of `last_run` and `a` must be false.
                Self::enforce_kary_nand(&[last_run.clone(), a.clone()])?;
            }
        }
        assert!(bits_iter.next().is_none());

        Ok(current_run)
    }

    pub fn select<T: CondSelectGadget<F>>(
        &self,
        first: &T,
        second: &T,
    ) -> Result<T, SynthesisError> {
        T::conditionally_select(&self, first, second)
    }
}

impl<F: Field> From<AllocatedBit<F>> for Boolean<F> {
    fn from(b: AllocatedBit<F>) -> Self {
        Boolean::Is(b)
    }
}

impl<F: Field> AllocVar<bool, F> for Boolean<F> {
    fn new_variable<T: Borrow<bool>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        if mode == AllocationMode::Constant {
            Ok(Boolean::Constant(*f()?.borrow()))
        } else {
            AllocatedBit::new_variable(cs, f, mode).map(Boolean::from)
        }
    }
}

impl<F: Field> EqGadget<F> for Boolean<F> {
    fn is_eq(&self, other: &Self) -> Result<Boolean<F>, SynthesisError> {
        // self | other | XNOR(self, other) | self == other
        // -----|-------|-------------------|--------------
        //   0  |   0   |         1         |      1
        //   0  |   1   |         0         |      0
        //   1  |   0   |         0         |      0
        //   1  |   1   |         1         |      1
        Ok(self.xor(other)?.not())
    }

    fn conditional_enforce_equal(
        &self,
        other: &Self,
        condition: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        use Boolean::*;
        let one = Variable::One;
        let difference = match (self, other) {
            // 1 == 1; 0 == 0
            (Constant(true), Constant(true)) | (Constant(false), Constant(false)) => return Ok(()),
            // false != true
            (Constant(_), Constant(_)) => return Err(SynthesisError::AssignmentMissing),
            // 1 - a
            (Constant(true), Is(a)) | (Is(a), Constant(true)) => lc!() + one - a.variable(),
            // a - 0 = a
            (Constant(false), Is(a)) | (Is(a), Constant(false)) => lc!() + a.variable(),
            // 1 - !a = 1 - (1 - a) = a
            (Constant(true), Not(a)) | (Not(a), Constant(true)) => lc!() + a.variable(),
            // !a - 0 = !a = 1 - a
            (Constant(false), Not(a)) | (Not(a), Constant(false)) => lc!() + one - a.variable(),
            // b - a,
            (Is(a), Is(b)) => lc!() + b.variable() - a.variable(),
            // !b - a = (1 - b) - a
            (Is(a), Not(b)) | (Not(b), Is(a)) => lc!() + one - b.variable() - a.variable(),
            // !b - !a = (1 - b) - (1 - a) = a - b,
            (Not(a), Not(b)) => lc!() + a.variable() - b.variable(),
        };

        if condition != &Constant(false) {
            let cs = self.cs().or(other.cs()).or(condition.cs()).unwrap();
            cs.enforce_constraint(lc!() + difference, condition.lc(), lc!())?;
        }
        Ok(())
    }

    fn conditional_enforce_not_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        use Boolean::*;
        let one = Variable::One;
        let difference = match (self, other) {
            // 1 != 0; 0 != 1
            (Constant(true), Constant(false)) | (Constant(false), Constant(true)) => return Ok(()),
            // false == false and true == true
            (Constant(_), Constant(_)) => return Err(SynthesisError::AssignmentMissing),
            // 1 - a
            (Constant(true), Is(a)) | (Is(a), Constant(true)) => lc!() + one - a.variable(),
            // a - 0 = a
            (Constant(false), Is(a)) | (Is(a), Constant(false)) => lc!() + a.variable(),
            // 1 - !a = 1 - (1 - a) = a
            (Constant(true), Not(a)) | (Not(a), Constant(true)) => lc!() + a.variable(),
            // !a - 0 = !a = 1 - a
            (Constant(false), Not(a)) | (Not(a), Constant(false)) => lc!() + one - a.variable(),
            // b - a,
            (Is(a), Is(b)) => lc!() + b.variable() - a.variable(),
            // !b - a = (1 - b) - a
            (Is(a), Not(b)) | (Not(b), Is(a)) => lc!() + one - b.variable() - a.variable(),
            // !b - !a = (1 - b) - (1 - a) = a - b,
            (Not(a), Not(b)) => lc!() + a.variable() - b.variable(),
        };

        if should_enforce != &Constant(false) {
            let cs = self
                .cs()
                .or(other.cs())
                .or(should_enforce.cs())
                .ok_or(SynthesisError::UnconstrainedVariable)?;
            cs.enforce_constraint(difference, should_enforce.lc(), should_enforce.lc())?;
        }
        Ok(())
    }
}

impl<F: Field> ToBytesGadget<F> for Boolean<F> {
    /// Outputs `1u8` if `self` is true, and `0u8` otherwise.
    fn to_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        let mut bits = vec![Boolean::constant(false); 7];
        bits.push(self.clone());
        bits.reverse();
        let value = self.value().map(|val| val as u8).ok();
        let byte = UInt8 { bits, value };
        Ok(vec![byte])
    }
}

impl<F: Field> CondSelectGadget<F> for Boolean<F> {
    fn conditionally_select(
        cond: &Boolean<F>,
        true_val: &Self,
        false_val: &Self,
    ) -> Result<Self, SynthesisError> {
        use Boolean::*;
        match cond {
            Constant(true) => Ok(true_val.clone()),
            Constant(false) => Ok(false_val.clone()),
            cond @ Not(_) => Self::conditionally_select(&cond.not(), false_val, true_val),
            cond @ Is(_) => match (true_val, false_val) {
                (x, &Constant(false)) => cond.and(x),
                (&Constant(false), x) => cond.not().and(x),
                (&Constant(true), x) => cond.or(x),
                (x, &Constant(true)) => cond.not().or(x),
                (a, b) => {
                    let cs = cond.cs().unwrap();
                    let result: Boolean<F> =
                        AllocatedBit::new_witness_without_booleanity_check(cs.clone(), || {
                            let cond = cond.value()?;
                            Ok(if cond { a.value()? } else { b.value()? })
                        })?
                        .into();
                    //
                    // a = self; b = other; c = cond;
                    //
                    // r = c * a + (1  - c) * b
                    // r = b + c * (a - b)
                    // c * (a - b) = r - b
                    //
                    // If a, b, cond are all boolean, so is r.
                    //
                    // self | other | cond | result
                    // -----|-------|----------------
                    //   0  |   0   |   1  |    0
                    //   0  |   1   |   1  |    0
                    //   1  |   0   |   1  |    1
                    //   1  |   1   |   1  |    1
                    //   0  |   0   |   0  |    0
                    //   0  |   1   |   0  |    1
                    //   1  |   0   |   0  |    0
                    //   1  |   1   |   0  |    1
                    cs.enforce_constraint(
                        cond.lc(),
                        lc!() + a.lc() - b.lc(),
                        lc!() + result.lc() - b.lc(),
                    )?;

                    Ok(result)
                }
            },
        }
    }
}

#[cfg(test)]
mod test {
    use super::{AllocatedBit, Boolean};
    use crate::prelude::*;
    use algebra::{bls12_381::Fr, BitIterator, Field, One, PrimeField, UniformRand, Zero};
    use r1cs_core::{ConstraintSystem, Namespace, SynthesisError};
    use rand::SeedableRng;
    use rand_xorshift::XorShiftRng;

    #[test]
    fn test_boolean_to_byte() -> Result<(), SynthesisError> {
        for val in [true, false].iter() {
            let cs = ConstraintSystem::<Fr>::new_ref();
            let a = Boolean::new_witness(cs.clone(), || Ok(*val))?;
            let bytes = a.to_bytes()?;
            assert_eq!(bytes.len(), 1);
            let byte = &bytes[0];
            assert_eq!(byte.value()?, *val as u8);

            for (i, bit) in byte.bits.iter().enumerate() {
                assert_eq!(bit.value()?, (byte.value()? >> i) & 1 == 1);
            }
        }
        Ok(())
    }

    #[test]
    fn test_xor() -> Result<(), SynthesisError> {
        for a_val in [false, true].iter().copied() {
            for b_val in [false, true].iter().copied() {
                let cs = ConstraintSystem::<Fr>::new_ref();
                let a = AllocatedBit::new_witness(cs.ns("a"), || Ok(a_val))?;
                let b = AllocatedBit::new_witness(cs.ns("b"), || Ok(b_val))?;
                let c = AllocatedBit::xor(&a, &b)?;
                assert_eq!(c.value()?, a_val ^ b_val);

                assert!(cs.is_satisfied().unwrap());
                assert_eq!(a.value()?, (a_val));
                assert_eq!(b.value()?, (b_val));
                assert_eq!(c.value()?, (a_val ^ b_val));
            }
        }
        Ok(())
    }

    #[test]
    fn test_or() -> Result<(), SynthesisError> {
        for a_val in [false, true].iter().copied() {
            for b_val in [false, true].iter().copied() {
                let cs = ConstraintSystem::<Fr>::new_ref();
                let a = AllocatedBit::new_witness(cs.ns("a"), || Ok(a_val))?;
                let b = AllocatedBit::new_witness(cs.ns("b"), || Ok(b_val))?;
                let c = AllocatedBit::or(&a, &b)?;
                assert_eq!(c.value()?, a_val | b_val);

                assert!(cs.is_satisfied().unwrap());
                assert_eq!(a.value()?, (a_val));
                assert_eq!(b.value()?, (b_val));
                assert_eq!(c.value()?, (a_val | b_val));
            }
        }
        Ok(())
    }

    #[test]
    fn test_and() -> Result<(), SynthesisError> {
        for a_val in [false, true].iter().copied() {
            for b_val in [false, true].iter().copied() {
                let cs = ConstraintSystem::<Fr>::new_ref();
                let a = AllocatedBit::new_witness(cs.ns("a"), || Ok(a_val))?;
                let b = AllocatedBit::new_witness(cs.ns("b"), || Ok(b_val))?;
                let c = AllocatedBit::and(&a, &b)?;
                assert_eq!(c.value()?, a_val & b_val);

                assert!(cs.is_satisfied().unwrap());
                assert_eq!(a.value()?, (a_val));
                assert_eq!(b.value()?, (b_val));
                assert_eq!(c.value()?, (a_val & b_val));
            }
        }
        Ok(())
    }

    #[test]
    fn test_and_not() -> Result<(), SynthesisError> {
        for a_val in [false, true].iter().copied() {
            for b_val in [false, true].iter().copied() {
                let cs = ConstraintSystem::<Fr>::new_ref();
                let a = AllocatedBit::new_witness(cs.ns("a"), || Ok(a_val))?;
                let b = AllocatedBit::new_witness(cs.ns("b"), || Ok(b_val))?;
                let c = AllocatedBit::and_not(&a, &b)?;
                assert_eq!(c.value()?, a_val & !b_val);

                assert!(cs.is_satisfied().unwrap());
                assert_eq!(a.value()?, (a_val));
                assert_eq!(b.value()?, (b_val));
                assert_eq!(c.value()?, (a_val & !b_val));
            }
        }
        Ok(())
    }

    #[test]
    fn test_nor() -> Result<(), SynthesisError> {
        for a_val in [false, true].iter().copied() {
            for b_val in [false, true].iter().copied() {
                let cs = ConstraintSystem::<Fr>::new_ref();
                let a = AllocatedBit::new_witness(cs.ns("a"), || Ok(a_val))?;
                let b = AllocatedBit::new_witness(cs.ns("b"), || Ok(b_val))?;
                let c = AllocatedBit::nor(&a, &b)?;
                assert_eq!(c.value()?, !a_val & !b_val);

                assert!(cs.is_satisfied().unwrap());
                assert_eq!(a.value()?, (a_val));
                assert_eq!(b.value()?, (b_val));
                assert_eq!(c.value()?, (!a_val & !b_val));
            }
        }
        Ok(())
    }

    #[test]
    fn test_enforce_equal() -> Result<(), SynthesisError> {
        for a_bool in [false, true].iter().cloned() {
            for b_bool in [false, true].iter().cloned() {
                for a_neg in [false, true].iter().cloned() {
                    for b_neg in [false, true].iter().cloned() {
                        let cs = ConstraintSystem::<Fr>::new_ref();

                        let mut a = Boolean::new_witness(cs.ns("a"), || Ok(a_bool))?;
                        let mut b = Boolean::new_witness(cs.ns("b"), || Ok(b_bool))?;

                        if a_neg {
                            a = a.not();
                        }
                        if b_neg {
                            b = b.not();
                        }

                        a.enforce_equal(&b)?;

                        assert_eq!(
                            cs.is_satisfied().unwrap(),
                            (a_bool ^ a_neg) == (b_bool ^ b_neg)
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_conditional_enforce_equal() -> Result<(), SynthesisError> {
        for a_bool in [false, true].iter().cloned() {
            for b_bool in [false, true].iter().cloned() {
                for a_neg in [false, true].iter().cloned() {
                    for b_neg in [false, true].iter().cloned() {
                        let cs = ConstraintSystem::<Fr>::new_ref();

                        // First test if constraint system is satisfied
                        // when we do want to enforce the condition.
                        let mut a = Boolean::new_witness(cs.clone(), || Ok(a_bool))?;
                        let mut b = Boolean::new_witness(cs.clone(), || Ok(b_bool))?;

                        if a_neg {
                            a = a.not();
                        }
                        if b_neg {
                            b = b.not();
                        }

                        a.conditional_enforce_equal(&b, &Boolean::constant(true))?;

                        assert_eq!(
                            cs.is_satisfied().unwrap(),
                            (a_bool ^ a_neg) == (b_bool ^ b_neg)
                        );

                        // Now test if constraint system is satisfied even
                        // when we don't want to enforce the condition.
                        let cs = ConstraintSystem::<Fr>::new_ref();

                        let mut a = Boolean::new_witness(cs.ns("a"), || Ok(a_bool))?;
                        let mut b = Boolean::new_witness(cs.ns("b"), || Ok(b_bool))?;

                        if a_neg {
                            a = a.not();
                        }
                        if b_neg {
                            b = b.not();
                        }

                        let false_cond = Boolean::new_witness(cs.ns("cond"), || Ok(false))?;
                        a.conditional_enforce_equal(&b, &false_cond)?;

                        assert!(cs.is_satisfied().unwrap());
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_boolean_negation() -> Result<(), SynthesisError> {
        let cs = ConstraintSystem::<Fr>::new_ref();

        let mut b = Boolean::new_witness(cs.clone(), || Ok(true))?;
        assert!(matches!(b, Boolean::Is(_)));

        b = b.not();
        assert!(matches!(b, Boolean::Not(_)));

        b = b.not();
        assert!(matches!(b, Boolean::Is(_)));

        b = Boolean::Constant(true);
        assert!(matches!(b, Boolean::Constant(true)));

        b = b.not();
        assert!(matches!(b, Boolean::Constant(false)));

        b = b.not();
        assert!(matches!(b, Boolean::Constant(true)));
        Ok(())
    }

    #[derive(Eq, PartialEq, Copy, Clone, Debug)]
    enum OpType {
        True,
        False,
        AllocatedTrue,
        AllocatedFalse,
        NegatedAllocatedTrue,
        NegatedAllocatedFalse,
    }

    const VARIANTS: [OpType; 6] = [
        OpType::True,
        OpType::False,
        OpType::AllocatedTrue,
        OpType::AllocatedFalse,
        OpType::NegatedAllocatedTrue,
        OpType::NegatedAllocatedFalse,
    ];

    fn construct<F: Field>(
        cs: impl Into<Namespace<F>>,
        operand: OpType,
        name: &'static str,
    ) -> Result<Boolean<F>, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs().ns(name);

        let b = match operand {
            OpType::True => Boolean::constant(true),
            OpType::False => Boolean::constant(false),
            OpType::AllocatedTrue => Boolean::new_witness(cs, || Ok(true))?,
            OpType::AllocatedFalse => Boolean::new_witness(cs, || Ok(false))?,
            OpType::NegatedAllocatedTrue => Boolean::new_witness(cs, || Ok(true))?.not(),
            OpType::NegatedAllocatedFalse => Boolean::new_witness(cs, || Ok(false))?.not(),
        };
        Ok(b)
    }

    #[test]
    fn test_boolean_xor() -> Result<(), SynthesisError> {
        for first_operand in VARIANTS.iter().cloned() {
            for second_operand in VARIANTS.iter().cloned() {
                let cs = ConstraintSystem::<Fr>::new_ref();

                let a = construct(cs.clone(), first_operand, "a")?;
                let b = construct(cs.clone(), second_operand, "b")?;
                let c = Boolean::xor(&a, &b)?;

                assert!(cs.is_satisfied().unwrap());

                match (first_operand, second_operand, c) {
                    (OpType::True, OpType::True, Boolean::Constant(false)) => (),
                    (OpType::True, OpType::False, Boolean::Constant(true)) => (),
                    (OpType::True, OpType::AllocatedTrue, Boolean::Not(_)) => (),
                    (OpType::True, OpType::AllocatedFalse, Boolean::Not(_)) => (),
                    (OpType::True, OpType::NegatedAllocatedTrue, Boolean::Is(_)) => (),
                    (OpType::True, OpType::NegatedAllocatedFalse, Boolean::Is(_)) => (),

                    (OpType::False, OpType::True, Boolean::Constant(true)) => (),
                    (OpType::False, OpType::False, Boolean::Constant(false)) => (),
                    (OpType::False, OpType::AllocatedTrue, Boolean::Is(_)) => (),
                    (OpType::False, OpType::AllocatedFalse, Boolean::Is(_)) => (),
                    (OpType::False, OpType::NegatedAllocatedTrue, Boolean::Not(_)) => (),
                    (OpType::False, OpType::NegatedAllocatedFalse, Boolean::Not(_)) => (),

                    (OpType::AllocatedTrue, OpType::True, Boolean::Not(_)) => (),
                    (OpType::AllocatedTrue, OpType::False, Boolean::Is(_)) => (),
                    (OpType::AllocatedTrue, OpType::AllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedTrue, OpType::AllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(v.value(), Ok(true));
                    }
                    (OpType::AllocatedTrue, OpType::NegatedAllocatedTrue, Boolean::Not(ref v)) => {
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedTrue, OpType::NegatedAllocatedFalse, Boolean::Not(ref v)) => {
                        assert_eq!(v.value(), Ok(true));
                    }
                    (OpType::AllocatedFalse, OpType::True, Boolean::Not(_)) => (),
                    (OpType::AllocatedFalse, OpType::False, Boolean::Is(_)) => (),
                    (OpType::AllocatedFalse, OpType::AllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }
                    (OpType::AllocatedFalse, OpType::AllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedFalse, OpType::NegatedAllocatedTrue, Boolean::Not(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }
                    (
                        OpType::AllocatedFalse,
                        OpType::NegatedAllocatedFalse,
                        Boolean::Not(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }

                    (OpType::NegatedAllocatedTrue, OpType::True, Boolean::Is(_)) => (),
                    (OpType::NegatedAllocatedTrue, OpType::False, Boolean::Not(_)) => (),
                    (OpType::NegatedAllocatedTrue, OpType::AllocatedTrue, Boolean::Not(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::NegatedAllocatedTrue, OpType::AllocatedFalse, Boolean::Not(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }
                    (
                        OpType::NegatedAllocatedTrue,
                        OpType::NegatedAllocatedTrue,
                        Boolean::Is(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (
                        OpType::NegatedAllocatedTrue,
                        OpType::NegatedAllocatedFalse,
                        Boolean::Is(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }

                    (OpType::NegatedAllocatedFalse, OpType::True, Boolean::Is(_)) => (),
                    (OpType::NegatedAllocatedFalse, OpType::False, Boolean::Not(_)) => (),
                    (OpType::NegatedAllocatedFalse, OpType::AllocatedTrue, Boolean::Not(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }
                    (
                        OpType::NegatedAllocatedFalse,
                        OpType::AllocatedFalse,
                        Boolean::Not(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (
                        OpType::NegatedAllocatedFalse,
                        OpType::NegatedAllocatedTrue,
                        Boolean::Is(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }
                    (
                        OpType::NegatedAllocatedFalse,
                        OpType::NegatedAllocatedFalse,
                        Boolean::Is(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }

                    _ => unreachable!(),
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_boolean_cond_select() -> Result<(), SynthesisError> {
        for condition in VARIANTS.iter().cloned() {
            for first_operand in VARIANTS.iter().cloned() {
                for second_operand in VARIANTS.iter().cloned() {
                    let cs = ConstraintSystem::<Fr>::new_ref();

                    let cond = construct(cs.clone(), condition, "cond")?;
                    let a = construct(cs.clone(), first_operand, "a")?;
                    let b = construct(cs.clone(), second_operand, "b")?;
                    let c = cond.select(&a, &b)?;

                    assert!(
                        cs.is_satisfied().unwrap(),
                        "failed with operands: cond: {:?}, a: {:?}, b: {:?}",
                        condition,
                        first_operand,
                        second_operand,
                    );
                    assert_eq!(
                        c.value()?,
                        if cond.value()? {
                            a.value()?
                        } else {
                            b.value()?
                        }
                    );
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_boolean_or() -> Result<(), SynthesisError> {
        for first_operand in VARIANTS.iter().cloned() {
            for second_operand in VARIANTS.iter().cloned() {
                let cs = ConstraintSystem::<Fr>::new_ref();

                let a = construct(cs.clone(), first_operand, "a")?;
                let b = construct(cs.clone(), second_operand, "b")?;
                let c = a.or(&b)?;

                assert!(cs.is_satisfied().unwrap());

                match (first_operand, second_operand, c.clone()) {
                    (OpType::True, OpType::True, Boolean::Constant(true)) => (),
                    (OpType::True, OpType::False, Boolean::Constant(true)) => (),
                    (OpType::True, OpType::AllocatedTrue, Boolean::Constant(true)) => (),
                    (OpType::True, OpType::AllocatedFalse, Boolean::Constant(true)) => (),
                    (OpType::True, OpType::NegatedAllocatedTrue, Boolean::Constant(true)) => (),
                    (OpType::True, OpType::NegatedAllocatedFalse, Boolean::Constant(true)) => (),

                    (OpType::False, OpType::True, Boolean::Constant(true)) => (),
                    (OpType::False, OpType::False, Boolean::Constant(false)) => (),
                    (OpType::False, OpType::AllocatedTrue, Boolean::Is(_)) => (),
                    (OpType::False, OpType::AllocatedFalse, Boolean::Is(_)) => (),
                    (OpType::False, OpType::NegatedAllocatedTrue, Boolean::Not(_)) => (),
                    (OpType::False, OpType::NegatedAllocatedFalse, Boolean::Not(_)) => (),

                    (OpType::AllocatedTrue, OpType::True, Boolean::Constant(true)) => (),
                    (OpType::AllocatedTrue, OpType::False, Boolean::Is(_)) => (),
                    (OpType::AllocatedTrue, OpType::AllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(v.value(), Ok(true));
                    }
                    (OpType::AllocatedTrue, OpType::AllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(v.value(), Ok(true));
                    }
                    (OpType::AllocatedTrue, OpType::NegatedAllocatedTrue, Boolean::Not(ref v)) => {
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedTrue, OpType::NegatedAllocatedFalse, Boolean::Not(ref v)) => {
                        assert_eq!(v.value(), Ok(false));
                    }

                    (OpType::AllocatedFalse, OpType::True, Boolean::Constant(true)) => (),
                    (OpType::AllocatedFalse, OpType::False, Boolean::Is(_)) => (),
                    (OpType::AllocatedFalse, OpType::AllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(v.value(), Ok(true));
                    }
                    (OpType::AllocatedFalse, OpType::AllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedFalse, OpType::NegatedAllocatedTrue, Boolean::Not(ref v)) => {
                        assert_eq!(v.value(), Ok(true));
                    }
                    (
                        OpType::AllocatedFalse,
                        OpType::NegatedAllocatedFalse,
                        Boolean::Not(ref v),
                    ) => {
                        assert_eq!(v.value(), Ok(false));
                    }

                    (OpType::NegatedAllocatedTrue, OpType::True, Boolean::Constant(true)) => (),
                    (OpType::NegatedAllocatedTrue, OpType::False, Boolean::Not(_)) => (),
                    (OpType::NegatedAllocatedTrue, OpType::AllocatedTrue, Boolean::Not(ref v)) => {
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::NegatedAllocatedTrue, OpType::AllocatedFalse, Boolean::Not(ref v)) => {
                        assert_eq!(v.value(), Ok(true));
                    }
                    (
                        OpType::NegatedAllocatedTrue,
                        OpType::NegatedAllocatedTrue,
                        Boolean::Not(ref v),
                    ) => {
                        assert_eq!(v.value(), Ok(true));
                    }
                    (
                        OpType::NegatedAllocatedTrue,
                        OpType::NegatedAllocatedFalse,
                        Boolean::Not(ref v),
                    ) => {
                        assert_eq!(v.value(), Ok(false));
                    }

                    (OpType::NegatedAllocatedFalse, OpType::True, Boolean::Constant(true)) => (),
                    (OpType::NegatedAllocatedFalse, OpType::False, Boolean::Not(_)) => (),
                    (OpType::NegatedAllocatedFalse, OpType::AllocatedTrue, Boolean::Not(ref v)) => {
                        assert_eq!(v.value(), Ok(false));
                    }
                    (
                        OpType::NegatedAllocatedFalse,
                        OpType::AllocatedFalse,
                        Boolean::Not(ref v),
                    ) => {
                        assert_eq!(v.value(), Ok(false));
                    }
                    (
                        OpType::NegatedAllocatedFalse,
                        OpType::NegatedAllocatedTrue,
                        Boolean::Not(ref v),
                    ) => {
                        assert_eq!(v.value(), Ok(false));
                    }
                    (
                        OpType::NegatedAllocatedFalse,
                        OpType::NegatedAllocatedFalse,
                        Boolean::Not(ref v),
                    ) => {
                        assert_eq!(v.value(), Ok(false));
                    }

                    _ => panic!(
                        "this should never be encountered, in case: (a = {:?}, b = {:?}, c = {:?})",
                        a, b, c
                    ),
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_boolean_and() -> Result<(), SynthesisError> {
        for first_operand in VARIANTS.iter().cloned() {
            for second_operand in VARIANTS.iter().cloned() {
                let cs = ConstraintSystem::<Fr>::new_ref();

                let a = construct(cs.clone(), first_operand, "a")?;
                let b = construct(cs.clone(), second_operand, "b")?;
                let c = a.and(&b)?;

                assert!(cs.is_satisfied().unwrap());

                match (first_operand, second_operand, c) {
                    (OpType::True, OpType::True, Boolean::Constant(true)) => (),
                    (OpType::True, OpType::False, Boolean::Constant(false)) => (),
                    (OpType::True, OpType::AllocatedTrue, Boolean::Is(_)) => (),
                    (OpType::True, OpType::AllocatedFalse, Boolean::Is(_)) => (),
                    (OpType::True, OpType::NegatedAllocatedTrue, Boolean::Not(_)) => (),
                    (OpType::True, OpType::NegatedAllocatedFalse, Boolean::Not(_)) => (),

                    (OpType::False, OpType::True, Boolean::Constant(false)) => (),
                    (OpType::False, OpType::False, Boolean::Constant(false)) => (),
                    (OpType::False, OpType::AllocatedTrue, Boolean::Constant(false)) => (),
                    (OpType::False, OpType::AllocatedFalse, Boolean::Constant(false)) => (),
                    (OpType::False, OpType::NegatedAllocatedTrue, Boolean::Constant(false)) => (),
                    (OpType::False, OpType::NegatedAllocatedFalse, Boolean::Constant(false)) => (),

                    (OpType::AllocatedTrue, OpType::True, Boolean::Is(_)) => (),
                    (OpType::AllocatedTrue, OpType::False, Boolean::Constant(false)) => (),
                    (OpType::AllocatedTrue, OpType::AllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }
                    (OpType::AllocatedTrue, OpType::AllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedTrue, OpType::NegatedAllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedTrue, OpType::NegatedAllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }

                    (OpType::AllocatedFalse, OpType::True, Boolean::Is(_)) => (),
                    (OpType::AllocatedFalse, OpType::False, Boolean::Constant(false)) => (),
                    (OpType::AllocatedFalse, OpType::AllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedFalse, OpType::AllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedFalse, OpType::NegatedAllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::AllocatedFalse, OpType::NegatedAllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }

                    (OpType::NegatedAllocatedTrue, OpType::True, Boolean::Not(_)) => (),
                    (OpType::NegatedAllocatedTrue, OpType::False, Boolean::Constant(false)) => (),
                    (OpType::NegatedAllocatedTrue, OpType::AllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (OpType::NegatedAllocatedTrue, OpType::AllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (
                        OpType::NegatedAllocatedTrue,
                        OpType::NegatedAllocatedTrue,
                        Boolean::Is(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (
                        OpType::NegatedAllocatedTrue,
                        OpType::NegatedAllocatedFalse,
                        Boolean::Is(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }

                    (OpType::NegatedAllocatedFalse, OpType::True, Boolean::Not(_)) => (),
                    (OpType::NegatedAllocatedFalse, OpType::False, Boolean::Constant(false)) => (),
                    (OpType::NegatedAllocatedFalse, OpType::AllocatedTrue, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }
                    (OpType::NegatedAllocatedFalse, OpType::AllocatedFalse, Boolean::Is(ref v)) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (
                        OpType::NegatedAllocatedFalse,
                        OpType::NegatedAllocatedTrue,
                        Boolean::Is(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::zero());
                        assert_eq!(v.value(), Ok(false));
                    }
                    (
                        OpType::NegatedAllocatedFalse,
                        OpType::NegatedAllocatedFalse,
                        Boolean::Is(ref v),
                    ) => {
                        assert_eq!(cs.assigned_value(v.variable()).unwrap(), Fr::one());
                        assert_eq!(v.value(), Ok(true));
                    }

                    _ => {
                        panic!(
                            "unexpected behavior at {:?} AND {:?}",
                            first_operand, second_operand
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_enforce_in_field() -> Result<(), SynthesisError> {
        {
            let cs = ConstraintSystem::<Fr>::new_ref();

            let mut bits = vec![];
            for b in BitIterator::new(Fr::characteristic()).skip(1) {
                bits.push(Boolean::new_witness(cs.clone(), || Ok(b))?);
            }

            Boolean::enforce_in_field(&bits)?;

            assert!(!cs.is_satisfied().unwrap());
        }

        let mut rng = XorShiftRng::seed_from_u64(1231275789u64);

        for _ in 0..1000 {
            let r = Fr::rand(&mut rng);
            let cs = ConstraintSystem::<Fr>::new_ref();

            let mut bits = vec![];
            for b in BitIterator::new(r.into_repr()).skip(1) {
                bits.push(Boolean::new_witness(cs.clone(), || Ok(b))?);
            }

            Boolean::enforce_in_field(&bits)?;

            assert!(cs.is_satisfied().unwrap());
        }
        Ok(())
    }

    #[test]
    fn test_enforce_nand() -> Result<(), SynthesisError> {
        {
            let cs = ConstraintSystem::<Fr>::new_ref();

            assert!(
                Boolean::enforce_kary_nand(&[Boolean::new_constant(cs.clone(), false)?]).is_ok()
            );
            assert!(
                Boolean::enforce_kary_nand(&[Boolean::new_constant(cs.clone(), true)?]).is_err()
            );
        }

        for i in 1..5 {
            // with every possible assignment for them
            for mut b in 0..(1 << i) {
                // with every possible negation
                for mut n in 0..(1 << i) {
                    let cs = ConstraintSystem::<Fr>::new_ref();

                    let mut expected = true;

                    let mut bits = vec![];
                    for _ in 0..i {
                        expected &= b & 1 == 1;

                        let bit = if n & 1 == 1 {
                            Boolean::new_witness(cs.clone(), || Ok(b & 1 == 1))?
                        } else {
                            Boolean::new_witness(cs.clone(), || Ok(b & 1 == 0))?.not()
                        };
                        bits.push(bit);

                        b >>= 1;
                        n >>= 1;
                    }

                    let expected = !expected;

                    Boolean::enforce_kary_nand(&bits)?;

                    if expected {
                        assert!(cs.is_satisfied().unwrap());
                    } else {
                        assert!(!cs.is_satisfied().unwrap());
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn test_kary_and() -> Result<(), SynthesisError> {
        // test different numbers of operands
        for i in 1..15 {
            // with every possible assignment for them
            for mut b in 0..(1 << i) {
                let cs = ConstraintSystem::<Fr>::new_ref();

                let mut expected = true;

                let mut bits = vec![];
                for _ in 0..i {
                    expected &= b & 1 == 1;
                    bits.push(Boolean::new_witness(cs.clone(), || Ok(b & 1 == 1))?);
                    b >>= 1;
                }

                let r = Boolean::kary_and(&bits)?;

                assert!(cs.is_satisfied().unwrap());

                if let Boolean::Is(ref r) = r {
                    assert_eq!(r.value()?, expected);
                }
            }
        }
        Ok(())
    }
}
