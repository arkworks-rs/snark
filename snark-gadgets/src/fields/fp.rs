use algebra::{bytes::ToBytes, Field, FpParameters, PairingEngine, PrimeField};
use snark::{ConstraintSystem, LinearCombination, SynthesisError};

use std::{borrow::Borrow, fmt::Debug};

use super::FieldGadget;
use crate::{
    boolean::{AllocatedBit, Boolean},
    uint8::UInt8,
};

use crate::utils::{
    AllocGadget, CondSelectGadget, ConditionalEqGadget, EqGadget, NEqGadget, ToBitsGadget,
    ToBytesGadget, TwoBitLookupGadget,
};

use crate::{
    Assignment,
    ConstraintVar::{self, *},
};

#[derive(Derivative)]
#[derivative(Debug(bound = "E: Debug"))]
pub struct FpGadget<E: PairingEngine> {
    pub value:    Option<E::Fr>,
    pub variable: ConstraintVar<E>,
}

impl<E: PairingEngine> FpGadget<E> {
    #[inline]
    pub fn from<CS: ConstraintSystem<E>>(mut cs: CS, value: &E::Fr) -> Self {
        Self::alloc(cs.ns(|| "from"), || Ok(*value)).unwrap()
    }
}

impl<E: PairingEngine> FieldGadget<E::Fr, E> for FpGadget<E> {
    type Variable = ConstraintVar<E>;

    #[inline]
    fn get_value(&self) -> Option<E::Fr> {
        self.value
    }

    #[inline]
    fn get_variable(&self) -> Self::Variable {
        self.variable.clone()
    }

    #[inline]
    fn zero<CS: ConstraintSystem<E>>(_cs: CS) -> Result<Self, SynthesisError> {
        let value = Some(E::Fr::zero());
        Ok(FpGadget {
            value,
            variable: ConstraintVar::zero(),
        })
    }

    #[inline]
    fn one<CS: ConstraintSystem<E>>(_cs: CS) -> Result<Self, SynthesisError> {
        let value = Some(E::Fr::one());
        Ok(FpGadget {
            value,
            variable: CS::one().into(),
        })
    }

    #[inline]
    fn add<CS: ConstraintSystem<E>>(
        &self,
        mut _cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let value = match (self.value, other.value) {
            (Some(val1), Some(val2)) => Some(val1 + &val2),
            (..) => None,
        };

        Ok(FpGadget {
            value,
            variable: &self.variable + &other.variable,
        })
    }

    fn double<CS: ConstraintSystem<E>>(&self, _cs: CS) -> Result<Self, SynthesisError> {
        let value = self.value.map(|val| val.double());
        let mut variable = self.variable.clone();
        variable.double_in_place();
        Ok(FpGadget { value, variable })
    }

    fn double_in_place<CS: ConstraintSystem<E>>(
        &mut self,
        _cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.value.as_mut().map(|val| val.double_in_place());
        self.variable.double_in_place();
        Ok(self)
    }

    #[inline]
    fn sub<CS: ConstraintSystem<E>>(
        &self,
        mut _cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let value = match (self.value, other.value) {
            (Some(val1), Some(val2)) => Some(val1 - &val2),
            (..) => None,
        };

        Ok(FpGadget {
            value,
            variable: &self.variable - &other.variable,
        })
    }

    #[inline]
    fn negate<CS: ConstraintSystem<E>>(&self, cs: CS) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.negate_in_place(cs)?;
        Ok(result)
    }

    #[inline]
    fn negate_in_place<CS: ConstraintSystem<E>>(
        &mut self,
        _cs: CS,
    ) -> Result<&mut Self, SynthesisError> {
        self.value.as_mut().map(|val| *val = -(*val));
        self.variable.negate_in_place();
        Ok(self)
    }

    #[inline]
    fn mul<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let product = Self::alloc(cs.ns(|| "mul"), || {
            Ok(self.value.get()? * &other.value.get()?)
        })?;
        cs.enforce(
            || "mul_constraint",
            |lc| &self.variable + lc,
            |lc| &other.variable + lc,
            |lc| &product.variable + lc,
        );
        Ok(product)
    }

    #[inline]
    fn add_constant<CS: ConstraintSystem<E>>(
        &self,
        _cs: CS,
        other: &E::Fr,
    ) -> Result<Self, SynthesisError> {
        let value = self.value.map(|val| val + other);
        Ok(FpGadget {
            value,
            variable: self.variable.clone() + (*other, CS::one()),
        })
    }

    #[inline]
    fn add_constant_in_place<CS: ConstraintSystem<E>>(
        &mut self,
        _cs: CS,
        other: &E::Fr,
    ) -> Result<&mut Self, SynthesisError> {
        self.value.as_mut().map(|val| *val += other);
        self.variable += (*other, CS::one());
        Ok(self)
    }

    #[inline]
    fn mul_by_constant<CS: ConstraintSystem<E>>(
        &self,
        cs: CS,
        other: &E::Fr,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.mul_by_constant_in_place(cs, other)?;
        Ok(result)
    }

    #[inline]
    fn mul_by_constant_in_place<CS: ConstraintSystem<E>>(
        &mut self,
        mut _cs: CS,
        other: &E::Fr,
    ) -> Result<&mut Self, SynthesisError> {
        self.value.as_mut().map(|val| *val *= other);
        self.variable *= *other;
        Ok(self)
    }

    #[inline]
    fn inverse<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Self, SynthesisError> {
        let inverse = Self::alloc(cs.ns(|| "inverse"), || {
            let result = self.value.get()?;
            let inv = result.inverse().expect("Inverse doesn't exist!");
            Ok(inv)
        })?;

        let one = CS::one();
        cs.enforce(
            || "inv_constraint",
            |lc| &self.variable + lc,
            |lc| &inverse.variable + lc,
            |lc| lc + one,
        );
        Ok(inverse)
    }

    fn frobenius_map<CS: ConstraintSystem<E>>(
        &self,
        _: CS,
        _: usize,
    ) -> Result<Self, SynthesisError> {
        Ok(self.clone())
    }

    fn frobenius_map_in_place<CS: ConstraintSystem<E>>(
        &mut self,
        _: CS,
        _: usize,
    ) -> Result<&mut Self, SynthesisError> {
        Ok(self)
    }

    fn mul_equals<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        other: &Self,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        cs.enforce(
            || "mul_constraint",
            |lc| &self.variable + lc,
            |lc| &other.variable + lc,
            |lc| &result.variable + lc,
        );
        Ok(())
    }

    fn square_equals<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        result: &Self,
    ) -> Result<(), SynthesisError> {
        cs.enforce(
            || "sqr_constraint",
            |lc| &self.variable + lc,
            |lc| &self.variable + lc,
            |lc| &result.variable + lc,
        );
        Ok(())
    }

    fn cost_of_mul() -> usize {
        1
    }

    fn cost_of_inv() -> usize {
        1
    }
}

impl<E: PairingEngine> PartialEq for FpGadget<E> {
    fn eq(&self, other: &Self) -> bool {
        !self.value.is_none() && !other.value.is_none() && self.value == other.value
    }
}

impl<E: PairingEngine> Eq for FpGadget<E> {}

impl<E: PairingEngine> EqGadget<E> for FpGadget<E> {}

impl<E: PairingEngine> ConditionalEqGadget<E> for FpGadget<E> {
    #[inline]
    fn conditional_enforce_equal<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        other: &Self,
        condition: &Boolean,
    ) -> Result<(), SynthesisError> {
        let difference = self.sub(cs.ns(|| "difference"), other)?;
        let one = CS::one();
        let one_const = E::Fr::one();
        cs.enforce(
            || "conditional_equals",
            |lc| &difference.variable + lc,
            |lc| lc + &condition.lc(one, one_const),
            |lc| lc,
        );
        Ok(())
    }

    fn cost() -> usize {
        1
    }
}

impl<E: PairingEngine> NEqGadget<E> for FpGadget<E> {
    #[inline]
    fn enforce_not_equal<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<(), SynthesisError> {
        let a_minus_b = self.sub(cs.ns(|| "A - B"), other)?;
        a_minus_b.inverse(cs.ns(|| "Enforce inverse exists"))?;
        Ok(())
    }

    fn cost() -> usize {
        1
    }
}

impl<E: PairingEngine> ToBitsGadget<E> for FpGadget<E> {
    /// Outputs the binary representation of the value in `self` in *big-endian*
    /// form.
    fn to_bits<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Vec<Boolean>, SynthesisError> {
        let num_bits = <E::Fr as PrimeField>::Params::MODULUS_BITS;
        use algebra::BitIterator;
        let bit_values = match self.value {
            Some(value) => {
                let mut field_char = BitIterator::new(E::Fr::characteristic());
                let mut tmp = Vec::with_capacity(num_bits as usize);
                let mut found_one = false;
                for b in BitIterator::new(value.into_repr()) {
                    // Skip leading bits
                    found_one |= field_char.next().unwrap();
                    if !found_one {
                        continue;
                    }

                    tmp.push(Some(b));
                }

                assert_eq!(tmp.len(), num_bits as usize);

                tmp
            },
            None => vec![None; num_bits as usize],
        };

        let mut bits = vec![];
        for (i, b) in bit_values.into_iter().enumerate() {
            bits.push(AllocatedBit::alloc(cs.ns(|| format!("bit {}", i)), || {
                b.get()
            })?);
        }

        let mut lc = LinearCombination::zero();
        let mut coeff = E::Fr::one();

        for bit in bits.iter().rev() {
            lc = lc + (coeff, bit.get_variable());

            coeff.double_in_place();
        }

        lc = &self.variable - lc;

        cs.enforce(|| "unpacking_constraint", |lc| lc, |lc| lc, |_| lc);

        Ok(bits.into_iter().map(Boolean::from).collect())
    }

    fn to_bits_strict<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Boolean>, SynthesisError> {
        let bits = self.to_bits(&mut cs)?;
        Boolean::enforce_in_field::<_, _, E::Fr>(&mut cs, &bits)?;

        Ok(bits)
    }
}

impl<E: PairingEngine> ToBytesGadget<E> for FpGadget<E> {
    fn to_bytes<CS: ConstraintSystem<E>>(&self, mut cs: CS) -> Result<Vec<UInt8>, SynthesisError> {
        let byte_values = match self.value {
            Some(value) => to_bytes![&value.into_repr()]?
                .into_iter()
                .map(Some)
                .collect::<Vec<_>>(),
            None => {
                let default = E::Fr::default();
                let default_len = to_bytes![&default].unwrap().len();
                vec![None; default_len]
            },
        };

        let bytes = UInt8::alloc_vec(cs.ns(|| "Alloc bytes"), &byte_values)?;

        let mut lc = LinearCombination::zero();
        let mut coeff = E::Fr::one();

        for bit in bytes
            .iter()
            .flat_map(|byte_gadget| byte_gadget.bits.clone())
        {
            match bit {
                Boolean::Is(bit) => {
                    lc = lc + (coeff, bit.get_variable());
                    coeff.double_in_place();
                },
                Boolean::Constant(_) | Boolean::Not(_) => unreachable!(),
            }
        }

        lc = &self.variable - lc;

        cs.enforce(|| "unpacking_constraint", |lc| lc, |lc| lc, |_| lc);

        Ok(bytes)
    }

    fn to_bytes_strict<CS: ConstraintSystem<E>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let bytes = self.to_bytes(&mut cs)?;
        Boolean::enforce_in_field::<_, _, E::Fr>(
            &mut cs,
            &bytes.iter()
                .flat_map(|byte_gadget| byte_gadget.into_bits_le())
                // This reverse maps the bits into big-endian form, as required by `enforce_in_field`.
                .rev()
                .collect::<Vec<_>>(),
        )?;

        Ok(bytes)
    }
}

impl<E: PairingEngine> CondSelectGadget<E> for FpGadget<E> {
    #[inline]
    fn conditionally_select<CS: ConstraintSystem<E>>(
        mut cs: CS,
        cond: &Boolean,
        first: &Self,
        second: &Self,
    ) -> Result<Self, SynthesisError> {
        if let Boolean::Constant(cond) = *cond {
            if cond {
                Ok(first.clone())
            } else {
                Ok(second.clone())
            }
        } else {
            let result = Self::alloc(cs.ns(|| ""), || {
                cond.get_value()
                    .and_then(|cond| if cond { first } else { second }.get_value())
                    .get()
            })?;
            // a = self; b = other; c = cond;
            //
            // r = c * a + (1  - c) * b
            // r = b + c * (a - b)
            // c * (a - b) = r - b
            let one = CS::one();
            cs.enforce(
                || "conditionally_select",
                |_| cond.lc(one, E::Fr::one()),
                |lc| (&first.variable - &second.variable) + lc,
                |lc| (&result.variable - &second.variable) + lc,
            );

            Ok(result)
        }
    }

    fn cost() -> usize {
        1
    }
}
/// Uses two bits to perform a lookup into a table
/// `b` is little-endian: `b[0]` is LSB.
impl<E: PairingEngine> TwoBitLookupGadget<E> for FpGadget<E> {
    type TableConstant = E::Fr;
    fn two_bit_lookup<CS: ConstraintSystem<E>>(
        mut cs: CS,
        b: &[Boolean],
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        debug_assert!(b.len() == 2);
        debug_assert!(c.len() == 4);

        let result = Self::alloc(cs.ns(|| "Allocate lookup result"), || {
            match (b[0].get_value().get()?, b[1].get_value().get()?) {
                (false, false) => Ok(c[0]),
                (false, true) => Ok(c[2]),
                (true, false) => Ok(c[1]),
                (true, true) => Ok(c[3]),
            }
        })?;
        let one = CS::one();
        cs.enforce(
            || "Enforce lookup",
            |lc| lc + b[1].lc(one, c[3] - &c[2] - &c[1] + &c[0]) + (c[1] - &c[0], one),
            |lc| lc + b[0].lc(one, E::Fr::one()),
            |lc| result.get_variable() + lc + (-c[0], one) + b[1].lc(one, c[0] - &c[2]),
        );

        Ok(result)
    }

    fn cost() -> usize {
        1
    }
}

impl<E: PairingEngine> Clone for FpGadget<E> {
    fn clone(&self) -> Self {
        Self {
            value:    self.value.clone(),
            variable: self.variable.clone(),
        }
    }
}

impl<E: PairingEngine> AllocGadget<E::Fr, E> for FpGadget<E> {
    #[inline]
    fn alloc<F, T, CS: ConstraintSystem<E>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<E::Fr>,
    {
        let mut value = None;
        let variable = cs.alloc(
            || "alloc",
            || {
                let tmp = *value_gen()?.borrow();
                value = Some(tmp);
                Ok(tmp)
            },
        )?;
        Ok(FpGadget {
            value,
            variable: Var(variable),
        })
    }

    #[inline]
    fn alloc_input<F, T, CS: ConstraintSystem<E>>(
        mut cs: CS,
        value_gen: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<E::Fr>,
    {
        let mut value = None;
        let variable = cs.alloc_input(
            || "alloc",
            || {
                let tmp = *value_gen()?.borrow();
                value = Some(tmp);
                Ok(tmp)
            },
        )?;
        Ok(FpGadget {
            value,
            variable: Var(variable),
        })
    }
}
