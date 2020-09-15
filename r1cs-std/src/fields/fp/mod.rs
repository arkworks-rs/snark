use algebra::{BigInteger, FpParameters, PrimeField};
use r1cs_core::{lc, ConstraintSystemRef, LinearCombination, Namespace, SynthesisError, Variable};

use core::borrow::Borrow;

use crate::fields::{FieldOpsBounds, FieldVar};
use crate::{prelude::*, Assignment, ToConstraintFieldGadget, Vec};

pub mod cmp;

#[derive(Debug, Clone)]
#[must_use]
pub struct AllocatedFp<F: PrimeField> {
    pub(crate) value: Option<F>,
    pub variable: Variable,
    pub cs: ConstraintSystemRef<F>,
}

impl<F: PrimeField> AllocatedFp<F> {
    pub fn new(value: Option<F>, variable: Variable, cs: ConstraintSystemRef<F>) -> Self {
        Self {
            value,
            variable,
            cs,
        }
    }
}

/// Represent variables corresponding to the field `F`.
#[derive(Clone, Debug)]
#[must_use]
pub enum FpVar<F: PrimeField> {
    Constant(F),
    Var(AllocatedFp<F>),
}

impl<F: PrimeField> R1CSVar<F> for FpVar<F> {
    type Value = F;

    fn cs(&self) -> Option<ConstraintSystemRef<F>> {
        match self {
            Self::Constant(_) => None,
            Self::Var(a) => Some(a.cs.clone()),
        }
    }

    fn value(&self) -> Result<Self::Value, SynthesisError> {
        match self {
            Self::Constant(v) => Ok(*v),
            Self::Var(v) => v.value(),
        }
    }
}

impl<F: PrimeField> From<Boolean<F>> for FpVar<F> {
    fn from(other: Boolean<F>) -> Self {
        if let Boolean::Constant(b) = other {
            Self::Constant(F::from(b as u8))
        } else {
            // `other` is a variable
            let cs = other.cs().unwrap();
            let variable = cs.new_lc(other.lc()).unwrap();
            Self::Var(AllocatedFp::new(
                other.value().ok().map(|b| F::from(b as u8)),
                variable,
                cs,
            ))
        }
    }
}

impl<F: PrimeField> From<AllocatedFp<F>> for FpVar<F> {
    fn from(other: AllocatedFp<F>) -> Self {
        Self::Var(other)
    }
}

impl<'a, F: PrimeField> FieldOpsBounds<'a, F, Self> for FpVar<F> {}
impl<'a, F: PrimeField> FieldOpsBounds<'a, F, FpVar<F>> for &'a FpVar<F> {}

impl<F: PrimeField> AllocatedFp<F> {
    pub fn from(other: Boolean<F>) -> Self {
        if let Some(cs) = other.cs() {
            let variable = cs.new_lc(other.lc()).unwrap();
            Self::new(other.value().ok().map(|b| F::from(b as u8)), variable, cs)
        } else {
            unreachable!("Cannot create a constant value")
        }
    }

    pub fn value(&self) -> Result<F, SynthesisError> {
        self.cs.assigned_value(self.variable).get()
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn add(&self, other: &Self) -> Self {
        let value = match (self.value, other.value) {
            (Some(val1), Some(val2)) => Some(val1 + &val2),
            (..) => None,
        };

        let variable = self
            .cs
            .new_lc(lc!() + self.variable + other.variable)
            .unwrap();
        AllocatedFp::new(value, variable, self.cs.clone())
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn sub(&self, other: &Self) -> Self {
        let value = match (self.value, other.value) {
            (Some(val1), Some(val2)) => Some(val1 - &val2),
            (..) => None,
        };

        let variable = self
            .cs
            .new_lc(lc!() + self.variable - other.variable)
            .unwrap();
        AllocatedFp::new(value, variable, self.cs.clone())
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn mul(&self, other: &Self) -> Self {
        let product = AllocatedFp::new_witness(self.cs.clone(), || {
            Ok(self.value.get()? * &other.value.get()?)
        })
        .unwrap();
        self.cs
            .enforce_constraint(
                lc!() + self.variable,
                lc!() + other.variable,
                lc!() + product.variable,
            )
            .unwrap();
        product
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn add_constant(&self, other: F) -> Self {
        if other.is_zero() {
            self.clone()
        } else {
            let value = self.value.map(|val| val + other);
            let variable = self
                .cs
                .new_lc(lc!() + self.variable + (other, Variable::One))
                .unwrap();
            AllocatedFp::new(value, variable, self.cs.clone())
        }
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn sub_constant(&self, other: F) -> Self {
        self.add_constant(-other)
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn mul_constant(&self, other: F) -> Self {
        if other.is_one() {
            self.clone()
        } else {
            let value = self.value.map(|val| val * other);
            let variable = self.cs.new_lc(lc!() + (other, self.variable)).unwrap();
            AllocatedFp::new(value, variable, self.cs.clone())
        }
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn double(&self) -> Result<Self, SynthesisError> {
        let value = self.value.map(|val| val.double());
        let variable = self.cs.new_lc(lc!() + self.variable + self.variable)?;
        Ok(Self::new(value, variable, self.cs.clone()))
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn negate(&self) -> Self {
        let mut result = self.clone();
        result.negate_in_place();
        result
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn negate_in_place(&mut self) -> &mut Self {
        self.value.as_mut().map(|val| *val = -(*val));
        self.variable = self.cs.new_lc(lc!() - self.variable).unwrap();
        self
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn square(&self) -> Result<Self, SynthesisError> {
        Ok(self.mul(self))
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn inverse(&self) -> Result<Self, SynthesisError> {
        let inverse = Self::new_witness(self.cs.clone(), || {
            Ok(self.value.get()?.inverse().unwrap_or(F::zero()))
        })?;

        self.cs.enforce_constraint(
            lc!() + self.variable,
            lc!() + inverse.variable,
            lc!() + Variable::One,
        )?;
        Ok(inverse)
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn frobenius_map(&self, _: usize) -> Result<Self, SynthesisError> {
        Ok(self.clone())
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn mul_equals(&self, other: &Self, result: &Self) -> Result<(), SynthesisError> {
        self.cs.enforce_constraint(
            lc!() + self.variable,
            lc!() + other.variable,
            lc!() + result.variable,
        )
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn square_equals(&self, result: &Self) -> Result<(), SynthesisError> {
        self.cs.enforce_constraint(
            lc!() + self.variable,
            lc!() + self.variable,
            lc!() + result.variable,
        )
    }

    /// Outputs the bit `self == other`.
    ///
    /// # Constraint cost
    ///
    /// Consumes three constraints
    #[tracing::instrument(target = "r1cs")]
    pub fn is_eq(&self, other: &Self) -> Result<Boolean<F>, SynthesisError> {
        Ok(self.is_neq(other)?.not())
    }

    /// Outputs the bit `self != other`.
    ///
    /// # Constraint cost
    ///
    /// Consumes three constraints
    #[tracing::instrument(target = "r1cs")]
    pub fn is_neq(&self, other: &Self) -> Result<Boolean<F>, SynthesisError> {
        let is_not_equal = Boolean::new_witness(self.cs.clone(), || {
            Ok(self.value.get()? != other.value.get()?)
        })?;
        let multiplier = self.cs.new_witness_variable(|| {
            if is_not_equal.value()? {
                (self.value.get()? - other.value.get()?).inverse().get()
            } else {
                Ok(F::one())
            }
        })?;

        // Completeness:
        // Case 1: self != other:
        // ----------------------
        //   constraint 1:
        //   (self - other) * multiplier = is_not_equal
        //   => (non_zero) * multiplier = 1 (satisfied, because multiplier = 1/(self - other)
        //
        //   constraint 2:
        //   (self - other) * not(is_not_equal) = 0
        //   => (non_zero) * not(1) = 0
        //   => (non_zero) * 0 = 0
        //
        // Case 2: self == other:
        // ----------------------
        //   constraint 1:
        //   (self - other) * multiplier = is_not_equal
        //   => 0 * multiplier = 0 (satisfied, because multiplier = 1
        //
        //   constraint 2:
        //   (self - other) * not(is_not_equal) = 0
        //   => 0 * not(0) = 0
        //   => 0 * 1 = 0
        //
        // --------------------------------------------------------------------
        //
        // Soundness:
        // Case 1: self != other, but is_not_equal = 0.
        // --------------------------------------------
        //   constraint 1:
        //   (self - other) * multiplier = is_not_equal
        //   => non_zero * multiplier = 0 (only satisfiable if multiplier == 0)
        //
        //   constraint 2:
        //   (self - other) * not(is_not_equal) = 0
        //   => (non_zero) * 1 = 0 (impossible)
        //
        // Case 2: self == other, but is_not_equal = 1.
        // --------------------------------------------
        //   constraint 1:
        //   (self - other) * multiplier = is_not_equal
        //   0 * multiplier = 1 (unsatisfiable)
        self.cs.enforce_constraint(
            lc!() + self.variable - other.variable,
            lc!() + multiplier,
            is_not_equal.lc(),
        )?;
        self.cs.enforce_constraint(
            lc!() + self.variable - other.variable,
            is_not_equal.not().lc(),
            lc!(),
        )?;
        Ok(is_not_equal)
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn conditional_enforce_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        self.cs.enforce_constraint(
            lc!() + self.variable - other.variable,
            lc!() + should_enforce.lc(),
            lc!(),
        )
    }

    #[tracing::instrument(target = "r1cs")]
    pub fn conditional_enforce_not_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        let multiplier = Self::new_witness(self.cs.clone(), || {
            if should_enforce.value()? {
                (self.value.get()? - other.value.get()?).inverse().get()
            } else {
                Ok(F::zero())
            }
        })?;

        self.cs.enforce_constraint(
            lc!() + self.variable - other.variable,
            lc!() + multiplier.variable,
            should_enforce.lc(),
        )?;
        Ok(())
    }
}

/****************************************************************************/
/****************************************************************************/

impl<F: PrimeField> ToBitsGadget<F> for AllocatedFp<F> {
    /// Outputs the unique bit-wise decomposition of `self` in *little-endian*
    /// form.
    #[tracing::instrument(target = "r1cs")]
    fn to_bits_le(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        let bits = self.to_non_unique_bits_le()?;
        Boolean::enforce_in_field_le(&bits)?;
        Ok(bits)
    }

    #[tracing::instrument(target = "r1cs")]
    fn to_non_unique_bits_le(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        let cs = self.cs.clone();
        use algebra::BitIteratorBE;
        let mut bits = if let Some(value) = self.value {
            let field_char = BitIteratorBE::new(F::characteristic());
            let bits: Vec<_> = BitIteratorBE::new(value.into_repr())
                .zip(field_char)
                .skip_while(|(_, c)| !c)
                .map(|(b, _)| Some(b))
                .collect();
            assert_eq!(bits.len(), F::Params::MODULUS_BITS as usize);
            bits
        } else {
            vec![None; F::Params::MODULUS_BITS as usize]
        };

        // Convert to little-endian
        bits.reverse();

        let bits: Vec<_> = bits
            .into_iter()
            .map(|b| Boolean::new_witness(cs.clone(), || b.get()))
            .collect::<Result<_, _>>()?;

        let mut lc = LinearCombination::zero();
        let mut coeff = F::one();

        for bit in bits.iter() {
            lc = &lc + bit.lc() * coeff;

            coeff.double_in_place();
        }

        lc = lc - &self.variable;

        cs.enforce_constraint(lc!(), lc!(), lc)?;

        Ok(bits)
    }
}

impl<F: PrimeField> ToBytesGadget<F> for AllocatedFp<F> {
    /// Outputs the unique byte decomposition of `self` in *little-endian*
    /// form.
    #[tracing::instrument(target = "r1cs")]
    fn to_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        let num_bits = F::BigInt::NUM_LIMBS * 64;
        let mut bits = self.to_bits_le()?;
        let remainder = core::iter::repeat(Boolean::constant(false)).take(num_bits - bits.len());
        bits.extend(remainder);
        let bytes = bits
            .chunks(8)
            .map(|chunk| UInt8::from_bits_le(chunk))
            .collect();
        Ok(bytes)
    }

    #[tracing::instrument(target = "r1cs")]
    fn to_non_unique_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        let num_bits = F::BigInt::NUM_LIMBS * 64;
        let mut bits = self.to_non_unique_bits_le()?;
        let remainder = core::iter::repeat(Boolean::constant(false)).take(num_bits - bits.len());
        bits.extend(remainder);
        let bytes = bits
            .chunks(8)
            .map(|chunk| UInt8::from_bits_le(chunk))
            .collect();
        Ok(bytes)
    }
}

impl<F: PrimeField> ToConstraintFieldGadget<F> for AllocatedFp<F> {
    #[tracing::instrument(target = "r1cs")]
    fn to_constraint_field(&self) -> Result<Vec<FpVar<F>>, SynthesisError> {
        Ok(vec![self.clone().into()])
    }
}

impl<F: PrimeField> CondSelectGadget<F> for AllocatedFp<F> {
    #[inline]
    #[tracing::instrument(target = "r1cs")]
    fn conditionally_select(
        cond: &Boolean<F>,
        true_val: &Self,
        false_val: &Self,
    ) -> Result<Self, SynthesisError> {
        match cond {
            Boolean::Constant(true) => Ok(true_val.clone()),
            Boolean::Constant(false) => Ok(false_val.clone()),
            _ => {
                let cs = cond.cs().unwrap();
                let result = Self::new_witness(cs.clone(), || {
                    cond.value()
                        .and_then(|c| if c { true_val } else { false_val }.value.get())
                })?;
                // a = self; b = other; c = cond;
                //
                // r = c * a + (1  - c) * b
                // r = b + c * (a - b)
                // c * (a - b) = r - b
                cs.enforce_constraint(
                    cond.lc(),
                    lc!() + true_val.variable - false_val.variable,
                    lc!() + result.variable - false_val.variable,
                )?;

                Ok(result)
            }
        }
    }
}

/// Uses two bits to perform a lookup into a table
/// `b` is little-endian: `b[0]` is LSB.
impl<F: PrimeField> TwoBitLookupGadget<F> for AllocatedFp<F> {
    type TableConstant = F;
    #[tracing::instrument(target = "r1cs")]
    fn two_bit_lookup(b: &[Boolean<F>], c: &[Self::TableConstant]) -> Result<Self, SynthesisError> {
        debug_assert_eq!(b.len(), 2);
        debug_assert_eq!(c.len(), 4);
        if let Some(cs) = b.cs() {
            let result = Self::new_witness(cs.clone(), || {
                let lsb = usize::from(b[0].value()?);
                let msb = usize::from(b[1].value()?);
                let index = lsb + (msb << 1);
                Ok(c[index])
            })?;
            let one = Variable::One;
            cs.enforce_constraint(
                lc!() + b[1].lc() * (c[3] - &c[2] - &c[1] + &c[0]) + (c[1] - &c[0], one),
                lc!() + b[0].lc(),
                lc!() + result.variable - (c[0], one) + b[1].lc() * (c[0] - &c[2]),
            )?;

            Ok(result)
        } else {
            unreachable!("must provide a way to obtain a ConstraintSystemRef")
        }
    }
}

impl<F: PrimeField> ThreeBitCondNegLookupGadget<F> for AllocatedFp<F> {
    type TableConstant = F;

    #[tracing::instrument(target = "r1cs")]
    fn three_bit_cond_neg_lookup(
        b: &[Boolean<F>],
        b0b1: &Boolean<F>,
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        debug_assert_eq!(b.len(), 3);
        debug_assert_eq!(c.len(), 4);

        if let Some(cs) = b.cs() {
            let result = Self::new_witness(cs.clone(), || {
                let lsb = usize::from(b[0].value()?);
                let msb = usize::from(b[1].value()?);
                let index = lsb + (msb << 1);
                let intermediate = c[index];

                let is_negative = b[2].value()?;
                let y = if is_negative {
                    -intermediate
                } else {
                    intermediate
                };
                Ok(y)
            })?;

            let y_lc = b0b1.lc() * (c[3] - &c[2] - &c[1] + &c[0])
                + b[0].lc() * (c[1] - &c[0])
                + b[1].lc() * (c[2] - &c[0])
                + (c[0], Variable::One);
            cs.enforce_constraint(
                y_lc.clone() + y_lc.clone(),
                b[2].lc(),
                y_lc.clone() - result.variable,
            )?;

            Ok(result)
        } else {
            unreachable!("must provide a way to obtain a ConstraintSystemRef")
        }
    }
}

impl<F: PrimeField> AllocVar<F, F> for AllocatedFp<F> {
    fn new_variable<T: Borrow<F>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();
        if mode == AllocationMode::Constant {
            let v = *f()?.borrow();
            let lc = cs.new_lc(lc!() + (v, Variable::One))?;
            Ok(Self::new(Some(v), lc, cs))
        } else {
            let mut value = None;
            let value_generator = || {
                value = Some(*f()?.borrow());
                value.ok_or(SynthesisError::AssignmentMissing)
            };
            let variable = if mode == AllocationMode::Input {
                cs.new_input_variable(value_generator)?
            } else {
                cs.new_witness_variable(value_generator)?
            };
            Ok(Self::new(value, variable, cs.clone()))
        }
    }
}

impl<F: PrimeField> FieldVar<F, F> for FpVar<F> {
    fn constant(f: F) -> Self {
        Self::Constant(f)
    }

    fn zero() -> Self {
        Self::Constant(F::zero())
    }

    fn one() -> Self {
        Self::Constant(F::one())
    }

    #[tracing::instrument(target = "r1cs")]
    fn double(&self) -> Result<Self, SynthesisError> {
        match self {
            Self::Constant(c) => Ok(Self::Constant(c.double())),
            Self::Var(v) => Ok(Self::Var(v.double()?)),
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn negate(&self) -> Result<Self, SynthesisError> {
        match self {
            Self::Constant(c) => Ok(Self::Constant(-*c)),
            Self::Var(v) => Ok(Self::Var(v.negate())),
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn square(&self) -> Result<Self, SynthesisError> {
        match self {
            Self::Constant(c) => Ok(Self::Constant(c.square())),
            Self::Var(v) => Ok(Self::Var(v.square()?)),
        }
    }

    /// Enforce that `self * other == result`.
    #[tracing::instrument(target = "r1cs")]
    fn mul_equals(&self, other: &Self, result: &Self) -> Result<(), SynthesisError> {
        use FpVar::*;
        match (self, other, result) {
            (Constant(_), Constant(_), Constant(_)) => Ok(()),
            (Constant(_), Constant(_), _) | (Constant(_), Var(_), _) | (Var(_), Constant(_), _) => {
                result.enforce_equal(&(self * other))
            } // this multiplication should be free
            (Var(v1), Var(v2), Var(v3)) => v1.mul_equals(v2, v3),
            (Var(v1), Var(v2), Constant(f)) => {
                let cs = v1.cs.clone();
                let v3 = AllocatedFp::new_constant(cs.clone(), f).unwrap();
                v1.mul_equals(v2, &v3)
            }
        }
    }

    /// Enforce that `self * self == result`.
    #[tracing::instrument(target = "r1cs")]
    fn square_equals(&self, result: &Self) -> Result<(), SynthesisError> {
        use FpVar::*;
        match (self, result) {
            (Constant(_), Constant(_)) => Ok(()),
            (Constant(f), Var(r)) => {
                let cs = r.cs.clone();
                let v = AllocatedFp::new_witness(cs, || Ok(f))?;
                v.square_equals(&r)
            }
            (Var(v), Constant(f)) => {
                let cs = v.cs.clone();
                let r = AllocatedFp::new_witness(cs, || Ok(f))?;
                v.square_equals(&r)
            }
            (Var(v1), Var(v2)) => v1.square_equals(v2),
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn inverse(&self) -> Result<Self, SynthesisError> {
        match self {
            FpVar::Var(v) => v.inverse().map(FpVar::Var),
            FpVar::Constant(f) => f.inverse().get().map(FpVar::Constant),
        }
    }

    /// Returns (self / denominator), but requires fewer constraints than
    /// self * denominator.inverse()
    /// It is up to the caller to ensure that denominator is non-zero,
    /// since in that case the result is unconstrained.
    #[tracing::instrument(target = "r1cs")]
    fn mul_by_inverse(&self, denominator: &Self) -> Result<Self, SynthesisError> {
        use FpVar::*;
        match (self, denominator) {
            (Constant(s), Constant(d)) => Ok(Constant(*s / *d)),
            (Var(s), Constant(d)) => Ok(Var(s.mul_constant(d.inverse().get()?))),
            (Constant(s), Var(d)) => Ok(Var(d.inverse()?.mul_constant(*s))),
            (Var(s), Var(d)) => Ok(Var(d.inverse()?.mul(s))),
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn frobenius_map(&self, power: usize) -> Result<Self, SynthesisError> {
        match self {
            FpVar::Var(v) => v.frobenius_map(power).map(FpVar::Var),
            FpVar::Constant(f) => {
                let mut f = *f;
                f.frobenius_map(power);
                Ok(FpVar::Constant(f))
            }
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn frobenius_map_in_place(&mut self, power: usize) -> Result<&mut Self, SynthesisError> {
        *self = self.frobenius_map(power)?;
        Ok(self)
    }
}

/****************************************************************************/
/****************************************************************************/

impl_ops!(
    FpVar<F>,
    F,
    Add,
    add,
    AddAssign,
    add_assign,
    |this: &'a FpVar<F>, other: &'a FpVar<F>| {
        use FpVar::*;
        match (this, other) {
            (Constant(c1), Constant(c2)) => Constant(*c1 + *c2),
            (Constant(c), Var(v)) | (Var(v), Constant(c)) => Var(v.add_constant(*c)),
            (Var(v1), Var(v2)) => Var(v1.add(v2)),
        }
    },
    |this: &'a FpVar<F>, other: F| { this + &FpVar::Constant(other) },
    F: PrimeField,
);

impl_ops!(
    FpVar<F>,
    F,
    Sub,
    sub,
    SubAssign,
    sub_assign,
    |this: &'a FpVar<F>, other: &'a FpVar<F>| {
        use FpVar::*;
        match (this, other) {
            (Constant(c1), Constant(c2)) => Constant(*c1 - *c2),
            (Var(v), Constant(c)) => Var(v.sub_constant(*c)),
            (Constant(c), Var(v)) => Var(v.sub_constant(*c).negate()),
            (Var(v1), Var(v2)) => Var(v1.sub(v2)),
        }
    },
    |this: &'a FpVar<F>, other: F| { this - &FpVar::Constant(other) },
    F: PrimeField
);

impl_ops!(
    FpVar<F>,
    F,
    Mul,
    mul,
    MulAssign,
    mul_assign,
    |this: &'a FpVar<F>, other: &'a FpVar<F>| {
        use FpVar::*;
        match (this, other) {
            (Constant(c1), Constant(c2)) => Constant(*c1 * *c2),
            (Constant(c), Var(v)) | (Var(v), Constant(c)) => Var(v.mul_constant(*c)),
            (Var(v1), Var(v2)) => Var(v1.mul(v2)),
        }
    },
    |this: &'a FpVar<F>, other: F| {
        if other.is_zero() {
            FpVar::zero()
        } else {
            this * &FpVar::Constant(other)
        }
    },
    F: PrimeField
);

/****************************************************************************/
/****************************************************************************/

impl<F: PrimeField> EqGadget<F> for FpVar<F> {
    #[tracing::instrument(target = "r1cs")]
    fn is_eq(&self, other: &Self) -> Result<Boolean<F>, SynthesisError> {
        match (self, other) {
            (Self::Constant(c1), Self::Constant(c2)) => Ok(Boolean::Constant(c1 == c2)),
            (Self::Constant(c), Self::Var(v)) | (Self::Var(v), Self::Constant(c)) => {
                let cs = v.cs.clone();
                let c = AllocatedFp::new_constant(cs, c)?;
                c.is_eq(v)
            }
            (Self::Var(v1), Self::Var(v2)) => v1.is_eq(v2),
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn conditional_enforce_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        match (self, other) {
            (Self::Constant(_), Self::Constant(_)) => Ok(()),
            (Self::Constant(c), Self::Var(v)) | (Self::Var(v), Self::Constant(c)) => {
                let cs = v.cs.clone();
                let c = AllocatedFp::new_constant(cs, c)?;
                c.conditional_enforce_equal(v, should_enforce)
            }
            (Self::Var(v1), Self::Var(v2)) => v1.conditional_enforce_equal(v2, should_enforce),
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn conditional_enforce_not_equal(
        &self,
        other: &Self,
        should_enforce: &Boolean<F>,
    ) -> Result<(), SynthesisError> {
        match (self, other) {
            (Self::Constant(_), Self::Constant(_)) => Ok(()),
            (Self::Constant(c), Self::Var(v)) | (Self::Var(v), Self::Constant(c)) => {
                let cs = v.cs.clone();
                let c = AllocatedFp::new_constant(cs, c)?;
                c.conditional_enforce_not_equal(v, should_enforce)
            }
            (Self::Var(v1), Self::Var(v2)) => v1.conditional_enforce_not_equal(v2, should_enforce),
        }
    }
}

impl<F: PrimeField> ToBitsGadget<F> for FpVar<F> {
    #[tracing::instrument(target = "r1cs")]
    fn to_bits_le(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        match self {
            Self::Constant(_) => self.to_non_unique_bits_le(),
            Self::Var(v) => v.to_bits_le(),
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn to_non_unique_bits_le(&self) -> Result<Vec<Boolean<F>>, SynthesisError> {
        use algebra::BitIteratorLE;
        match self {
            Self::Constant(c) => Ok(BitIteratorLE::without_trailing_zeros(&c.into_repr())
                .map(Boolean::constant)
                .collect::<Vec<_>>()),
            Self::Var(v) => v.to_non_unique_bits_le(),
        }
    }
}

impl<F: PrimeField> ToBytesGadget<F> for FpVar<F> {
    /// Outputs the unique byte decomposition of `self` in *little-endian*
    /// form.
    #[tracing::instrument(target = "r1cs")]
    fn to_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        match self {
            Self::Constant(c) => Ok(UInt8::constant_vec(&to_bytes![c].unwrap())),
            Self::Var(v) => v.to_bytes(),
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn to_non_unique_bytes(&self) -> Result<Vec<UInt8<F>>, SynthesisError> {
        match self {
            Self::Constant(c) => Ok(UInt8::constant_vec(&to_bytes![c].unwrap())),
            Self::Var(v) => v.to_non_unique_bytes(),
        }
    }
}

impl<F: PrimeField> ToConstraintFieldGadget<F> for FpVar<F> {
    #[tracing::instrument(target = "r1cs")]
    fn to_constraint_field(&self) -> Result<Vec<FpVar<F>>, SynthesisError> {
        Ok(vec![self.clone()])
    }
}

impl<F: PrimeField> CondSelectGadget<F> for FpVar<F> {
    #[tracing::instrument(target = "r1cs")]
    fn conditionally_select(
        cond: &Boolean<F>,
        true_value: &Self,
        false_value: &Self,
    ) -> Result<Self, SynthesisError> {
        match cond {
            Boolean::Constant(true) => Ok(true_value.clone()),
            Boolean::Constant(false) => Ok(false_value.clone()),
            _ => {
                match (true_value, false_value) {
                    (Self::Constant(t), Self::Constant(f)) => {
                        let is = AllocatedFp::from(cond.clone());
                        let not = AllocatedFp::from(cond.not());
                        // cond * t + (1 - cond) * f
                        Ok(is.mul_constant(*t).add(&not.mul_constant(*f)).into())
                    }
                    (_, _) => {
                        let cs = cond.cs().unwrap();
                        let true_value = match true_value {
                            Self::Constant(f) => AllocatedFp::new_constant(cs.clone(), f)?,
                            Self::Var(v) => v.clone(),
                        };
                        let false_value = match false_value {
                            Self::Constant(f) => AllocatedFp::new_constant(cs.clone(), f)?,
                            Self::Var(v) => v.clone(),
                        };
                        cond.select(&true_value, &false_value).map(Self::Var)
                    }
                }
            }
        }
    }
}

/// Uses two bits to perform a lookup into a table
/// `b` is little-endian: `b[0]` is LSB.
impl<F: PrimeField> TwoBitLookupGadget<F> for FpVar<F> {
    type TableConstant = F;

    #[tracing::instrument(target = "r1cs")]
    fn two_bit_lookup(b: &[Boolean<F>], c: &[Self::TableConstant]) -> Result<Self, SynthesisError> {
        debug_assert_eq!(b.len(), 2);
        debug_assert_eq!(c.len(), 4);
        if b.cs().is_some() {
            AllocatedFp::two_bit_lookup(b, c).map(Self::Var)
        } else {
            let lsb = usize::from(b[0].value()?);
            let msb = usize::from(b[1].value()?);
            let index = lsb + (msb << 1);
            Ok(Self::Constant(c[index]))
        }
    }
}

impl<F: PrimeField> ThreeBitCondNegLookupGadget<F> for FpVar<F> {
    type TableConstant = F;

    #[tracing::instrument(target = "r1cs")]
    fn three_bit_cond_neg_lookup(
        b: &[Boolean<F>],
        b0b1: &Boolean<F>,
        c: &[Self::TableConstant],
    ) -> Result<Self, SynthesisError> {
        debug_assert_eq!(b.len(), 3);
        debug_assert_eq!(c.len(), 4);

        if b.cs().or(b0b1.cs()).is_some() {
            AllocatedFp::three_bit_cond_neg_lookup(b, b0b1, c).map(Self::Var)
        } else {
            let lsb = usize::from(b[0].value()?);
            let msb = usize::from(b[1].value()?);
            let index = lsb + (msb << 1);
            let intermediate = c[index];

            let is_negative = b[2].value()?;
            let y = if is_negative {
                -intermediate
            } else {
                intermediate
            };
            Ok(Self::Constant(y))
        }
    }
}

impl<F: PrimeField> AllocVar<F, F> for FpVar<F> {
    fn new_variable<T: Borrow<F>>(
        cs: impl Into<Namespace<F>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        if mode == AllocationMode::Constant {
            Ok(Self::Constant(*f()?.borrow()))
        } else {
            AllocatedFp::new_variable(cs, f, mode).map(Self::Var)
        }
    }
}
