pub mod nonnative_group_gadget;

use crate::prelude::*;
use algebra::{Field, Group};
use r1cs_core::{ConstraintSystem, SynthesisError};

use std::{borrow::Borrow, fmt::Debug};

pub trait NonNativeGroupGadget<G: Group, ConstraintF: Field, SimulationF: Field>:
    Sized
    + ToBytesGadget<ConstraintF>
    + EqGadget<ConstraintF>
    + ToBitsGadget<ConstraintF>
    + CondSelectGadget<ConstraintF>
    + AllocGadget<G, ConstraintF>
    + ConstantGadget<G, ConstraintF>
    + Clone
    + Debug
{
    type Value: Debug;
    type Variable;

    fn add<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError>;

    fn zero<CS: ConstraintSystem<ConstraintF>>(cs: CS) -> Result<Self, SynthesisError>;

    fn is_zero<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Boolean, SynthesisError>;

    fn double_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        cs: CS,
    ) -> Result<(), SynthesisError>;

    fn negate<CS: ConstraintSystem<ConstraintF>>(&self, cs: CS) -> Result<Self, SynthesisError>;

    fn add_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        other: &G,
    ) -> Result<Self, SynthesisError>;

    fn sub_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &G,
    ) -> Result<Self, SynthesisError> {
        let neg_other = -(*other);
        self.add_constant(cs.ns(|| "Self - other"), &neg_other)
    }

    /// Variable base exponentiation.
    /// Inputs must be specified in *little-endian* form.
    /// If the addition law is incomplete for the identity element,
    /// `result` must not be the identity element.
    fn mul_bits<'a, CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        result: &Self,
        bits: impl Iterator<Item = &'a Boolean>,
    ) -> Result<Self, SynthesisError> {
        let mut power = self.clone();
        let mut result = result.clone();
        for (i, bit) in bits.enumerate() {
            let new_encoded = result.add(&mut cs.ns(|| format!("Add {}-th power", i)), &power)?;
            result = Self::conditionally_select(
                &mut cs.ns(|| format!("Select {}", i)),
                bit.borrow(),
                &new_encoded,
                &result,
            )?;
            power.double_in_place(&mut cs.ns(|| format!("{}-th Doubling", i)))?;
        }
        Ok(result)
    }


    /// Fixed base exponentiation, slighlty different interface from
    /// `precomputed_base_scalar_mul`. Inputs must be specified in
    /// *little-endian* form. If the addition law is incomplete for
    /// the identity element, `result` must not be the identity element.
    fn mul_bits_fixed_base<'a, CS: ConstraintSystem<ConstraintF>>(
        base: &'a G,
        mut cs: CS,
        result: &Self,
        bits: &[Boolean],
    ) -> Result<Self, SynthesisError> {
        let base_g = Self::from_value(cs.ns(|| "hardcode base"), base);
        base_g.mul_bits(cs, result, bits.into_iter())
    }

}


