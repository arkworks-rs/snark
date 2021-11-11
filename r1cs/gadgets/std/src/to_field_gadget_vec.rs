use algebra::{curves::models::SWModelParameters, curves::models::TEModelParameters, PrimeField};

use crate::fields::fp::FpGadget;
use crate::{
    fields::FieldGadget,
    groups::curves::short_weierstrass::{
        short_weierstrass_jacobian::AffineGadget as SWJAffineGadget,
        short_weierstrass_projective::AffineGadget as SWPAffineGadget,
    },
    groups::curves::twisted_edwards::AffineGadget as TEAffineGadget,
};
use r1cs_core::{ConstraintSystem, SynthesisError as Error};

/// Types that can be converted to a vector of elements that implement the `Field Gadget` trait.
pub trait ToConstraintFieldGadget<ConstraintF: PrimeField> {
    type FieldGadget: FieldGadget<ConstraintF, ConstraintF>;

    fn to_field_gadget_elements<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
    ) -> Result<Vec<Self::FieldGadget>, Error>;
}

impl<ConstraintF: PrimeField> ToConstraintFieldGadget<ConstraintF> for FpGadget<ConstraintF> {
    type FieldGadget = Self;

    fn to_field_gadget_elements<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Self::FieldGadget>, Error> {
        Ok(vec![self.clone()])
    }
}

impl<ConstraintF: PrimeField> ToConstraintFieldGadget<ConstraintF> for [FpGadget<ConstraintF>] {
    type FieldGadget = FpGadget<ConstraintF>;

    #[inline]
    fn to_field_gadget_elements<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Self::FieldGadget>, Error> {
        Ok(self.to_vec())
    }
}

impl<ConstraintF: PrimeField> ToConstraintFieldGadget<ConstraintF> for () {
    type FieldGadget = FpGadget<ConstraintF>;

    #[inline]
    fn to_field_gadget_elements<CS: ConstraintSystem<ConstraintF>>(
        &self,
        _cs: CS,
    ) -> Result<Vec<Self::FieldGadget>, Error> {
        Ok(Vec::new())
    }
}

impl<M, ConstraintF, FG> ToConstraintFieldGadget<ConstraintF>
    for SWPAffineGadget<M, ConstraintF, FG>
where
    M: SWModelParameters,
    ConstraintF: PrimeField,
    FG: FieldGadget<M::BaseField, ConstraintF>
        + ToConstraintFieldGadget<ConstraintF, FieldGadget = FpGadget<ConstraintF>>,
{
    type FieldGadget = FpGadget<ConstraintF>;

    #[inline]
    fn to_field_gadget_elements<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Self::FieldGadget>, Error> {
        let mut x_fe = self.x.to_field_gadget_elements(cs.ns(|| "x"))?;
        let y_fe = self.y.to_field_gadget_elements(cs.ns(|| "y"))?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}

impl<M, ConstraintF, FG> ToConstraintFieldGadget<ConstraintF>
    for SWJAffineGadget<M, ConstraintF, FG>
where
    M: SWModelParameters,
    ConstraintF: PrimeField,
    FG: FieldGadget<M::BaseField, ConstraintF>
        + ToConstraintFieldGadget<ConstraintF, FieldGadget = FpGadget<ConstraintF>>,
{
    type FieldGadget = FpGadget<ConstraintF>;

    #[inline]
    fn to_field_gadget_elements<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Self::FieldGadget>, Error> {
        let mut x_fe = self.x.to_field_gadget_elements(cs.ns(|| "x"))?;
        let y_fe = self.y.to_field_gadget_elements(cs.ns(|| "y"))?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}

impl<M, ConstraintF, FG> ToConstraintFieldGadget<ConstraintF> for TEAffineGadget<M, ConstraintF, FG>
where
    M: TEModelParameters,
    ConstraintF: PrimeField,
    FG: FieldGadget<M::BaseField, ConstraintF>
        + ToConstraintFieldGadget<ConstraintF, FieldGadget = FpGadget<ConstraintF>>,
{
    type FieldGadget = FpGadget<ConstraintF>;

    #[inline]
    fn to_field_gadget_elements<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<Self::FieldGadget>, Error> {
        let mut x_fe = self.x.to_field_gadget_elements(cs.ns(|| "x"))?;
        let y_fe = self.y.to_field_gadget_elements(cs.ns(|| "y"))?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}
