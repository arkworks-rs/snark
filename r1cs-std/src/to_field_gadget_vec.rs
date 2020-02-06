use algebra::{
    curves::models::SWModelParameters,
    PrimeField,
};

use crate::{
    fields::FieldGadget,
    groups::curves::short_weierstrass::short_weierstrass_projective::AffineGadget as SWPAffineGadget,
};
use crate::fields::fp::FpGadget;

type Error = r1cs_core::SynthesisError;

/// Types that can be converted to a vector of elements that implement the `Field Gadget` trait.
pub trait ToConstraintFieldGadget<ConstraintF: PrimeField> {

    type FieldGadget: FieldGadget<ConstraintF, ConstraintF>;

    fn to_field_gadget_elements(&self) -> Result<Vec<Self::FieldGadget>, Error>;
}

impl<ConstraintF: PrimeField> ToConstraintFieldGadget<ConstraintF> for FpGadget<ConstraintF> {
    type FieldGadget = Self;

    fn to_field_gadget_elements(&self) -> Result<Vec<Self::FieldGadget>, Error> {
        Ok(vec![self.clone()])
    }
}

impl<ConstraintF: PrimeField> ToConstraintFieldGadget<ConstraintF> for [FpGadget<ConstraintF>] {
    type FieldGadget = FpGadget<ConstraintF>;

    #[inline]
    fn to_field_gadget_elements(&self) -> Result<Vec<Self::FieldGadget>, Error> {
        Ok(self.to_vec())
    }
}

impl<ConstraintF: PrimeField> ToConstraintFieldGadget<ConstraintF> for () {
    type FieldGadget = FpGadget<ConstraintF>;

    #[inline]
    fn to_field_gadget_elements(&self) -> Result<Vec<Self::FieldGadget>, Error> {
        Ok(Vec::new())
    }
}

impl<M, ConstraintF, FG> ToConstraintFieldGadget<ConstraintF> for SWPAffineGadget<M, ConstraintF, FG>
    where
        M:              SWModelParameters,
        ConstraintF:    PrimeField,
        FG:             FieldGadget<M::BaseField, ConstraintF> +
        ToConstraintFieldGadget<ConstraintF, FieldGadget = FpGadget<ConstraintF>>,
{
    type FieldGadget = FpGadget<ConstraintF>;

    #[inline]
    fn to_field_gadget_elements(&self) -> Result<Vec<Self::FieldGadget>, Error> {
        let mut x_fe = self.x.to_field_gadget_elements()?;
        let y_fe = self.y.to_field_gadget_elements()?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}