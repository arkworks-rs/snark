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

//TODO: Maybe generalize to FieldGadget<F, ConstraintF> instead of FpGadget<ConstraintF> ?
/// Types that can be converted to a vector of `FpGadget` elements.
pub trait ToConstraintFieldGadget<ConstraintF: PrimeField> {
    fn to_field_gadget_elements(&self) -> Result<Vec<FpGadget<ConstraintF>>, Error>;
}

impl<ConstraintF: PrimeField> ToConstraintFieldGadget<ConstraintF> for FpGadget<ConstraintF> {
    fn to_field_gadget_elements(&self) -> Result<Vec<FpGadget<ConstraintF>>, Error> {
        Ok(vec![self.clone()])
    }
}

impl<ConstraintF: PrimeField> ToConstraintFieldGadget<ConstraintF> for [FpGadget<ConstraintF>] {
    #[inline]
    fn to_field_gadget_elements(&self) -> Result<Vec<FpGadget<ConstraintF>>, Error> {
        Ok(self.to_vec())
    }
}

impl<ConstraintF: PrimeField> ToConstraintFieldGadget<ConstraintF> for () {
    #[inline]
    fn to_field_gadget_elements(&self) -> Result<Vec<FpGadget<ConstraintF>>, Error> {
        Ok(Vec::new())
    }
}

impl<M, ConstraintF, FG> ToConstraintFieldGadget<ConstraintF> for SWPAffineGadget<M, ConstraintF, FG>
    where
        M:              SWModelParameters,
        ConstraintF:    PrimeField,
        FG:             FieldGadget<M::BaseField, ConstraintF> + ToConstraintFieldGadget<ConstraintF>,
{
    #[inline]
    fn to_field_gadget_elements(&self) -> Result<Vec<FpGadget<ConstraintF>>, Error> {
        let mut x_fe = self.x.to_field_gadget_elements()?;
        let y_fe = self.y.to_field_gadget_elements()?;
        x_fe.extend_from_slice(&y_fe);
        Ok(x_fe)
    }
}