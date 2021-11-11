use algebra::{
    fields::{Fp2Parameters, Fp2ParamsWrapper, QuadExtParameters},
    PrimeField, SquareRootField,
};
use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{fields::fp::FpGadget, prelude::*};

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    QuadExtParametersGadget<ConstraintF> for Fp2ParamsWrapper<P>
{
    type BaseFieldGadget = FpGadget<ConstraintF>;

    fn mul_base_field_gadget_by_nonresidue<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        fe: &Self::BaseFieldGadget,
    ) -> Result<Self::BaseFieldGadget, SynthesisError> {
        fe.mul_by_constant(cs, &Self::NONRESIDUE)
    }

    fn mul_base_field_gadget_by_frobenius_coeff<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        c1: &mut Self::BaseFieldGadget,
        power: usize,
    ) -> Result<(), SynthesisError> {
        c1.mul_by_constant_in_place(cs, &Self::FROBENIUS_COEFF_C1[power % 2])?;
        Ok(())
    }
}

pub type Fp2Gadget<P, ConstraintF> = QuadExtFieldGadget<Fp2ParamsWrapper<P>, ConstraintF>;

impl<P: Fp2Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    Fp2Gadget<P, ConstraintF>
{
    #[inline]
    pub fn mul_assign_by_base_field_gadget<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        fe: &FpGadget<P::Fp>,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.mul_in_place(cs.ns(|| "compute new_c0"), &fe)?;
        self.c1.mul_in_place(cs.ns(|| "compute new_c1"), &fe)?;
        Ok(self)
    }

    #[inline]
    pub fn mul_by_base_field_constant_in_place<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        fe: &P::Fp,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.mul_by_constant_in_place(cs.ns(|| "c0"), fe)?;
        self.c1.mul_by_constant_in_place(cs.ns(|| "c1"), fe)?;
        Ok(self)
    }

    #[inline]
    pub fn mul_by_base_field_constant<CS: ConstraintSystem<ConstraintF>>(
        &self,
        cs: CS,
        fe: &P::Fp,
    ) -> Result<Self, SynthesisError> {
        let mut result = self.clone();
        result.mul_by_base_field_constant_in_place(cs, fe)?;
        Ok(result)
    }
}
