use algebra::{
    fields::{Fp3Parameters, Fp3ParamsWrapper},
    PrimeField, SquareRootField,
};
use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{fields::fp::FpGadget, prelude::*};

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    CubicExtParametersGadget<ConstraintF> for Fp3ParamsWrapper<P>
{
    type BaseFieldGadget = FpGadget<ConstraintF>;

    fn mul_base_field_gadget_by_nonresidue<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        fe: &Self::BaseFieldGadget,
    ) -> Result<Self::BaseFieldGadget, SynthesisError> {
        fe.mul_by_constant(cs, &P::NONRESIDUE)
    }

    fn mul_base_field_gadget_by_frobenius_coeff<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        c1: &mut Self::BaseFieldGadget,
        c2: &mut Self::BaseFieldGadget,
        power: usize,
    ) -> Result<(), SynthesisError> {
        c1.mul_by_constant_in_place(cs.ns(|| "c1_power"), &P::FROBENIUS_COEFF_FP3_C1[power % 3])?;
        c2.mul_by_constant_in_place(cs.ns(|| "c2_power"), &P::FROBENIUS_COEFF_FP3_C2[power % 3])?;

        Ok(())
    }
}

pub type Fp3Gadget<P, ConstraintF> = CubicExtFieldGadget<Fp3ParamsWrapper<P>, ConstraintF>;

impl<P: Fp3Parameters<Fp = ConstraintF>, ConstraintF: PrimeField + SquareRootField>
    Fp3Gadget<P, ConstraintF>
{
    /// Multiply a Fp3Gadget by a Fp gadget.
    #[inline]
    pub fn mul_assign_by_base_field_gadget<CS: ConstraintSystem<ConstraintF>>(
        &mut self,
        mut cs: CS,
        fe: &FpGadget<P::Fp>,
    ) -> Result<&mut Self, SynthesisError> {
        self.c0.mul_in_place(cs.ns(|| "c0"), fe)?;
        self.c1.mul_in_place(cs.ns(|| "c1"), fe)?;
        self.c2.mul_in_place(cs.ns(|| "c2"), fe)?;
        Ok(self)
    }
}
