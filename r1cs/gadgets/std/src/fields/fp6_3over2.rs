use algebra::{
    fields::{
        fp2::Fp2Parameters,
        fp6_3over2::{Fp6Parameters, Fp6ParamsWrapper},
        SquareRootField,
    },
    PrimeField,
};
use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{fields::fp2::Fp2Gadget, prelude::*};

impl<P, ConstraintF: PrimeField + SquareRootField> CubicExtParametersGadget<ConstraintF>
    for Fp6ParamsWrapper<P>
where
    P: Fp6Parameters,
    P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    type BaseFieldGadget = Fp2Gadget<P::Fp2Params, ConstraintF>;

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
        c1.mul_by_constant_in_place(cs.ns(|| "c1_power"), &P::FROBENIUS_COEFF_FP6_C1[power % 6])?;
        c2.mul_by_constant_in_place(cs.ns(|| "c2_power"), &P::FROBENIUS_COEFF_FP6_C2[power % 6])?;

        Ok(())
    }
}

pub type Fp6Gadget<P, ConstraintF> = CubicExtFieldGadget<Fp6ParamsWrapper<P>, ConstraintF>;

impl<P, ConstraintF: PrimeField + SquareRootField> Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    #[inline]
    pub fn mul_by_0_c1_0<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        c1: &Fp2Gadget<P::Fp2Params, ConstraintF>,
    ) -> Result<Self, SynthesisError> {
        // Karatsuba multiplication
        // v0 = a0 * b0 = 0

        // v1 = a1 * b1
        let v1 = self.c1.mul(cs.ns(|| "first mul"), c1)?;

        // v2 = a2 * b2 = 0

        let a1_plus_a2 = self.c1.add(cs.ns(|| "a1 + a2"), &self.c2)?;
        let b1_plus_b2 = c1.clone();

        let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;

        // c0 = (NONRESIDUE * ((a1 + a2)*(b1 + b2) - v1 - v2)) + v0
        //    = NONRESIDUE * ((a1 + a2) * b1 - v1)
        let c0 = a1_plus_a2
            .mul(cs.ns(|| "second mul"), &b1_plus_b2)?
            .sub(cs.ns(|| "first sub"), &v1)?
            .mul_by_constant(cs.ns(|| "mul_by_nonresidue"), &P::NONRESIDUE)?;

        // c1 = (a0 + a1) * (b0 + b1) - v0 - v1 + NONRESIDUE * v2
        //    = (a0 + a1) * b1 - v1
        let c1 = a0_plus_a1
            .mul(cs.ns(|| "third mul"), &c1)?
            .sub(cs.ns(|| "second sub"), &v1)?;
        // c2 = (a0 + a2) * (b0 + b2) - v0 - v2 + v1
        //    = v1
        let c2 = v1;
        Ok(Self::new(c0, c1, c2))
    }

    pub fn mul_by_c0_c1_0<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        c0: &Fp2Gadget<P::Fp2Params, ConstraintF>,
        c1: &Fp2Gadget<P::Fp2Params, ConstraintF>,
    ) -> Result<Self, SynthesisError> {
        let v0 = self.c0.mul(cs.ns(|| "v0"), c0)?;
        let v1 = self.c1.mul(cs.ns(|| "v1"), c1)?;
        // v2 = 0.

        let a1_plus_a2 = self.c1.add(cs.ns(|| "a1 + a2"), &self.c2)?;
        let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;
        let a0_plus_a2 = self.c0.add(cs.ns(|| "a0 + a2"), &self.c2)?;

        let b1_plus_b2 = c1.clone();
        let b0_plus_b1 = c0.add(cs.ns(|| "b0 + b1"), &c1)?;
        let b0_plus_b2 = c0.clone();

        let c0 = {
            let cs = &mut cs.ns(|| "c0");
            a1_plus_a2
                .mul(cs.ns(|| "(a1 + a2) * (b1 + b2)"), &b1_plus_b2)?
                .sub(cs.ns(|| "sub v1"), &v1)?
                .mul_by_constant(cs.ns(|| "First mul_by_nonresidue"), &P::NONRESIDUE)?
                .add(cs.ns(|| "add v0"), &v0)?
        };

        let c1 = {
            let cs = &mut cs.ns(|| "c1");
            a0_plus_a1
                .mul(cs.ns(|| "(a0 + a1) * (b0 + b1)"), &b0_plus_b1)?
                .sub(cs.ns(|| "sub v0"), &v0)?
                .sub(cs.ns(|| "sub v1"), &v1)?
        };

        let c2 = {
            a0_plus_a2
                .mul(cs.ns(|| "(a0 + a2) * (b0 + b2)"), &b0_plus_b2)?
                .sub(cs.ns(|| "sub v0"), &v0)?
                .add(cs.ns(|| "add v1"), &v1)?
        };

        Ok(Self::new(c0, c1, c2))
    }
}
