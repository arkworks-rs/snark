use algebra::{
    fields::{
        fp4::{Fp4Parameters, Fp4ParamsWrapper},
        Field, Fp2Parameters,
    },
    Fp2, Fp2ParamsWrapper, PrimeField, SquareRootField,
};
use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{fields::fp2::Fp2Gadget, prelude::*};

impl<P, ConstraintF: PrimeField + SquareRootField> QuadExtParametersGadget<ConstraintF>
    for Fp4ParamsWrapper<P>
where
    P: Fp4Parameters,
    P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    type BaseFieldGadget = Fp2Gadget<P::Fp2Params, ConstraintF>;

    fn mul_base_field_gadget_by_nonresidue<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        fe: &Self::BaseFieldGadget,
    ) -> Result<Self::BaseFieldGadget, SynthesisError> {
        let new_c0 =
            Fp2ParamsWrapper::<P::Fp2Params>::mul_base_field_gadget_by_nonresidue(cs, &fe.c1)?;
        let new_c1 = fe.c0.clone();
        Ok(Self::BaseFieldGadget::new(new_c0, new_c1))
    }

    fn mul_base_field_gadget_by_frobenius_coeff<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        c1: &mut Self::BaseFieldGadget,
        power: usize,
    ) -> Result<(), SynthesisError> {
        c1.c0.mul_by_constant_in_place(
            cs.ns(|| "c1_c0_power"),
            &P::FROBENIUS_COEFF_FP4_C1[power % 4],
        )?;
        c1.c1.mul_by_constant_in_place(
            cs.ns(|| "c1_c1_power"),
            &P::FROBENIUS_COEFF_FP4_C1[power % 4],
        )?;

        Ok(())
    }

    fn cyclotomic_square_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        fe: &QuadExtFieldGadget<Self, ConstraintF>,
    ) -> Result<QuadExtFieldGadget<Self, ConstraintF>, SynthesisError> {
        let c1_squared = fe.c1.square(cs.ns(|| "c1^2"))?;
        let c1_squared_nr =
            Self::mul_base_field_gadget_by_nonresidue(cs.ns(|| "nr * c1^2"), &c1_squared)?;
        let one = Fp2::<P::Fp2Params>::one();

        let c0 = {
            let c1_squared_nr_doubled = c1_squared_nr.double(cs.ns(|| "2(nr*c1^2)"))?;
            c1_squared_nr_doubled.add_constant(cs.ns(|| "2(nr*c1^2) + 1"), &one)?
        };

        let c1 = {
            let c1_plus_c0 = fe.c0.add(cs.ns(|| "c1 + c0"), &fe.c1)?;
            let c1_plus_c0_squared = c1_plus_c0.square(cs.ns(|| "(c1 + c0)^2"))?;
            c1_plus_c0_squared
                .sub(cs.ns(|| "(c1 + c0)^2 - nr*c1^2"), &c1_squared_nr)?
                .sub(cs.ns(|| "(c1 + c0)^2 - nr*c1^2 - c1^2"), &c1_squared)?
                .sub_constant(cs.ns(|| "(c1 + c0)^2 - nr*c1^2 - c1^2 - 1"), &one)?
        };
        Ok(QuadExtFieldGadget::<Self, ConstraintF>::new(c0, c1))
    }
}

pub type Fp4Gadget<P, ConstraintF> = QuadExtFieldGadget<Fp4ParamsWrapper<P>, ConstraintF>;

impl<P, ConstraintF: PrimeField + SquareRootField> Fp4Gadget<P, ConstraintF>
where
    P: Fp4Parameters,
    P::Fp2Params: Fp2Parameters<Fp = ConstraintF>,
{
    #[inline]
    //Mul by an element of the form c0: (a, 0) c1:(c, d)
    pub fn mul_by_023<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let v0 = {
            let v0_c0 = self
                .c0
                .c0
                .mul(cs.ns(|| "self.c0.c0 * other.c0.c0"), &other.c0.c0)?;
            let v0_c1 = self
                .c0
                .c1
                .mul(cs.ns(|| "self.c0.c1 * other.c0.c0"), &other.c0.c0)?;
            Fp2Gadget::<P::Fp2Params, ConstraintF>::new(v0_c0, v0_c1)
        };
        let v1 = self.c1.mul(cs.ns(|| "self.c1 * other.c1"), &other.c1)?;
        let c0 = {
            let non_residue_times_v1 = Fp4ParamsWrapper::<P>::mul_base_field_gadget_by_nonresidue(
                cs.ns(|| "v1 mul_by_nr"),
                &v1,
            )?;
            v0.add(cs.ns(|| "v0 + beta * v1"), &non_residue_times_v1)?
        };
        let c1 = {
            let a0_plus_a1 = self.c0.add(cs.ns(|| "a0 + a1"), &self.c1)?;
            let b0_plus_b1 = other.c0.add(cs.ns(|| "b0 + b1"), &other.c1)?;
            let a0_plus_a1_times_b0_plus_b1 =
                a0_plus_a1.mul(&mut cs.ns(|| "(a0 + a1) * (b0 + b1)"), &b0_plus_b1)?;
            a0_plus_a1_times_b0_plus_b1
                .sub(cs.ns(|| "res - v0"), &v0)?
                .sub(cs.ns(|| "res - v0 - v1"), &v1)?
        };

        Ok(Self::new(c0, c1))
    }
}
