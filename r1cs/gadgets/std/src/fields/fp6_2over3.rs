use algebra::{
    fields::{
        fp3::{Fp3Parameters, Fp3ParamsWrapper},
        fp6_2over3::{Fp6Parameters, Fp6ParamsWrapper},
    },
    PrimeField, SquareRootField,
};
use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{fields::fp3::Fp3Gadget, prelude::*};

impl<P, ConstraintF: PrimeField + SquareRootField> QuadExtParametersGadget<ConstraintF>
    for Fp6ParamsWrapper<P>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    type BaseFieldGadget = Fp3Gadget<P::Fp3Params, ConstraintF>;

    fn mul_base_field_gadget_by_nonresidue<CS: ConstraintSystem<ConstraintF>>(
        cs: CS,
        fe: &Self::BaseFieldGadget,
    ) -> Result<Self::BaseFieldGadget, SynthesisError> {
        let new_c0 =
            Fp3ParamsWrapper::<P::Fp3Params>::mul_base_field_gadget_by_nonresidue(cs, &fe.c2)?;
        let new_c1 = fe.c0.clone();
        let new_c2 = fe.c1.clone();
        Ok(Fp3Gadget::<P::Fp3Params, ConstraintF>::new(
            new_c0, new_c1, new_c2,
        ))
    }

    fn mul_base_field_gadget_by_frobenius_coeff<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        c1: &mut Self::BaseFieldGadget,
        power: usize,
    ) -> Result<(), SynthesisError> {
        c1.c0
            .mul_by_constant_in_place(cs.ns(|| "mul1"), &P::FROBENIUS_COEFF_FP6_C1[power % 6])?;
        c1.c1
            .mul_by_constant_in_place(cs.ns(|| "mul2"), &P::FROBENIUS_COEFF_FP6_C1[power % 6])?;
        c1.c2
            .mul_by_constant_in_place(cs.ns(|| "mul3"), &P::FROBENIUS_COEFF_FP6_C1[power % 6])?;
        Ok(())
    }

    fn cyclotomic_square_gadget<CS: ConstraintSystem<ConstraintF>>(
        mut cs: CS,
        fe: &QuadExtFieldGadget<Self, ConstraintF>,
    ) -> Result<QuadExtFieldGadget<Self, ConstraintF>, SynthesisError> {
        let mut result = QuadExtFieldGadget::<Self, ConstraintF>::zero(cs.ns(|| "alloc result"))?;
        let fp2_nr = <P::Fp3Params as Fp3Parameters>::NONRESIDUE;

        let z0 = &fe.c0.c0;
        let z4 = &fe.c0.c1;
        let z3 = &fe.c0.c2;
        let z2 = &fe.c1.c0;
        let z1 = &fe.c1.c1;
        let z5 = &fe.c1.c2;

        // t0 + t1*y = (z0 + z1*y)^2 = a^2
        let tmp = z0.mul(cs.ns(|| "first mul"), &z1)?;
        let t0 = {
            // (z0 + &z1) * &(z0 + &(fp2_nr * &z1)) - &tmp - &(tmp * &fp2_nr);
            let mut cs = cs.ns(|| "t0");
            let tmp1 = z0.add(cs.ns(|| "tmp1"), &z1)?;
            let tmp2 = z1
                .mul_by_constant(cs.ns(|| "tmp2.1"), &fp2_nr)?
                .add(cs.ns(|| "tmp2.2"), &z0)?;
            let tmp4 = tmp
                .mul_by_constant(cs.ns(|| "tmp4.1"), &fp2_nr)?
                .add(cs.ns(|| "tmp4.2"), &tmp)?;
            tmp1.mul(cs.ns(|| "tmp3.1"), &tmp2)?
                .sub(cs.ns(|| "tmp3.2"), &tmp4)?
        };
        let t1 = tmp.double(cs.ns(|| "t1"))?;

        // t2 + t3*y = (z2 + z3*y)^2 = b^2
        let tmp = z2.mul(cs.ns(|| "second mul"), &z3)?;
        let t2 = {
            // (z2 + &z3) * &(z2 + &(fp2_nr * &z3)) - &tmp - &(tmp * &fp2_nr);
            let mut cs = cs.ns(|| "t2");
            let tmp1 = z2.add(cs.ns(|| "tmp1"), &z3)?;
            let tmp2 = z3
                .mul_by_constant(cs.ns(|| "tmp2.1"), &fp2_nr)?
                .add(cs.ns(|| "tmp2.2"), &z2)?;
            let tmp4 = tmp
                .mul_by_constant(cs.ns(|| "tmp4.1"), &fp2_nr)?
                .add(cs.ns(|| "tmp4.2"), &tmp)?;
            tmp1.mul(cs.ns(|| "tmp3.1"), &tmp2)?
                .sub(cs.ns(|| "tmp3.2"), &tmp4)?
        };
        let t3 = tmp.double(cs.ns(|| "t3"))?;

        // t4 + t5*y = (z4 + z5*y)^2 = c^2
        let tmp = z4.mul(cs.ns(|| "third mul"), &z5)?;
        let t4 = {
            // (z4 + &z5) * &(z4 + &(fp2_nr * &z5)) - &tmp - &(tmp * &fp2_nr);
            let mut cs = cs.ns(|| "t4");
            let tmp1 = z4.add(cs.ns(|| "tmp1"), &z5)?;
            let tmp2 = z5
                .mul_by_constant(cs.ns(|| "tmp2.1"), &fp2_nr)?
                .add(cs.ns(|| "tmp2.2"), &z4)?;
            let tmp4 = tmp
                .mul_by_constant(cs.ns(|| "tmp4.1"), &fp2_nr)?
                .add(cs.ns(|| "tmp4.2"), &tmp)?;
            tmp1.mul(cs.ns(|| "tmp3.1"), &tmp2)?
                .sub(cs.ns(|| "tmp3.2"), &tmp4)?
        };
        let t5 = tmp.double(cs.ns(|| "t5"))?;

        // for A

        // z0 = 3 * t0 - 2 * z0
        result.c0.c0 = {
            let mut cs = cs.ns(|| "result.c0.c0");
            t0.sub(cs.ns(|| "1"), &z0)?
                .double(cs.ns(|| "2"))?
                .add(cs.ns(|| "3"), &t0)?
        };

        // z1 = 3 * t1 + 2 * z1
        result.c1.c1 = {
            let mut cs = cs.ns(|| "result.c1.c1");
            t1.add(cs.ns(|| "1"), &z1)?
                .double(cs.ns(|| "2"))?
                .add(cs.ns(|| "3"), &t1)?
        };

        // for B

        // z2 = 3 * (xi * t5) + 2 * z2
        result.c1.c0 = {
            let mut cs = cs.ns(|| "result.c1.c0");
            let tmp = t5.mul_by_constant(cs.ns(|| "1"), &fp2_nr)?;
            z2.add(cs.ns(|| "2"), &tmp)?
                .double(cs.ns(|| "3"))?
                .add(cs.ns(|| "4"), &tmp)?
        };

        // z3 = 3 * t4 - 2 * z3
        result.c0.c2 = {
            let mut cs = cs.ns(|| "result.c0.c2");
            t4.sub(cs.ns(|| "1"), &z3)?
                .double(cs.ns(|| "2"))?
                .add(cs.ns(|| "3"), &t4)?
        };

        // for C

        // z4 = 3 * t2 - 2 * z4
        result.c0.c1 = {
            let mut cs = cs.ns(|| "result.c0.c1");
            t2.sub(cs.ns(|| "1"), &z4)?
                .double(cs.ns(|| "2"))?
                .add(cs.ns(|| "3"), &t2)?
        };

        // z5 = 3 * t3 + 2 * z5
        result.c1.c2 = {
            let mut cs = cs.ns(|| "result.c1.c2");
            t3.add(cs.ns(|| "1"), &z5)?
                .double(cs.ns(|| "2"))?
                .add(cs.ns(|| "3"), &t3)?
        };

        Ok(result)
    }
}

pub type Fp6Gadget<P, ConstraintF> = QuadExtFieldGadget<Fp6ParamsWrapper<P>, ConstraintF>;

impl<P, ConstraintF: PrimeField + SquareRootField> Fp6Gadget<P, ConstraintF>
where
    P: Fp6Parameters,
    P::Fp3Params: Fp3Parameters<Fp = ConstraintF>,
{
    #[inline]
    pub fn mul_by_2345<CS: ConstraintSystem<ConstraintF>>(
        &self,
        mut cs: CS,
        other: &Self,
    ) -> Result<Self, SynthesisError> {
        let v0 = {
            let t = Fp3ParamsWrapper::<P::Fp3Params>::mul_base_field_gadget_by_nonresidue(
                cs.ns(|| "other.c0.c2 * nr"),
                &other.c0.c2,
            )?;
            let c0 = self.c0.c1.mul(cs.ns(|| "compute v0_c0"), &t)?;
            let c1 = self.c0.c2.mul(cs.ns(|| "compute v0_c1"), &t)?;
            let c2 = self.c0.c0.mul(cs.ns(|| "compute v0_c2"), &other.c0.c2)?;
            Fp3Gadget::<P::Fp3Params, ConstraintF>::new(c0, c1, c2)
        };
        let v1 = self.c1.mul(cs.ns(|| "compute v1"), &other.c1)?;
        let beta_v1 =
            Fp6ParamsWrapper::<P>::mul_base_field_gadget_by_nonresidue(cs.ns(|| "v1*nr"), &v1)?;

        let c0 = v0.add(cs.ns(|| "compute result c0"), &beta_v1)?;
        let c1 = {
            let self_c0_plus_c1 = self.c0.add(cs.ns(|| "self.c0 + self.c1"), &self.c1)?;
            let other_c0_plus_c1 = other.c0.add(cs.ns(|| "other.c0 + other.c1"), &other.c1)?;
            self_c0_plus_c1
                .mul(
                    cs.ns(|| "(self.c0 + self.c1)*(other.c0 + other.c1)"),
                    &other_c0_plus_c1,
                )?
                .sub(
                    cs.ns(|| "(self.c0 + self.c1)*(other.c0 + other.c1) - v0"),
                    &v0,
                )?
                .sub(
                    cs.ns(|| "(self.c0 + self.c1)*(other.c0 + other.c1) - v0 - v1"),
                    &v1,
                )?
        };
        Ok(Self::new(c0, c1))
    }
}
