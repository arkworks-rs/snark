use algebra::{
    curves::mnt4::{
        g2::{AteAdditionCoefficients, AteDoubleCoefficients},
        G1Prepared, G2Prepared, MNT4Parameters,
    },
    Field,
};
use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{
    fields::{fp::FpGadget, fp2::Fp2Gadget, FieldGadget},
    groups::curves::short_weierstrass::AffineGadget,
    pairing::mnt4::PairingGadget,
    prelude::*,
    Vec,
};
use core::borrow::Borrow;

pub type G1Gadget<P> = AffineGadget<
    <P as MNT4Parameters>::G1Parameters,
    <P as MNT4Parameters>::Fp,
    FpGadget<<P as MNT4Parameters>::Fp>,
>;

pub type G2Gadget<P> =
    AffineGadget<<P as MNT4Parameters>::G2Parameters, <P as MNT4Parameters>::Fp, Fp2G<P>>;

#[derive(Derivative)]
#[derivative(Clone(bound = "P: MNT4Parameters"), Debug(bound = "P: MNT4Parameters"))]
pub struct G1PreparedGadget<P: MNT4Parameters> {
    pub x: FpGadget<P::Fp>,
    pub y: FpGadget<P::Fp>,
    pub x_twist: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub y_twist: Fp2Gadget<P::Fp2Params, P::Fp>,
}

impl<P: MNT4Parameters> G1PreparedGadget<P> {
    pub fn get_value(&self) -> Option<G1Prepared<P>> {
        match (
            self.x.get_value(),
            self.y.get_value(),
            self.x_twist.get_value(),
            self.y_twist.get_value(),
        ) {
            (Some(x), Some(y), Some(x_twist), Some(y_twist)) => Some(G1Prepared {
                x,
                y,
                x_twist,
                y_twist,
            }),
            _ => None,
        }
    }

    pub fn from_affine<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        q: &G1Gadget<P>,
    ) -> Result<Self, SynthesisError> {
        let x_twist = Fp2Gadget::new(
            q.x.mul_by_constant(cs.ns(|| "g1.x * twist.c0"), &P::TWIST.c0)?,
            q.x.mul_by_constant(cs.ns(|| "g1.x * twist.c1"), &P::TWIST.c1)?,
        );
        let y_twist = Fp2Gadget::new(
            q.y.mul_by_constant(cs.ns(|| "g1.y * twist.c0"), &P::TWIST.c0)?,
            q.y.mul_by_constant(cs.ns(|| "g1.y * twist.c1"), &P::TWIST.c1)?,
        );
        Ok(G1PreparedGadget {
            x: q.x.clone(),
            y: q.y.clone(),
            x_twist,
            y_twist,
        })
    }
}

impl<P: MNT4Parameters> AllocGadget<G1Prepared<P>, P::Fp> for G1PreparedGadget<P> {
    fn alloc_constant<T, CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<G1Prepared<P>>,
    {
        let obj = t.borrow();
        let x = FpGadget::<P::Fp>::alloc_constant(&mut cs.ns(|| "alloc_x"), &obj.x)?;
        let y = FpGadget::<P::Fp>::alloc_constant(&mut cs.ns(|| "alloc_y"), &obj.y)?;
        let x_twist = Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(
            &mut cs.ns(|| "alloc_x_twist"),
            &obj.x_twist,
        )?;
        let y_twist = Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(
            &mut cs.ns(|| "alloc_y_twist"),
            &obj.y_twist,
        )?;
        Ok(G1PreparedGadget {
            x,
            y,
            x_twist,
            y_twist,
        })
    }

    fn alloc<F, T, CS: ConstraintSystem<P::Fp>>(_cs: CS, _f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<G1Prepared<P>>,
    {
        todo!()
    }

    fn alloc_input<F, T, CS: ConstraintSystem<P::Fp>>(
        _cs: CS,
        _f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<G1Prepared<P>>,
    {
        todo!()
    }
}

impl<P: MNT4Parameters> ToBytesGadget<P::Fp> for G1PreparedGadget<P> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x = self.x.to_bytes(&mut cs.ns(|| "x to bytes"))?;
        let mut y = self.y.to_bytes(&mut cs.ns(|| "y to bytes"))?;
        let mut x_twist = self.x_twist.to_bytes(&mut cs.ns(|| "x_twist to bytes"))?;
        let mut y_twist = self.y_twist.to_bytes(&mut cs.ns(|| "y_twist to bytes"))?;

        x.append(&mut y);
        x.append(&mut x_twist);
        x.append(&mut y_twist);
        Ok(x)
    }

    fn to_non_unique_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x = self.x.to_non_unique_bytes(&mut cs.ns(|| "x to bytes"))?;
        let mut y = self.y.to_non_unique_bytes(&mut cs.ns(|| "y to bytes"))?;
        let mut x_twist = self
            .x_twist
            .to_non_unique_bytes(&mut cs.ns(|| "x_twist to bytes"))?;
        let mut y_twist = self
            .y_twist
            .to_non_unique_bytes(&mut cs.ns(|| "y_twist to bytes"))?;

        x.append(&mut y);
        x.append(&mut x_twist);
        x.append(&mut y_twist);
        Ok(x)
    }
}

type Fp2G<P> = Fp2Gadget<<P as MNT4Parameters>::Fp2Params, <P as MNT4Parameters>::Fp>;
#[derive(Derivative)]
#[derivative(Clone(bound = "P: MNT4Parameters"), Debug(bound = "P: MNT4Parameters"))]
pub struct G2PreparedGadget<P: MNT4Parameters> {
    pub x: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub y: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub x_over_twist: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub y_over_twist: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub double_coefficients: Vec<AteDoubleCoefficientsGadget<P>>,
    pub addition_coefficients: Vec<AteAdditionCoefficientsGadget<P>>,
}

impl<P: MNT4Parameters> AllocGadget<G2Prepared<P>, P::Fp> for G2PreparedGadget<P> {
    fn alloc_constant<T, CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<G2Prepared<P>>,
    {
        let obj = t.borrow();
        let x = Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(&mut cs.ns(|| "alloc_x"), &obj.x)?;
        let y = Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(&mut cs.ns(|| "alloc_y"), &obj.y)?;
        let x_over_twist = Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(
            &mut cs.ns(|| "alloc_x_over_twist"),
            &obj.x_over_twist,
        )?;
        let y_over_twist = Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(
            &mut cs.ns(|| "alloc_y_over_twist"),
            &obj.y_over_twist,
        )?;
        let mut double_coefficients = Vec::<AteDoubleCoefficientsGadget<P>>::new();
        for (i, item) in obj.double_coefficients.iter().enumerate() {
            double_coefficients.push(AteDoubleCoefficientsGadget::<P>::alloc_constant(
                &mut cs.ns(|| format!("alloc_double_coefficients_{}", i)),
                item,
            )?);
        }
        let mut addition_coefficients = Vec::<AteAdditionCoefficientsGadget<P>>::new();
        for (i, item) in obj.addition_coefficients.iter().enumerate() {
            addition_coefficients.push(AteAdditionCoefficientsGadget::<P>::alloc_constant(
                &mut cs.ns(|| format!("alloc_addition_coefficients_{}", i)),
                item,
            )?);
        }

        Ok(G2PreparedGadget {
            x,
            y,
            x_over_twist,
            y_over_twist,
            double_coefficients,
            addition_coefficients,
        })
    }

    fn alloc<F, T, CS: ConstraintSystem<P::Fp>>(_cs: CS, _f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<G2Prepared<P>>,
    {
        todo!()
    }

    fn alloc_input<F, T, CS: ConstraintSystem<P::Fp>>(
        _cs: CS,
        _f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<G2Prepared<P>>,
    {
        todo!()
    }
}

impl<P: MNT4Parameters> ToBytesGadget<P::Fp> for G2PreparedGadget<P> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x = self.x.to_bytes(&mut cs.ns(|| "x to bytes"))?;
        let mut y = self.y.to_bytes(&mut cs.ns(|| "y to bytes"))?;
        let mut x_over_twist = self
            .x_over_twist
            .to_bytes(&mut cs.ns(|| "x_over_twist to bytes"))?;
        let mut y_over_twist = self
            .y_over_twist
            .to_bytes(&mut cs.ns(|| "y_over_twist to bytes"))?;

        x.append(&mut y);
        x.append(&mut x_over_twist);
        x.append(&mut y_over_twist);

        for (i, coeff) in self.double_coefficients.iter().enumerate() {
            x.extend_from_slice(&coeff.to_bytes(cs.ns(|| format!("double_coefficients {}", i)))?);
        }
        for (i, coeff) in self.addition_coefficients.iter().enumerate() {
            x.extend_from_slice(&coeff.to_bytes(cs.ns(|| format!("addition_coefficients {}", i)))?);
        }
        Ok(x)
    }

    fn to_non_unique_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut x = self.x.to_non_unique_bytes(&mut cs.ns(|| "x to bytes"))?;
        let mut y = self.y.to_non_unique_bytes(&mut cs.ns(|| "y to bytes"))?;
        let mut x_over_twist = self
            .x_over_twist
            .to_non_unique_bytes(&mut cs.ns(|| "x_over_twist to bytes"))?;
        let mut y_over_twist = self
            .y_over_twist
            .to_non_unique_bytes(&mut cs.ns(|| "y_over_twist to bytes"))?;

        x.append(&mut y);
        x.append(&mut x_over_twist);
        x.append(&mut y_over_twist);

        for (i, coeff) in self.double_coefficients.iter().enumerate() {
            x.extend_from_slice(
                &coeff.to_non_unique_bytes(cs.ns(|| format!("double_coefficients {}", i)))?,
            );
        }
        for (i, coeff) in self.addition_coefficients.iter().enumerate() {
            x.extend_from_slice(
                &coeff.to_non_unique_bytes(cs.ns(|| format!("addition_coefficients {}", i)))?,
            );
        }
        Ok(x)
    }
}

impl<P: MNT4Parameters> G2PreparedGadget<P> {
    pub fn get_value(&self) -> Option<G2Prepared<P>> {
        match (
            self.x.get_value(),
            self.y.get_value(),
            self.x_over_twist.get_value(),
            self.y_over_twist.get_value(),
            self.double_coefficients
                .iter()
                .map(|coeff| coeff.get_value())
                .collect::<Option<Vec<AteDoubleCoefficients<P>>>>(),
            self.addition_coefficients
                .iter()
                .map(|coeff| coeff.get_value())
                .collect::<Option<Vec<AteAdditionCoefficients<P>>>>(),
        ) {
            (
                Some(x),
                Some(y),
                Some(x_over_twist),
                Some(y_over_twist),
                Some(double_coefficients),
                Some(addition_coefficients),
            ) => Some(G2Prepared {
                x,
                y,
                x_over_twist,
                y_over_twist,
                double_coefficients,
                addition_coefficients,
            }),
            _ => None,
        }
    }

    pub fn from_affine<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        q: &G2Gadget<P>,
    ) -> Result<Self, SynthesisError> {
        let twist_inv = P::TWIST.inverse().unwrap();

        let mut g2p = G2PreparedGadget {
            x: q.x.clone(),
            y: q.y.clone(),
            x_over_twist: q.x.mul_by_constant(cs.ns(|| "x over twist"), &twist_inv)?,
            y_over_twist: q.y.mul_by_constant(cs.ns(|| "y over twist"), &twist_inv)?,
            double_coefficients: vec![],
            addition_coefficients: vec![],
        };

        let fp2_one = Fp2G::<P>::one(cs.ns(|| "one"))?;
        let mut r = G2ProjectiveExtendedGadget {
            x: q.x.clone(),
            y: q.y.clone(),
            z: fp2_one.clone(),
            t: fp2_one,
        };

        for (idx, value) in P::ATE_LOOP_COUNT.iter().rev().enumerate() {
            let mut tmp = *value;
            let skip_extraneous_bits = 64 - value.leading_zeros();
            let mut v = Vec::with_capacity(16);
            for i in 0..64 {
                if idx == 0 && (i == 0 || i >= skip_extraneous_bits) {
                    continue;
                }
                v.push(tmp & 1 == 1);
                tmp >>= 1;
            }

            let mut cs = cs.ns(|| format!("ate loop iteration {}", idx));

            for (j, bit) in v.iter().rev().enumerate() {
                let (r2, coeff) = PairingGadget::<P>::doubling_step_for_flipped_miller_loop(
                    cs.ns(|| format!("doubling step {}", j)),
                    &r,
                )?;
                g2p.double_coefficients.push(coeff);
                r = r2;

                if *bit {
                    let (r2, coeff) =
                        PairingGadget::<P>::mixed_addition_step_for_flipped_miller_loop(
                            cs.ns(|| format!("mixed addition step {}", j)),
                            &q.x,
                            &q.y,
                            &r,
                        )?;
                    g2p.addition_coefficients.push(coeff);
                    r = r2;
                }

                tmp >>= 1;
            }
        }

        if P::ATE_IS_LOOP_COUNT_NEG {
            let rz_inv = r.z.inverse(cs.ns(|| "inverse r.z"))?;
            let rz2_inv = rz_inv.square(cs.ns(|| "rz_inv^2"))?;
            let rz3_inv = rz_inv.mul(cs.ns(|| "rz_inv * rz_inv^2"), &rz2_inv)?;

            let minus_r_affine_x = r.x.mul(cs.ns(|| "r.x * rz2_inv"), &rz2_inv)?;
            let minus_r_affine_y =
                r.y.negate(cs.ns(|| "-r.y"))?
                    .mul(cs.ns(|| "-r.y * rz3_inv"), &rz3_inv)?;

            let add_result = PairingGadget::<P>::mixed_addition_step_for_flipped_miller_loop(
                cs.ns(|| "mixed_addition step"),
                &minus_r_affine_x,
                &minus_r_affine_y,
                &r,
            )?;
            g2p.addition_coefficients.push(add_result.1);
        }

        Ok(g2p)
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = "P: MNT4Parameters"), Debug(bound = "P: MNT4Parameters"))]
pub struct AteDoubleCoefficientsGadget<P: MNT4Parameters> {
    pub c_h: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub c_4c: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub c_j: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub c_l: Fp2Gadget<P::Fp2Params, P::Fp>,
}

impl<P: MNT4Parameters> AllocGadget<AteDoubleCoefficients<P>, P::Fp>
    for AteDoubleCoefficientsGadget<P>
{
    fn alloc_constant<T, CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<AteDoubleCoefficients<P>>,
    {
        let obj = t.borrow();
        let c_h =
            Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(&mut cs.ns(|| "alloc_c_h"), &obj.c_h)?;
        let c_4c = Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(
            &mut cs.ns(|| "alloc_c_4c"),
            &obj.c_4c,
        )?;
        let c_j =
            Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(&mut cs.ns(|| "alloc_c_j"), &obj.c_j)?;
        let c_l =
            Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(&mut cs.ns(|| "alloc_c_l"), &obj.c_l)?;

        Ok(AteDoubleCoefficientsGadget {
            c_h,
            c_4c,
            c_j,
            c_l,
        })
    }

    fn alloc<F, T, CS: ConstraintSystem<P::Fp>>(_cs: CS, _f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<AteDoubleCoefficients<P>>,
    {
        todo!()
    }

    fn alloc_input<F, T, CS: ConstraintSystem<P::Fp>>(
        _cs: CS,
        _f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<AteDoubleCoefficients<P>>,
    {
        todo!()
    }
}

impl<P: MNT4Parameters> ToBytesGadget<P::Fp> for AteDoubleCoefficientsGadget<P> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c_h = self.c_h.to_bytes(&mut cs.ns(|| "c_h to bytes"))?;
        let mut c_4c = self.c_4c.to_bytes(&mut cs.ns(|| "c_4c to bytes"))?;
        let mut c_j = self.c_j.to_bytes(&mut cs.ns(|| "c_j to bytes"))?;
        let mut c_l = self.c_l.to_bytes(&mut cs.ns(|| "c_l to bytes"))?;

        c_h.append(&mut c_4c);
        c_h.append(&mut c_j);
        c_h.append(&mut c_l);
        Ok(c_h)
    }

    fn to_non_unique_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c_h = self
            .c_h
            .to_non_unique_bytes(&mut cs.ns(|| "c_h to bytes"))?;
        let mut c_4c = self
            .c_4c
            .to_non_unique_bytes(&mut cs.ns(|| "c_4c to bytes"))?;
        let mut c_j = self
            .c_j
            .to_non_unique_bytes(&mut cs.ns(|| "c_j to bytes"))?;
        let mut c_l = self
            .c_l
            .to_non_unique_bytes(&mut cs.ns(|| "c_l to bytes"))?;

        c_h.append(&mut c_4c);
        c_h.append(&mut c_j);
        c_h.append(&mut c_l);
        Ok(c_h)
    }
}

impl<P: MNT4Parameters> AteDoubleCoefficientsGadget<P> {
    pub fn get_value(&self) -> Option<AteDoubleCoefficients<P>> {
        match (
            self.c_h.get_value(),
            self.c_4c.get_value(),
            self.c_j.get_value(),
            self.c_l.get_value(),
        ) {
            (Some(c_h), Some(c_4c), Some(c_j), Some(c_l)) => Some(AteDoubleCoefficients {
                c_h,
                c_4c,
                c_j,
                c_l,
            }),
            _ => None,
        }
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = "P: MNT4Parameters"), Debug(bound = "P: MNT4Parameters"))]
pub struct AteAdditionCoefficientsGadget<P: MNT4Parameters> {
    pub c_l1: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub c_rz: Fp2Gadget<P::Fp2Params, P::Fp>,
}

impl<P: MNT4Parameters> AllocGadget<AteAdditionCoefficients<P>, P::Fp>
    for AteAdditionCoefficientsGadget<P>
{
    fn alloc_constant<T, CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        t: T,
    ) -> Result<Self, SynthesisError>
    where
        T: Borrow<AteAdditionCoefficients<P>>,
    {
        let obj = t.borrow();
        let c_l1 = Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(
            &mut cs.ns(|| "alloc_c_l1"),
            &obj.c_l1,
        )?;
        let c_rz = Fp2Gadget::<P::Fp2Params, P::Fp>::alloc_constant(
            &mut cs.ns(|| "alloc_c_rz"),
            &obj.c_rz,
        )?;
        Ok(AteAdditionCoefficientsGadget { c_l1, c_rz })
    }

    fn alloc<F, T, CS: ConstraintSystem<P::Fp>>(_cs: CS, _f: F) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<AteAdditionCoefficients<P>>,
    {
        todo!()
    }

    fn alloc_input<F, T, CS: ConstraintSystem<P::Fp>>(
        _cs: CS,
        _f: F,
    ) -> Result<Self, SynthesisError>
    where
        F: FnOnce() -> Result<T, SynthesisError>,
        T: Borrow<AteAdditionCoefficients<P>>,
    {
        todo!()
    }
}

impl<P: MNT4Parameters> ToBytesGadget<P::Fp> for AteAdditionCoefficientsGadget<P> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c_l1 = self.c_l1.to_bytes(&mut cs.ns(|| "c_l1 to bytes"))?;
        let mut c_rz = self.c_rz.to_bytes(&mut cs.ns(|| "c_rz to bytes"))?;

        c_l1.append(&mut c_rz);
        Ok(c_l1)
    }

    fn to_non_unique_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut c_l1 = self
            .c_l1
            .to_non_unique_bytes(&mut cs.ns(|| "c_l1 to bytes"))?;
        let mut c_rz = self
            .c_rz
            .to_non_unique_bytes(&mut cs.ns(|| "c_rz to bytes"))?;

        c_l1.append(&mut c_rz);
        Ok(c_l1)
    }
}

impl<P: MNT4Parameters> AteAdditionCoefficientsGadget<P> {
    pub fn get_value(&self) -> Option<AteAdditionCoefficients<P>> {
        match (self.c_l1.get_value(), self.c_rz.get_value()) {
            (Some(c_l1), Some(c_rz)) => Some(AteAdditionCoefficients { c_l1, c_rz }),
            _ => None,
        }
    }
}

pub struct G2ProjectiveExtendedGadget<P: MNT4Parameters> {
    pub x: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub y: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub z: Fp2Gadget<P::Fp2Params, P::Fp>,
    pub t: Fp2Gadget<P::Fp2Params, P::Fp>,
}
