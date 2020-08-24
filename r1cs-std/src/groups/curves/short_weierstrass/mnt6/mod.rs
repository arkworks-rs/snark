use algebra::{
    curves::mnt6::{
        g2::{AteAdditionCoefficients, AteDoubleCoefficients},
        G1Prepared, G2Prepared, MNT6Parameters,
    },
    Field,
};
use r1cs_core::{Namespace, SynthesisError};

use crate::{
    fields::{fp::FpVar, fp3::Fp3Var, FieldVar},
    groups::curves::short_weierstrass::ProjectiveVar,
    pairing::mnt6::PairingVar,
    prelude::*,
    Vec,
};
use core::borrow::Borrow;

pub type G1Var<P> =
    ProjectiveVar<<P as MNT6Parameters>::G1Parameters, FpVar<<P as MNT6Parameters>::Fp>>;

pub type G2Var<P> = ProjectiveVar<<P as MNT6Parameters>::G2Parameters, Fp3G<P>>;

#[derive(Derivative)]
#[derivative(Clone(bound = "P: MNT6Parameters"), Debug(bound = "P: MNT6Parameters"))]
pub struct G1PreparedVar<P: MNT6Parameters> {
    pub x: FpVar<P::Fp>,
    pub y: FpVar<P::Fp>,
    pub x_twist: Fp3Var<P::Fp3Params>,
    pub y_twist: Fp3Var<P::Fp3Params>,
}

impl<P: MNT6Parameters> G1PreparedVar<P> {
    pub fn value(&self) -> Result<G1Prepared<P>, SynthesisError> {
        let x = self.x.value()?;
        let y = self.y.value()?;
        let x_twist = self.x_twist.value()?;
        let y_twist = self.y_twist.value()?;
        Ok(G1Prepared {
            x,
            y,
            x_twist,
            y_twist,
        })
    }

    pub fn from_group_var(q: &G1Var<P>) -> Result<Self, SynthesisError> {
        let q = q.to_affine()?;
        let zero = FpVar::<P::Fp>::zero();
        let x_twist = Fp3Var::new(q.x.clone(), zero.clone(), zero.clone()) * P::TWIST;
        let y_twist = Fp3Var::new(q.y.clone(), zero.clone(), zero.clone()) * P::TWIST;
        let result = G1PreparedVar {
            x: q.x.clone(),
            y: q.y.clone(),
            x_twist,
            y_twist,
        };
        Ok(result)
    }
}

impl<P: MNT6Parameters> AllocVar<G1Prepared<P>, P::Fp> for G1PreparedVar<P> {
    fn new_variable<T: Borrow<G1Prepared<P>>>(
        cs: impl Into<Namespace<P::Fp>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let g1_prep = f().map(|b| *b.borrow());

        let x = FpVar::new_variable(cs.ns("x"), || g1_prep.map(|g| g.x), mode)?;
        let y = FpVar::new_variable(cs.ns("y"), || g1_prep.map(|g| g.y), mode)?;
        let x_twist = Fp3Var::new_variable(cs.ns("x_twist"), || g1_prep.map(|g| g.x_twist), mode)?;
        let y_twist = Fp3Var::new_variable(cs.ns("y_twist"), || g1_prep.map(|g| g.y_twist), mode)?;
        Ok(Self {
            x,
            y,
            x_twist,
            y_twist,
        })
    }
}

impl<P: MNT6Parameters> ToBytesGadget<P::Fp> for G1PreparedVar<P> {
    #[inline]
    fn to_bytes(&self) -> Result<Vec<UInt8<P::Fp>>, SynthesisError> {
        let mut x = self.x.to_bytes()?;
        let mut y = self.y.to_bytes()?;
        let mut x_twist = self.x_twist.to_bytes()?;
        let mut y_twist = self.y_twist.to_bytes()?;

        x.append(&mut y);
        x.append(&mut x_twist);
        x.append(&mut y_twist);
        Ok(x)
    }

    fn to_non_unique_bytes(&self) -> Result<Vec<UInt8<P::Fp>>, SynthesisError> {
        let mut x = self.x.to_non_unique_bytes()?;
        let mut y = self.y.to_non_unique_bytes()?;
        let mut x_twist = self.x_twist.to_non_unique_bytes()?;
        let mut y_twist = self.y_twist.to_non_unique_bytes()?;

        x.append(&mut y);
        x.append(&mut x_twist);
        x.append(&mut y_twist);
        Ok(x)
    }
}

type Fp3G<P> = Fp3Var<<P as MNT6Parameters>::Fp3Params>;
#[derive(Derivative)]
#[derivative(Clone(bound = "P: MNT6Parameters"), Debug(bound = "P: MNT6Parameters"))]
pub struct G2PreparedVar<P: MNT6Parameters> {
    pub x: Fp3Var<P::Fp3Params>,
    pub y: Fp3Var<P::Fp3Params>,
    pub x_over_twist: Fp3Var<P::Fp3Params>,
    pub y_over_twist: Fp3Var<P::Fp3Params>,
    pub double_coefficients: Vec<AteDoubleCoefficientsVar<P>>,
    pub addition_coefficients: Vec<AteAdditionCoefficientsVar<P>>,
}

impl<P: MNT6Parameters> AllocVar<G2Prepared<P>, P::Fp> for G2PreparedVar<P> {
    fn new_variable<T: Borrow<G2Prepared<P>>>(
        cs: impl Into<Namespace<P::Fp>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let g2_prep = f().map(|b| b.borrow().clone());
        let g2 = g2_prep.as_ref().map_err(|e| *e);

        let x = Fp3Var::new_variable(cs.ns("x"), || g2.map(|g| g.x), mode)?;
        let y = Fp3Var::new_variable(cs.ns("y"), || g2.map(|g| g.y), mode)?;
        let x_over_twist =
            Fp3Var::new_variable(cs.ns("x_over_twist"), || g2.map(|g| g.x_over_twist), mode)?;
        let y_over_twist =
            Fp3Var::new_variable(cs.ns("y_over_twist"), || g2.map(|g| g.y_over_twist), mode)?;
        let double_coefficients = Vec::new_variable(
            cs.ns("double coeffs"),
            || g2.map(|g| g.double_coefficients.clone()),
            mode,
        )?;
        let addition_coefficients = Vec::new_variable(
            cs.ns("add coeffs"),
            || g2.map(|g| g.addition_coefficients.clone()),
            mode,
        )?;
        Ok(Self {
            x,
            y,
            x_over_twist,
            y_over_twist,
            double_coefficients,
            addition_coefficients,
        })
    }
}

impl<P: MNT6Parameters> ToBytesGadget<P::Fp> for G2PreparedVar<P> {
    #[inline]
    fn to_bytes(&self) -> Result<Vec<UInt8<P::Fp>>, SynthesisError> {
        let mut x = self.x.to_bytes()?;
        let mut y = self.y.to_bytes()?;
        let mut x_over_twist = self.x_over_twist.to_bytes()?;
        let mut y_over_twist = self.y_over_twist.to_bytes()?;

        x.append(&mut y);
        x.append(&mut x_over_twist);
        x.append(&mut y_over_twist);

        for coeff in self.double_coefficients.iter() {
            x.extend_from_slice(&coeff.to_bytes()?);
        }
        for coeff in self.addition_coefficients.iter() {
            x.extend_from_slice(&coeff.to_bytes()?);
        }
        Ok(x)
    }

    fn to_non_unique_bytes(&self) -> Result<Vec<UInt8<P::Fp>>, SynthesisError> {
        let mut x = self.x.to_non_unique_bytes()?;
        let mut y = self.y.to_non_unique_bytes()?;
        let mut x_over_twist = self.x_over_twist.to_non_unique_bytes()?;
        let mut y_over_twist = self.y_over_twist.to_non_unique_bytes()?;

        x.append(&mut y);
        x.append(&mut x_over_twist);
        x.append(&mut y_over_twist);

        for coeff in self.double_coefficients.iter() {
            x.extend_from_slice(&coeff.to_non_unique_bytes()?);
        }
        for coeff in self.addition_coefficients.iter() {
            x.extend_from_slice(&coeff.to_non_unique_bytes()?);
        }
        Ok(x)
    }
}

impl<P: MNT6Parameters> G2PreparedVar<P> {
    pub fn value(&self) -> Result<G2Prepared<P>, SynthesisError> {
        let x = self.x.value()?;
        let y = self.y.value()?;
        let x_over_twist = self.x_over_twist.value()?;
        let y_over_twist = self.y_over_twist.value()?;
        let double_coefficients = self
            .double_coefficients
            .iter()
            .map(|coeff| coeff.value())
            .collect::<Result<Vec<_>, SynthesisError>>()?;
        let addition_coefficients = self
            .addition_coefficients
            .iter()
            .map(|coeff| coeff.value())
            .collect::<Result<Vec<_>, SynthesisError>>()?;
        Ok(G2Prepared {
            x,
            y,
            x_over_twist,
            y_over_twist,
            double_coefficients,
            addition_coefficients,
        })
    }

    pub fn from_group_var(q: &G2Var<P>) -> Result<Self, SynthesisError> {
        let q = q.to_affine()?;
        let twist_inv = P::TWIST.inverse().unwrap();

        let mut g2p = G2PreparedVar {
            x: q.x.clone(),
            y: q.y.clone(),
            x_over_twist: &q.x * twist_inv,
            y_over_twist: &q.y * twist_inv,
            double_coefficients: vec![],
            addition_coefficients: vec![],
        };

        let mut r = G2ProjectiveExtendedVar {
            x: q.x.clone(),
            y: q.y.clone(),
            z: Fp3G::<P>::one(),
            t: Fp3G::<P>::one(),
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

            for bit in v.iter().rev() {
                let (r2, coeff) = PairingVar::<P>::doubling_step_for_flipped_miller_loop(&r)?;
                g2p.double_coefficients.push(coeff);
                r = r2;

                if *bit {
                    let (r2, coeff) = PairingVar::<P>::mixed_addition_step_for_flipped_miller_loop(
                        &q.x, &q.y, &r,
                    )?;
                    g2p.addition_coefficients.push(coeff);
                    r = r2;
                }

                tmp >>= 1;
            }
        }

        if P::ATE_IS_LOOP_COUNT_NEG {
            let rz_inv = r.z.inverse()?;
            let rz2_inv = rz_inv.square()?;
            let rz3_inv = &rz_inv * &rz2_inv;

            let minus_r_affine_x = &r.x * &rz2_inv;
            let minus_r_affine_y = r.y.negate()? * &rz3_inv;

            let add_result = PairingVar::<P>::mixed_addition_step_for_flipped_miller_loop(
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
#[derivative(Clone(bound = "P: MNT6Parameters"), Debug(bound = "P: MNT6Parameters"))]
pub struct AteDoubleCoefficientsVar<P: MNT6Parameters> {
    pub c_h: Fp3Var<P::Fp3Params>,
    pub c_4c: Fp3Var<P::Fp3Params>,
    pub c_j: Fp3Var<P::Fp3Params>,
    pub c_l: Fp3Var<P::Fp3Params>,
}

impl<P: MNT6Parameters> AllocVar<AteDoubleCoefficients<P>, P::Fp> for AteDoubleCoefficientsVar<P> {
    fn new_variable<T: Borrow<AteDoubleCoefficients<P>>>(
        cs: impl Into<Namespace<P::Fp>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let c_prep = f().map(|c| c.borrow().clone());
        let c = c_prep.as_ref().map_err(|e| *e);

        let c_h = Fp3Var::new_variable(cs.ns("c_h"), || c.map(|c| c.c_h), mode)?;
        let c_4c = Fp3Var::new_variable(cs.ns("c_4c"), || c.map(|c| c.c_4c), mode)?;
        let c_j = Fp3Var::new_variable(cs.ns("c_j"), || c.map(|c| c.c_j), mode)?;
        let c_l = Fp3Var::new_variable(cs.ns("c_l"), || c.map(|c| c.c_l), mode)?;
        Ok(Self {
            c_h,
            c_4c,
            c_j,
            c_l,
        })
    }
}

impl<P: MNT6Parameters> ToBytesGadget<P::Fp> for AteDoubleCoefficientsVar<P> {
    #[inline]
    fn to_bytes(&self) -> Result<Vec<UInt8<P::Fp>>, SynthesisError> {
        let mut c_h = self.c_h.to_bytes()?;
        let mut c_4c = self.c_4c.to_bytes()?;
        let mut c_j = self.c_j.to_bytes()?;
        let mut c_l = self.c_l.to_bytes()?;

        c_h.append(&mut c_4c);
        c_h.append(&mut c_j);
        c_h.append(&mut c_l);
        Ok(c_h)
    }

    fn to_non_unique_bytes(&self) -> Result<Vec<UInt8<P::Fp>>, SynthesisError> {
        let mut c_h = self.c_h.to_non_unique_bytes()?;
        let mut c_4c = self.c_4c.to_non_unique_bytes()?;
        let mut c_j = self.c_j.to_non_unique_bytes()?;
        let mut c_l = self.c_l.to_non_unique_bytes()?;

        c_h.append(&mut c_4c);
        c_h.append(&mut c_j);
        c_h.append(&mut c_l);
        Ok(c_h)
    }
}

impl<P: MNT6Parameters> AteDoubleCoefficientsVar<P> {
    pub fn value(&self) -> Result<AteDoubleCoefficients<P>, SynthesisError> {
        let c_h = self.c_h.value()?;
        let c_4c = self.c_4c.value()?;
        let c_j = self.c_j.value()?;
        let c_l = self.c_l.value()?;
        Ok(AteDoubleCoefficients {
            c_h,
            c_4c,
            c_j,
            c_l,
        })
    }
}

#[derive(Derivative)]
#[derivative(Clone(bound = "P: MNT6Parameters"), Debug(bound = "P: MNT6Parameters"))]
pub struct AteAdditionCoefficientsVar<P: MNT6Parameters> {
    pub c_l1: Fp3Var<P::Fp3Params>,
    pub c_rz: Fp3Var<P::Fp3Params>,
}

impl<P: MNT6Parameters> AllocVar<AteAdditionCoefficients<P>, P::Fp>
    for AteAdditionCoefficientsVar<P>
{
    fn new_variable<T: Borrow<AteAdditionCoefficients<P>>>(
        cs: impl Into<Namespace<P::Fp>>,
        f: impl FnOnce() -> Result<T, SynthesisError>,
        mode: AllocationMode,
    ) -> Result<Self, SynthesisError> {
        let ns = cs.into();
        let cs = ns.cs();

        let c_prep = f().map(|c| c.borrow().clone());
        let c = c_prep.as_ref().map_err(|e| *e);

        let c_l1 = Fp3Var::new_variable(cs.ns("c_l1"), || c.map(|c| c.c_l1), mode)?;
        let c_rz = Fp3Var::new_variable(cs.ns("c_rz"), || c.map(|c| c.c_rz), mode)?;
        Ok(Self { c_l1, c_rz })
    }
}

impl<P: MNT6Parameters> ToBytesGadget<P::Fp> for AteAdditionCoefficientsVar<P> {
    #[inline]
    fn to_bytes(&self) -> Result<Vec<UInt8<P::Fp>>, SynthesisError> {
        let mut c_l1 = self.c_l1.to_bytes()?;
        let mut c_rz = self.c_rz.to_bytes()?;

        c_l1.append(&mut c_rz);
        Ok(c_l1)
    }

    fn to_non_unique_bytes(&self) -> Result<Vec<UInt8<P::Fp>>, SynthesisError> {
        let mut c_l1 = self.c_l1.to_non_unique_bytes()?;
        let mut c_rz = self.c_rz.to_non_unique_bytes()?;

        c_l1.append(&mut c_rz);
        Ok(c_l1)
    }
}

impl<P: MNT6Parameters> AteAdditionCoefficientsVar<P> {
    pub fn value(&self) -> Result<AteAdditionCoefficients<P>, SynthesisError> {
        let c_l1 = self.c_l1.value()?;
        let c_rz = self.c_rz.value()?;
        Ok(AteAdditionCoefficients { c_l1, c_rz })
    }
}

pub struct G2ProjectiveExtendedVar<P: MNT6Parameters> {
    pub x: Fp3Var<P::Fp3Params>,
    pub y: Fp3Var<P::Fp3Params>,
    pub z: Fp3Var<P::Fp3Params>,
    pub t: Fp3Var<P::Fp3Params>,
}
