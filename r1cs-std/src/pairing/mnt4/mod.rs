use r1cs_core::{ConstraintSystem, SynthesisError};

use super::PairingGadget as PG;

use crate::{
    fields::{fp::FpGadget, fp2::Fp2Gadget, fp4::Fp4Gadget, FieldGadget},
    groups::mnt4::{
        AteAdditionCoefficientsGadget, AteDoubleCoefficientsGadget, G1Gadget, G1PreparedGadget,
        G2Gadget, G2PreparedGadget, G2ProjectiveExtendedGadget,
    },
};
use algebra::{
    curves::mnt4::{MNT4Parameters, MNT4},
    fields::BitIterator,
};
use core::marker::PhantomData;

pub struct PairingGadget<P: MNT4Parameters>(PhantomData<P>);

type Fp2G<P> = Fp2Gadget<<P as MNT4Parameters>::Fp2Params, <P as MNT4Parameters>::Fp>;
type Fp4G<P> = Fp4Gadget<<P as MNT4Parameters>::Fp4Params, <P as MNT4Parameters>::Fp>;
pub type GTGadget<P> = Fp4G<P>;

impl<P: MNT4Parameters> PairingGadget<P> {
    pub(crate) fn doubling_step_for_flipped_miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        r: &G2ProjectiveExtendedGadget<P>,
    ) -> Result<
        (
            G2ProjectiveExtendedGadget<P>,
            AteDoubleCoefficientsGadget<P>,
        ),
        SynthesisError,
    > {
        let a = r.t.square(cs.ns(|| "r.t^2"))?;
        let b = r.x.square(cs.ns(|| "r.x^2"))?;
        let c = r.y.square(cs.ns(|| "r.y^2"))?;
        let d = c.square(cs.ns(|| "c^2"))?;
        let mut e = r.x.add(cs.ns(|| "r.x + c"), &c)?;
        e.square_in_place(cs.ns(|| "(r.x + c)^2"))?;
        e.sub_in_place(cs.ns(|| "(r.x + c)^2 - b"), &b)?;
        e.sub_in_place(cs.ns(|| "(r.x + c)^2 - b - d"), &d)?;

        let mut f = b.double(cs.ns(|| "b + b"))?;
        f.add_in_place(cs.ns(|| "b + b + b"), &b)?;
        let twist_a = a.mul_by_constant(cs.ns(|| "TWIST_COEFF_A * a"), &P::TWIST_COEFF_A)?;
        f.add_in_place(cs.ns(|| "(b + b + b) + (TWIST_COEFF_A * a)"), &twist_a)?;
        let g = f.square(cs.ns(|| "f^2"))?;

        let d_eight = d
            .double(cs.ns(|| "2 * d"))?
            .double(cs.ns(|| "4 * d"))?
            .double(cs.ns(|| "8 * d"))?;

        let e2 = e.double(cs.ns(|| "2 * e"))?;
        let e4 = e2.double(cs.ns(|| "4 * e"))?;
        let x = g.sub(cs.ns(|| "- (e + e + e + e) + g"), &e4)?;

        let mut y = e2.sub(cs.ns(|| "e + e - x"), &x)?;
        y.mul_in_place(cs.ns(|| "f * (e + e - x)"), &f)?;
        y.sub_in_place(cs.ns(|| "- d_eight + f * (e + e - x)"), &d_eight)?;
        let mut z = r.y.add(cs.ns(|| "r.y + r.z"), &r.z)?;
        z.square_in_place(cs.ns(|| "(r.y + r.z)^2"))?;
        z.sub_in_place(cs.ns(|| "(r.y + r.z)^2 - c"), &c)?;
        let z2 = r.z.square(cs.ns(|| "r.z^2"))?;
        z.sub_in_place(cs.ns(|| "(r.y + r.z)^2 - c - r.z^2"), &z2)?;
        let t = z.square(cs.ns(|| "z^2"))?;

        let r2 = G2ProjectiveExtendedGadget { x, y, z, t };

        let c_h =
            r2.z.add(cs.ns(|| "r2.z + r.t"), &r.t)?
                .square(cs.ns(|| "(r2.z + r.t)^2"))?
                .sub(cs.ns(|| "(r2.z + r.t)^2 - r2.t"), &r2.t)?
                .sub(cs.ns(|| "(r2.z + r.t)^2 - r2.t - a"), &a)?;
        let c_4c = c.double(cs.ns(|| "2 * c"))?.double(cs.ns(|| "4 * c"))?;
        let mut c_j = f.add(cs.ns(|| "f + r.t"), &r.t)?;
        c_j.square_in_place(cs.ns(|| "(f + r.t)^2"))?;
        c_j.sub_in_place(cs.ns(|| "(f + r.t)^2 - g"), &g)?;
        c_j.sub_in_place(cs.ns(|| "(f + r.t)^2 - g - a"), &a)?;
        let mut c_l = f.add(cs.ns(|| "f + r.x"), &r.x)?;
        c_l.square_in_place(cs.ns(|| "(f + r.x)^2"))?;
        c_l.sub_in_place(cs.ns(|| "(f + r.x)^2 - g"), &g)?;
        c_l.sub_in_place(cs.ns(|| "(f + r.x)^2 - g - b"), &b)?;
        let coeff = AteDoubleCoefficientsGadget {
            c_h,
            c_4c,
            c_j,
            c_l,
        };

        Ok((r2, coeff))
    }

    pub(crate) fn mixed_addition_step_for_flipped_miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        x: &Fp2G<P>,
        y: &Fp2G<P>,
        r: &G2ProjectiveExtendedGadget<P>,
    ) -> Result<
        (
            G2ProjectiveExtendedGadget<P>,
            AteAdditionCoefficientsGadget<P>,
        ),
        SynthesisError,
    > {
        let a = y.square(cs.ns(|| "y^2"))?;
        let b = r.t.mul(cs.ns(|| "r.t * x"), &x)?;
        let mut d = r.z.add(cs.ns(|| "r.z + y"), &y)?;
        d.square_in_place(cs.ns(|| "(r.z + y)^2"))?;
        d.sub_in_place(cs.ns(|| "(r.z + y)^2 - a"), &a)?;
        d.sub_in_place(cs.ns(|| "(r.z + y)^2 - a - r.t"), &r.t)?;
        d.mul_in_place(cs.ns(|| "((r.z + y)^2 - a - r.t) * r.t"), &r.t)?;
        let h = b.sub(cs.ns(|| "b - r.x"), &r.x)?;
        let i = h.square(cs.ns(|| "h^2"))?;
        let e = i.double(cs.ns(|| "2 * i"))?.double(cs.ns(|| "4 * i"))?;
        let j = h.mul(cs.ns(|| "h * e"), &e)?;
        let v = r.x.mul(cs.ns(|| "r.x * e"), &e)?;
        let ry2 = r.y.double(cs.ns(|| "r.y + r.y"))?;
        let l1 = d.sub(cs.ns(|| "d - (r.y + r.y)"), &ry2)?;

        let v2 = v.double(cs.ns(|| "v + v"))?;
        let x = l1
            .square(cs.ns(|| "l1^2"))?
            .sub(cs.ns(|| "l1^2 - j"), &j)?
            .sub(cs.ns(|| "l1^2 - j - (v + v)"), &v2)?;
        let v_minus_x = v.sub(cs.ns(|| "v - x"), &x)?;
        let j_ry2 = j.mul(cs.ns(|| "j * (r.y + r.y)"), &ry2)?;
        let y = l1
            .mul(cs.ns(|| "l1 * (v - x)"), &v_minus_x)?
            .sub(cs.ns(|| "l1 * (v - x) - (j * (r.y + r.y)"), &j_ry2)?;
        let mut z = r.z.add(cs.ns(|| "r.z + h"), &h)?;
        z.square_in_place(cs.ns(|| "(r.z + h)^2"))?;
        z.sub_in_place(cs.ns(|| "(r.z + h)^2 - r.t"), &r.t)?;
        z.sub_in_place(cs.ns(|| "(r.z + h)^2 - r.t - i"), &i)?;
        let t = z.square(cs.ns(|| "z^2"))?;

        let r2 = G2ProjectiveExtendedGadget {
            x,
            y,
            z: z.clone(),
            t,
        };
        let coeff = AteAdditionCoefficientsGadget { c_l1: l1, c_rz: z };

        Ok((r2, coeff))
    }

    pub fn ate_miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        p: &G1PreparedGadget<P>,
        q: &G2PreparedGadget<P>,
    ) -> Result<Fp4G<P>, SynthesisError> {
        let mut l1_coeff = Fp2G::<P>::new(p.x.clone(), FpGadget::<P::Fp>::zero(cs.ns(|| "zero"))?);
        l1_coeff.sub_in_place(cs.ns(|| "l1_coeff"), &q.x_over_twist)?;

        let mut f = Fp4G::<P>::one(cs.ns(|| "one"))?;

        let mut dbl_idx: usize = 0;
        let mut add_idx: usize = 0;

        let mut found_one = false;

        for (j, bit) in BitIterator::new(P::ATE_LOOP_COUNT).enumerate() {
            // code below gets executed for all bits (EXCEPT the MSB itself) of
            // mnt6_param_p (skipping leading zeros) in MSB to LSB order
            if !found_one && bit {
                found_one = true;
                continue;
            } else if !found_one {
                continue;
            }

            let mut cs = cs.ns(|| format!("bit {}", j));

            let dc = &q.double_coefficients[dbl_idx];
            dbl_idx += 1;

            let c_j_x_twist = dc.c_j.mul(cs.ns(|| "dc.c_j * p.x_twist"), &p.x_twist)?;
            let c0 = dc.c_l.sub(cs.ns(|| "-dc.c_4c + dc.c_l"), &dc.c_4c)?.sub(
                cs.ns(|| "-dc.c_4c - (dc.c_j * p.x_twist) + dc.c_l"),
                &c_j_x_twist,
            )?;
            let c1 = dc.c_h.mul(cs.ns(|| "dc.c_h * p.y_twist"), &p.y_twist)?;
            let g_rr_at_p = Fp4G::<P>::new(c0, c1);

            f = f
                .square(cs.ns(|| "f^2"))?
                .mul(cs.ns(|| "f^2 * g_rr_at_p"), &g_rr_at_p)?;

            if bit {
                let ac = &q.addition_coefficients[add_idx];
                add_idx += 1;

                let l1_coeff_c_l1 = l1_coeff.mul(cs.ns(|| "l1_coeff * ac.c_l1"), &ac.c_l1)?;
                let g_rq_at_p = Fp4G::<P>::new(
                    ac.c_rz.mul(cs.ns(|| "ac.c_rz * p.y_twist"), &p.y_twist)?,
                    q.y_over_twist
                        .mul(cs.ns(|| "q.y_over_twist * ac.c_rz"), &ac.c_rz)?
                        .add(
                            cs.ns(|| "q.y_over_twist * ac.c_rz + (l1_coeff * ac.c_l1)"),
                            &l1_coeff_c_l1,
                        )?
                        .negate(cs.ns(|| "-(q.y_over_twist * ac.c_rz + (l1_coeff * ac.c_l1))"))?,
                );
                f.mul_in_place(cs.ns(|| "f *= g_rq_at_p"), &g_rq_at_p)?;
            }
        }

        if P::ATE_IS_LOOP_COUNT_NEG {
            let ac = &q.addition_coefficients[add_idx];

            let l1_coeff_c_l1 = l1_coeff.mul(cs.ns(|| "l1_coeff * ac.c_l1"), &ac.c_l1)?;
            let g_rnegr_at_p = Fp4G::<P>::new(
                ac.c_rz.mul(cs.ns(|| "ac.c_rz * p.y_twist"), &p.y_twist)?,
                q.y_over_twist
                    .mul(cs.ns(|| "q.y_over_twist * ac.c_rz"), &ac.c_rz)?
                    .add(
                        cs.ns(|| "q.y_over_twist * ac.c_rz + (l1_coeff * ac.c_l1)"),
                        &l1_coeff_c_l1,
                    )?
                    .negate(cs.ns(|| "-(q.y_over_twist * ac.c_rz + (l1_coeff * ac.c_l1))"))?,
            );
            f = f
                .mul(cs.ns(|| "f * g_rnegr_at_p"), &g_rnegr_at_p)?
                .inverse(cs.ns(|| "inverse f"))?;
        }

        Ok(f)
    }

    pub fn final_exponentiation<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        value: &Fp4G<P>,
    ) -> Result<GTGadget<P>, SynthesisError> {
        let value_inv = value.inverse(cs.ns(|| "value inverse"))?;
        let value_to_first_chunk = Self::final_exponentiation_first_chunk(
            cs.ns(|| "value_to_first_chunk"),
            value,
            &value_inv,
        )?;
        let value_inv_to_first_chunk = Self::final_exponentiation_first_chunk(
            cs.ns(|| "value_inv_to_first_chunk"),
            &value_inv,
            value,
        )?;
        Self::final_exponentiation_last_chunk(
            cs.ns(|| "final_exp_last_chunk"),
            &value_to_first_chunk,
            &value_inv_to_first_chunk,
        )
    }

    fn final_exponentiation_first_chunk<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        elt: &Fp4G<P>,
        elt_inv: &Fp4G<P>,
    ) -> Result<Fp4G<P>, SynthesisError> {
        // (q^2-1)

        // elt_q2 = elt^(q^2)
        let mut elt_q2 = elt.clone();
        elt_q2.frobenius_map_in_place(cs.ns(|| "frobenius 2"), 2)?;
        // elt_q2_over_elt = elt^(q^2-1)
        elt_q2.mul(cs.ns(|| "elt_q2 * elt_inv"), elt_inv)
    }

    fn final_exponentiation_last_chunk<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        elt: &Fp4G<P>,
        elt_inv: &Fp4G<P>,
    ) -> Result<Fp4G<P>, SynthesisError> {
        let elt_clone = elt.clone();
        let elt_inv_clone = elt_inv.clone();

        let mut elt_q = elt.clone();
        elt_q.frobenius_map_in_place(cs.ns(|| "frobenius 1"), 1)?;

        let w1_part = elt_q.cyclotomic_exp(cs.ns(|| "w1_part"), &P::FINAL_EXPONENT_LAST_CHUNK_1)?;
        let w0_part;
        if P::FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG {
            w0_part = elt_inv_clone
                .cyclotomic_exp(cs.ns(|| "w0_part"), &P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0)?;
        } else {
            w0_part = elt_clone
                .cyclotomic_exp(cs.ns(|| "w0_part"), &P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0)?;
        }

        w1_part.mul(cs.ns(|| "w1_part * w0_part"), &w0_part)
    }
}

impl<P: MNT4Parameters> PG<MNT4<P>, P::Fp> for PairingGadget<P> {
    type G1Gadget = G1Gadget<P>;
    type G2Gadget = G2Gadget<P>;
    type G1PreparedGadget = G1PreparedGadget<P>;
    type G2PreparedGadget = G2PreparedGadget<P>;
    type GTGadget = GTGadget<P>;

    fn miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        ps: &[Self::G1PreparedGadget],
        qs: &[Self::G2PreparedGadget],
    ) -> Result<Self::GTGadget, SynthesisError> {
        let mut result = Fp4G::<P>::one(cs.ns(|| "one"))?;
        for (i, (p, q)) in ps.iter().zip(qs.iter()).enumerate() {
            let miller =
                Self::ate_miller_loop(cs.ns(|| format!("ate miller loop iteration {}", i)), p, q)?;
            result.mul_in_place(
                cs.ns(|| format!("mul ate miller loop iteration {}", i)),
                &miller,
            )?;
        }

        Ok(result)
    }

    fn final_exponentiation<CS: ConstraintSystem<P::Fp>>(
        cs: CS,
        r: &Self::GTGadget,
    ) -> Result<Self::GTGadget, SynthesisError> {
        Self::final_exponentiation(cs, r)
    }

    fn prepare_g1<CS: ConstraintSystem<P::Fp>>(
        cs: CS,
        p: &Self::G1Gadget,
    ) -> Result<Self::G1PreparedGadget, SynthesisError> {
        Self::G1PreparedGadget::from_affine(cs, p)
    }

    fn prepare_g2<CS: ConstraintSystem<P::Fp>>(
        cs: CS,
        q: &Self::G2Gadget,
    ) -> Result<Self::G2PreparedGadget, SynthesisError> {
        Self::G2PreparedGadget::from_affine(cs, q)
    }
}
