use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{
    fields::{fp6_2over3::Fp6Gadget, FieldGadget},
    groups::curves::short_weierstrass::mnt::mnt6::{
        G1Gadget, G1PreparedGadget, G2Gadget, G2PreparedGadget,
    },
};

use crate::pairing::PairingGadget;
use algebra::curves::models::mnt6::{MNT6Parameters, MNT6p};
use std::marker::PhantomData;

pub struct MNT6PairingGadget<P: MNT6Parameters>(PhantomData<P>);

impl<P: MNT6Parameters> PairingGadget<MNT6p<P>, P::Fp> for MNT6PairingGadget<P> {
    type G1Gadget = G1Gadget<P>;
    type G2Gadget = G2Gadget<P>;
    type G1PreparedGadget = G1PreparedGadget<P>;
    type G2PreparedGadget = G2PreparedGadget<P>;
    type GTGadget = Fp6Gadget<P::Fp6Params, P::Fp>;

    fn miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        p: &[Self::G1PreparedGadget],
        q: &[Self::G2PreparedGadget],
    ) -> Result<Self::GTGadget, SynthesisError> {
        let mut result = Self::GTGadget::one(cs.ns(|| "one"))?;
        let it = p.iter().zip(q.iter());

        for (i, (ps, qs)) in it.into_iter().enumerate() {
            let mut cs = cs.ns(|| format!("Pair_{}", i));

            let mut f = Self::GTGadget::one(cs.ns(|| "f"))?;

            let mut idx: usize = 0;

            for (j, &n) in P::WNAF.iter().rev().enumerate() {
                let mut cs = cs.ns(|| format!("Iteration_{}", j));

                let c = &qs.coeffs[idx];
                idx += 1;

                //Double step
                //Compute g_rr_at_p_c0
                let g_rr_at_p_c0 = ps.clone().p_y_twist_squared;

                let mut t = c
                    .gamma
                    .mul_by_constant(cs.ns(|| "double compute gamma_twist"), &P::TWIST)?;
                t.mul_assign_by_base_field_gadget(
                    cs.ns(|| "double gamma_twist * ps.p.x"),
                    &ps.p.x,
                )?;
                let g_rr_at_p_c1 = c
                    .gamma_x
                    .sub(cs.ns(|| "gamma_x - r_y"), &c.r_y)?
                    .sub(cs.ns(|| "gamma_x - r_y - t"), &t)?;

                //Compute g_rr_at_p
                let g_rr_at_p = Self::GTGadget::new(g_rr_at_p_c0.clone(), g_rr_at_p_c1);

                //Compute new_f
                f = f
                    .square(cs.ns(|| "f^2"))?
                    .mul_by_2345(cs.ns(|| "double compute f"), &g_rr_at_p)?;

                if n != 0 {
                    //Addition Step
                    let c = &qs.coeffs[idx];
                    idx += 1;

                    let g_rq_at_p_c0 = ps.clone().p_y_twist_squared;

                    //Compute g_rq_at_p_c1
                    let neg_q_y = qs.q.y.negate(cs.ns(|| "- q.y"))?;
                    let q_y = if n > 0 { qs.clone().q.y } else { neg_q_y };

                    let mut t = c
                        .gamma
                        .mul_by_constant(cs.ns(|| "add compute gamma_twist"), &P::TWIST)?;
                    t.mul_assign_by_base_field_gadget(
                        cs.ns(|| "add gamma_twist * ps.p.x"),
                        &ps.p.x,
                    )?;
                    let g_rq_at_p_c1 = c
                        .gamma_x
                        .sub(cs.ns(|| "gamma_x - q_y"), &q_y)?
                        .sub(cs.ns(|| "gamma_x - q_y - t"), &t)?;

                    //Compute g_rq_at_p
                    let g_rq_at_p = Self::GTGadget::new(g_rq_at_p_c0, g_rq_at_p_c1);

                    //Compute new f
                    f = f.mul_by_2345(cs.ns(|| "add compute f"), &g_rq_at_p)?;
                }
            }
            if P::ATE_IS_LOOP_COUNT_NEG {
                f = f.unitary_inverse(cs.ns(|| "f unitary inverse"))?;
            }
            result.mul_in_place(cs.ns(|| format!("mul_assign_{}", i)), &f)?;
        }
        Ok(result)
    }

    fn final_exponentiation<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        value: &Self::GTGadget,
    ) -> Result<Self::GTGadget, SynthesisError> {
        let value_inv = value.inverse(cs.ns(|| "value_inverse"))?;

        //Final exp first chunk
        //use the Frobenius map a to compute value^{(q^3-1)(q-1)}
        let elt = {
            let elt_q3_over_elt = value
                .clone()
                .frobenius_map(cs.ns(|| "elt^(q^3)"), 3)?
                .mul(cs.ns(|| "elt^(q^3-1)"), &value_inv)?;
            elt_q3_over_elt
                .frobenius_map(cs.ns(|| "elt^((q^3-1) * q)"), 1)?
                .mul(cs.ns(|| "elt^((q^3-1)*(q+1)"), &elt_q3_over_elt)?
        };

        //Final exp last chunk (q^2 -q +1)/r = m_1*q + m_0, m_0 can be signed.
        //compute elt^q
        let elt_q = elt
            .clone()
            .frobenius_map(cs.ns(|| "elt_q_frobenius_1"), 1)?;

        let w1_part =
            elt_q.cyclotomic_exp(cs.ns(|| "compute w1"), P::FINAL_EXPONENT_LAST_CHUNK_1)?;

        let w0_part = if P::FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG {
            // we need the inverse of elt in this case, by recomputing first chunk exp
            let elt_inv = {
                let elt_inv_q3_over_elt_inv = value_inv
                    .frobenius_map(cs.ns(|| "elt_inv^(q^3)"), 3)?
                    .mul(cs.ns(|| "elt_inv^(q^3-1)"), &value_inv)?;
                elt_inv_q3_over_elt_inv
                    .frobenius_map(cs.ns(|| "elt_inv^((q^3-1) * q)"), 1)?
                    .mul(cs.ns(|| "elt_inv^((q^3-1)*(q+1)"), &elt_inv_q3_over_elt_inv)?
            };
            elt_inv.cyclotomic_exp(
                cs.ns(|| "compute w0"),
                P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0,
            )
        } else {
            elt.cyclotomic_exp(
                cs.ns(|| "compute w0"),
                P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0,
            )
        }?;

        w1_part.mul(cs.ns(|| "w0 * w1"), &w0_part)
    }

    fn prepare_g1<CS: ConstraintSystem<P::Fp>>(
        cs: CS,
        q: &Self::G1Gadget,
    ) -> Result<Self::G1PreparedGadget, SynthesisError> {
        Self::G1PreparedGadget::from_affine(cs, q)
    }

    fn prepare_g2<CS: ConstraintSystem<P::Fp>>(
        cs: CS,
        q: &Self::G2Gadget,
    ) -> Result<Self::G2PreparedGadget, SynthesisError> {
        Self::G2PreparedGadget::from_affine(cs, q)
    }
}
