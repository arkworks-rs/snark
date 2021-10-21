use algebra::{
    curves::bn::{BnParameters, G1Prepared, TwistType},
    fields::Field,
    ProjectiveCurve,
};
use r1cs_core::{ConstraintSystem, SynthesisError};

use crate::{
    fields::{fp::FpGadget, fp2::Fp2Gadget, FieldGadget},
    groups::curves::short_weierstrass::AffineGadget,
    prelude::*,
};
use std::fmt::Debug;

pub type G1Gadget<P> = AffineGadget<
    <P as BnParameters>::G1Parameters,
    <P as BnParameters>::Fp,
    FpGadget<<P as BnParameters>::Fp>,
>;

pub type G2Gadget<P> =
    AffineGadget<<P as BnParameters>::G2Parameters, <P as BnParameters>::Fp, Fp2G<P>>;

#[derive(Derivative)]
#[derivative(
    Clone(bound = "G1Gadget<P>: Clone"),
    Debug(bound = "G1Gadget<P>: Debug")
)]
pub struct G1PreparedGadget<P: BnParameters>(pub G1Gadget<P>);

impl<P: BnParameters> G1PreparedGadget<P> {
    pub fn get_value(&self) -> Option<G1Prepared<P>> {
        Some(G1Prepared::from(self.0.get_value().unwrap().into_affine()))
    }

    pub fn from_affine<CS: ConstraintSystem<P::Fp>>(
        _cs: CS,
        q: &G1Gadget<P>,
    ) -> Result<Self, SynthesisError> {
        Ok(G1PreparedGadget(q.clone()))
    }
}

impl<P: BnParameters> ToBytesGadget<P::Fp> for G1PreparedGadget<P> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.0.to_bytes(&mut cs.ns(|| "g_alpha to bytes"))
    }

    fn to_bytes_strict<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        self.0.to_bytes_strict(&mut cs.ns(|| "g_alpha to bytes"))
    }
}

type Fp2G<P> = Fp2Gadget<<P as BnParameters>::Fp2Params, <P as BnParameters>::Fp>;
type LCoeff<P> = (Fp2G<P>, Fp2G<P>);
#[derive(Derivative)]
#[derivative(
    Clone(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Clone"),
    Debug(bound = "Fp2Gadget<P::Fp2Params, P::Fp>: Debug")
)]
pub struct G2PreparedGadget<P: BnParameters> {
    pub ell_coeffs: Vec<LCoeff<P>>,
}

impl<P: BnParameters> ToBytesGadget<P::Fp> for G2PreparedGadget<P> {
    #[inline]
    fn to_bytes<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut bytes = Vec::new();
        for (i, coeffs) in self.ell_coeffs.iter().enumerate() {
            let mut cs = cs.ns(|| format!("Iteration {}", i));
            bytes.extend_from_slice(&coeffs.0.to_bytes(&mut cs.ns(|| "c0"))?);
            bytes.extend_from_slice(&coeffs.1.to_bytes(&mut cs.ns(|| "c1"))?);
        }
        Ok(bytes)
    }

    fn to_bytes_strict<CS: ConstraintSystem<P::Fp>>(
        &self,
        mut cs: CS,
    ) -> Result<Vec<UInt8>, SynthesisError> {
        let mut bytes = Vec::new();
        for (i, coeffs) in self.ell_coeffs.iter().enumerate() {
            let mut cs = cs.ns(|| format!("Iteration {}", i));
            bytes.extend_from_slice(&coeffs.0.to_bytes_strict(&mut cs.ns(|| "c0"))?);
            bytes.extend_from_slice(&coeffs.1.to_bytes_strict(&mut cs.ns(|| "c1"))?);
        }
        Ok(bytes)
    }
}

fn mul_by_char<P: BnParameters, CS: ConstraintSystem<P::Fp>>(
    mut cs: CS,
    q: &G2Gadget<P>,
) -> Result<G2Gadget<P>, SynthesisError> {
    let mut s = q.clone();
    s.x.frobenius_map_in_place(cs.ns(|| "s.x.frobenius_map_1"), 1)?;
    s.x.mul_by_constant_in_place(cs.ns(|| "s.x *= TWIST_MUL_BY_Q_X"), &P::TWIST_MUL_BY_Q_X)?;
    s.y.frobenius_map_in_place(cs.ns(|| "s.y.frobenius_map_1"), 1)?;
    s.y.mul_by_constant_in_place(cs.ns(|| "s.y *= TWIST_MUL_BY_Q_Y"), &P::TWIST_MUL_BY_Q_Y)?;
    Ok(s)
}

impl<P: BnParameters> G2PreparedGadget<P> {
    pub fn from_affine<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        q: &G2Gadget<P>,
    ) -> Result<Self, SynthesisError> {
        let two_inv = P::Fp::one().double().inverse().unwrap();
        let zero = G2Gadget::<P>::zero(cs.ns(|| "zero"))?;
        q.enforce_not_equal(cs.ns(|| "enforce not zero"), &zero)?;
        let mut ell_coeffs = vec![];
        let mut r = q.clone();
        let negq = q.negate(cs.ns(|| "-q"))?;

        for i in (1..P::ATE_LOOP_COUNT.len()).rev() {
            let mut cs = cs.ns(|| format!("Iteration {}", i));
            ell_coeffs.push(Self::double(cs.ns(|| "double"), &mut r, &two_inv)?);

            let bit = P::ATE_LOOP_COUNT[i - 1];

            match bit {
                1 => ell_coeffs.push(Self::add(cs.ns(|| "add_q"), &mut r, &q)?),
                -1 => ell_coeffs.push(Self::add(cs.ns(|| "add_neg_q"), &mut r, &negq)?),
                _ => continue,
            }
        }

        let q1 = mul_by_char::<P, _>(cs.ns(|| "q1 = q * char"), &q)?;
        let mut q2 = mul_by_char::<P, _>(cs.ns(|| " q2 = q1 * char"), &q1)?;

        if P::ATE_LOOP_COUNT_IS_NEGATIVE {
            r.y = r.y.negate(cs.ns(|| "-r.y"))?;
        }

        q2.y = q2.y.negate(cs.ns(|| "-q2.y"))?;

        ell_coeffs.push(Self::add(cs.ns(|| "add_last_1"), &mut r, &q1)?);
        ell_coeffs.push(Self::add(cs.ns(|| "add_last_2"), &mut r, &q2)?);

        Ok(Self { ell_coeffs })
    }

    fn double<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        r: &mut G2Gadget<P>,
        two_inv: &P::Fp,
    ) -> Result<LCoeff<P>, SynthesisError> {
        let a = r.y.inverse(cs.ns(|| "Inverse"))?;
        let mut b = r.x.square(cs.ns(|| "square x"))?;
        let b_tmp = b.clone();
        b.mul_by_base_field_constant_in_place(cs.ns(|| "mul by two_inv"), two_inv)?;
        b.add_in_place(cs.ns(|| "compute b"), &b_tmp)?;

        let c = a.mul(cs.ns(|| "compute c"), &b)?;
        let d = r.x.double(cs.ns(|| "compute d"))?;
        let x3 = c.square(cs.ns(|| "c^2"))?.sub(cs.ns(|| "sub d"), &d)?;
        let e = c
            .mul(cs.ns(|| "c*r.x"), &r.x)?
            .sub(cs.ns(|| "sub r.y"), &r.y)?;
        let c_x3 = c.mul(cs.ns(|| "c*x_3"), &x3)?;
        let y3 = e.sub(cs.ns(|| "e = c * x3"), &c_x3)?;
        let mut f = c;
        f.negate_in_place(cs.ns(|| "c = -c"))?;
        r.x = x3;
        r.y = y3;
        match P::TWIST_TYPE {
            TwistType::M => Ok((e, f)),
            TwistType::D => Ok((f, e)),
        }
    }

    fn add<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        r: &mut G2Gadget<P>,
        q: &G2Gadget<P>,
    ) -> Result<LCoeff<P>, SynthesisError> {
        let a =
            q.x.sub(cs.ns(|| "q.x - r.x"), &r.x)?
                .inverse(cs.ns(|| "calc a"))?;
        let b = q.y.sub(cs.ns(|| "q.y - r.y"), &r.y)?;
        let c = a.mul(cs.ns(|| "compute c"), &b)?;
        let d = r.x.add(cs.ns(|| "r.x + q.x"), &q.x)?;
        let x3 = c.square(cs.ns(|| "c^2"))?.sub(cs.ns(|| "sub d"), &d)?;

        let e =
            r.x.sub(cs.ns(|| "r.x - x3"), &x3)?
                .mul(cs.ns(|| "c * (r.x - x3)"), &c)?;
        let y3 = e.sub(cs.ns(|| "calc y3"), &r.y)?;
        let g = c
            .mul(cs.ns(|| "c*r.x"), &r.x)?
            .sub(cs.ns(|| "calc g"), &r.y)?;
        let mut f = c;
        f.negate_in_place(cs.ns(|| "c = -c"))?;
        r.x = x3;
        r.y = y3;
        match P::TWIST_TYPE {
            TwistType::M => Ok((g, f)),
            TwistType::D => Ok((f, g)),
        }
    }
}
