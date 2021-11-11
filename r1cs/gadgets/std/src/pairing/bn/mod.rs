use r1cs_core::{ConstraintSystem, SynthesisError};

use super::PairingGadget as PG;

use crate::{
    fields::{
        fp::FpGadget, fp12::Fp12Gadget, fp2::Fp2Gadget,
        quadratic_extension::QuadExtParametersGadget, FieldGadget,
    },
    groups::bn::{G1Gadget, G1PreparedGadget, G2Gadget, G2PreparedGadget},
};
use algebra::{curves::bn::*, fields::fp12_2over3over2::Fp12ParamsWrapper};
use std::marker::PhantomData;

pub struct PairingGadget<P: BnParameters>(PhantomData<P>);

type Fp2G<P> = Fp2Gadget<<P as BnParameters>::Fp2Params, <P as BnParameters>::Fp>;

impl<P: BnParameters> PairingGadget<P> {
    // Evaluate the line function at point p.
    fn ell<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        f: &mut Fp12Gadget<P::Fp12Params, P::Fp>,
        coeffs: &(Fp2G<P>, Fp2G<P>),
        p: &G1Gadget<P>,
    ) -> Result<(), SynthesisError> {
        let zero = FpGadget::<P::Fp>::zero(cs.ns(|| "fpg zero"))?;

        match P::TWIST_TYPE {
            TwistType::M => {
                let c0 = coeffs.0.clone();
                let mut c1 = coeffs.1.clone();
                let c2 = Fp2G::<P>::new(p.y.clone(), zero);

                c1.c0 = c1.c0.mul(cs.ns(|| "mul c1.c0"), &p.x)?;
                c1.c1 = c1.c1.mul(cs.ns(|| "mul c1.c1"), &p.x)?;
                *f = f.mul_by_014(cs.ns(|| "sparse mul f"), &c0, &c1, &c2)?;
                Ok(())
            }
            TwistType::D => {
                let c0 = Fp2G::<P>::new(p.y.clone(), zero);
                let mut c1 = coeffs.0.clone();
                let c2 = coeffs.1.clone();

                c1.c0 = c1.c0.mul(cs.ns(|| "mul c1.c0"), &p.x)?;
                c1.c1 = c1.c1.mul(cs.ns(|| "mul c1.c1"), &p.x)?;
                *f = f.mul_by_034(cs.ns(|| "sparse mul f"), &c0, &c1, &c2)?;
                Ok(())
            }
        }
    }

    fn exp_by_neg_x<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        f: &Fp12Gadget<P::Fp12Params, P::Fp>,
    ) -> Result<Fp12Gadget<P::Fp12Params, P::Fp>, SynthesisError> {
        let mut result = f.cyclotomic_exp(cs.ns(|| "exp_by_neg_x"), P::X)?;
        if !P::X_IS_NEGATIVE {
            result.conjugate_in_place(cs.ns(|| "conjugate"))?;
        }
        Ok(result)
    }
}

impl<P: BnParameters> PG<Bn<P>, P::Fp> for PairingGadget<P> {
    type G1Gadget = G1Gadget<P>;
    type G2Gadget = G2Gadget<P>;
    type G1PreparedGadget = G1PreparedGadget<P>;
    type G2PreparedGadget = G2PreparedGadget<P>;
    type GTGadget = Fp12Gadget<P::Fp12Params, P::Fp>;

    fn miller_loop<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        ps: &[Self::G1PreparedGadget],
        qs: &[Self::G2PreparedGadget],
    ) -> Result<Self::GTGadget, SynthesisError> {
        let mut pairs = vec![];
        for (p, q) in ps.iter().zip(qs.iter()) {
            pairs.push((p, q.ell_coeffs.iter()));
        }
        let mut f = Self::GTGadget::one(cs.ns(|| "one"))?;

        for i in (1..P::ATE_LOOP_COUNT.len()).rev() {
            let mut cs = cs.ns(|| format!("Iteration {}", i));
            if i != P::ATE_LOOP_COUNT.len() - 1 {
                f.square_in_place(cs.ns(|| "square"))?;
            }

            for (k, &mut (p, ref mut coeffs)) in pairs.iter_mut().enumerate() {
                let cs = cs.ns(|| format!("Double input {}", k));
                Self::ell(cs, &mut f, coeffs.next().unwrap(), &p.0)?;
            }

            let bit = P::ATE_LOOP_COUNT[i - 1];

            match bit {
                1 => {
                    for (k, &mut (p, ref mut coeffs)) in pairs.iter_mut().enumerate() {
                        let cs = cs.ns(|| format!("Addition input {}", k));
                        Self::ell(cs, &mut f, &coeffs.next().unwrap(), &p.0)?;
                    }
                }
                -1 => {
                    for (k, &mut (p, ref mut coeffs)) in pairs.iter_mut().enumerate() {
                        let cs = cs.ns(|| format!("Addition input {}", k));
                        Self::ell(cs, &mut f, &coeffs.next().unwrap(), &p.0)?;
                    }
                }
                _ => continue,
            }
        }

        if P::ATE_LOOP_COUNT_IS_NEGATIVE {
            f.conjugate_in_place(cs.ns(|| "f conjugate"))?;
        }

        for (i, &mut (p, ref mut coeffs)) in pairs.iter_mut().enumerate() {
            Self::ell(
                cs.ns(|| format!("Last addition step 1_{}", i)),
                &mut f,
                coeffs.next().unwrap(),
                &p.0,
            )?;
        }

        for (i, &mut (p, ref mut coeffs)) in pairs.iter_mut().enumerate() {
            Self::ell(
                cs.ns(|| format!("Last addition step 2_{}", i)),
                &mut f,
                coeffs.next().unwrap(),
                &p.0,
            )?;
        }

        Ok(f)
    }

    fn final_exponentiation<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        f: &Self::GTGadget,
    ) -> Result<Self::GTGadget, SynthesisError> {
        // Computing the final exponentation following
        // https://eprint.iacr.org/2016/130.pdf.
        // We don't use their "faster" formula because it is difficult to make
        // it work for curves with odd `P::X`.
        // Hence we implement the slower algorithm from Table 1 below.

        let f1 = f.frobenius_map(cs.ns(|| "frobmap 1"), 6)?;

        f.inverse(cs.ns(|| "inverse 1")).and_then(|mut f2| {
            // f2 = f^(-1);
            // r = f^(p^6 - 1)
            let mut r = f1;
            r.mul_in_place(cs.ns(|| "r = f1 * f2"), &f2)?;

            // f2 = f^(p^6 - 1)
            f2 = r.clone();
            // r = f^((p^6 - 1)(p^2))
            r.frobenius_map_in_place(cs.ns(|| "frobmap 2"), 2)?;

            // r = f^((p^6 - 1)(p^2) + (p^6 - 1))
            // r = f^((p^6 - 1)(p^2 + 1))
            r.mul_in_place(cs.ns(|| "mul 0"), &f2)?;

            // Hard part follows Laura Fuentes-Castaneda et al. "Faster hashing to G2"
            // by computing:
            //
            // result = elt^(q^3 * (12*z^3 + 6z^2 + 4z - 1) +
            //               q^2 * (12*z^3 + 6z^2 + 6z) +
            //               q   * (12*z^3 + 6z^2 + 4z) +
            //               1   * (12*z^3 + 12z^2 + 6z + 1))
            // which equals
            //
            // result = elt^( 2z * ( 6z^2 + 3z + 1 ) * (q^4 - q^2 + 1)/r ).

            let y0 = Self::exp_by_neg_x(cs.ns(|| "exp_by_neg_x_1"), &r)?;
            let y1 = Fp12ParamsWrapper::<P::Fp12Params>::cyclotomic_square_gadget(
                cs.ns(|| "square_1"),
                &y0,
            )?;
            let y2 = Fp12ParamsWrapper::<P::Fp12Params>::cyclotomic_square_gadget(
                cs.ns(|| "square_2"),
                &y1,
            )?;
            let mut y3 = y2.mul(cs.ns(|| "y3 = y2 * y1"), &y1)?;
            let y4 = Self::exp_by_neg_x(cs.ns(|| "exp_by_neg_x_2"), &y3)?;
            let y5 = Fp12ParamsWrapper::<P::Fp12Params>::cyclotomic_square_gadget(
                cs.ns(|| "square_3"),
                &y4,
            )?;
            let mut y6 = Self::exp_by_neg_x(cs.ns(|| "exp_by_neg_x_3"), &y5)?;
            y3.conjugate_in_place(cs.ns(|| "conjugate 1"))?;
            y6.conjugate_in_place(cs.ns(|| "conjugate_2"))?;
            let y7 = y6.mul(cs.ns(|| "y7 = y6 * y4"), &y4)?;
            let mut y8 = y7.mul(cs.ns(|| "y8 = y7 * y3"), &y3)?;
            let y9 = y8.mul(cs.ns(|| "y9 = y8 * y1"), &y1)?;
            let y10 = y8.mul(cs.ns(|| "y10 = y8 * y4"), &y4)?;
            let y11 = y10.mul(cs.ns(|| "y11 = y10 * r"), &r)?;
            let mut y12 = y9.clone();
            y12.frobenius_map_in_place(cs.ns(|| "frobmap 3"), 1)?;
            let y13 = y12.mul(cs.ns(|| "y13 = y12 * y11"), &y11)?;
            y8.frobenius_map_in_place(cs.ns(|| "frobmap 4"), 2)?;
            let y14 = y8.mul(cs.ns(|| "y14 = y8 * y13"), &y13)?;
            r.conjugate_in_place(cs.ns(|| "conjugate_3"))?;
            let mut y15 = r.mul(cs.ns(|| "y15 = r * y9"), &y9)?;
            y15.frobenius_map_in_place(cs.ns(|| "frobmap 5"), 3)?;
            let y16 = y15.mul(cs.ns(|| "y16 = y15 * y14"), &y14)?;

            Ok(y16)
        })
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
