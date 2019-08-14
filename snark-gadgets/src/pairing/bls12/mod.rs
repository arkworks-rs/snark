use r1cs_core::{ConstraintSystem, SynthesisError};

use super::PairingGadget as PG;

use crate::{
    fields::{fp::FpGadget, fp12::Fp12Gadget, fp2::Fp2Gadget, FieldGadget},
    groups::bls12::{G1Gadget, G1PreparedGadget, G2Gadget, G2PreparedGadget},
};
use algebra::{
    curves::{
        bls12::{
            Bls12, Bls12Parameters, G1Affine, G1Prepared, G1Projective, G2Affine, G2Prepared,
            G2Projective, TwistType,
        },
        models::ModelParameters,
        PairingCurve,
    },
    fields::{fp12_2over3over2::Fp12, BitIterator},
};
use std::marker::PhantomData;

pub mod bls12_377;

pub struct PairingGadget<P: Bls12Parameters>(PhantomData<P>);

type Fp2G<P> = Fp2Gadget<<P as Bls12Parameters>::Fp2Params, <P as Bls12Parameters>::Fp>;

impl<P: Bls12Parameters> PairingGadget<P> {
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
                let c2 = Fp2G::<P>::new(p.y.clone(), zero.clone());

                c1.c0 = c1.c0.mul(cs.ns(|| "mul c1.c0"), &p.x)?;
                c1.c1 = c1.c1.mul(cs.ns(|| "mul c1.c1"), &p.x)?;
                *f = f.mul_by_014(cs.ns(|| "sparse mul f"), &c0, &c1, &c2)?;
                Ok(())
            },
            TwistType::D => {
                let c0 = Fp2G::<P>::new(p.y.clone(), zero.clone());
                let mut c1 = coeffs.0.clone();
                let c2 = coeffs.1.clone();

                c1.c0 = c1.c0.mul(cs.ns(|| "mul c1.c0"), &p.x)?;
                c1.c1 = c1.c1.mul(cs.ns(|| "mul c1.c1"), &p.x)?;
                *f = f.mul_by_034(cs.ns(|| "sparse mul f"), &c0, &c1, &c2)?;
                Ok(())
            },
        }
    }

    fn exp_by_x<CS: ConstraintSystem<P::Fp>>(
        mut cs: CS,
        f: &Fp12Gadget<P::Fp12Params, P::Fp>,
    ) -> Result<Fp12Gadget<P::Fp12Params, P::Fp>, SynthesisError> {
        let mut result = f.cyclotomic_exp(cs.ns(|| "exp_by_x"), P::X)?;
        if P::X_IS_NEGATIVE {
            result.conjugate_in_place(cs.ns(|| "conjugate"))?;
        }
        Ok(result)
    }
}

impl<P: Bls12Parameters> PG<Bls12<P>, P::Fp> for PairingGadget<P>
where
    G1Affine<P>: PairingCurve<
        BaseField = <P::G1Parameters as ModelParameters>::BaseField,
        ScalarField = <P::G1Parameters as ModelParameters>::ScalarField,
        Projective = G1Projective<P>,
        PairWith = G2Affine<P>,
        Prepared = G1Prepared<P>,
        PairingResult = Fp12<P::Fp12Params>,
    >,
    G2Affine<P>: PairingCurve<
        BaseField = <P::G2Parameters as ModelParameters>::BaseField,
        ScalarField = <P::G1Parameters as ModelParameters>::ScalarField,
        Projective = G2Projective<P>,
        PairWith = G1Affine<P>,
        Prepared = G2Prepared<P>,
        PairingResult = Fp12<P::Fp12Params>,
    >,
{
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
        for (p, q) in ps.into_iter().zip(qs.into_iter()) {
            pairs.push((p, q.ell_coeffs.iter()));
        }
        let mut f = Self::GTGadget::one(cs.ns(|| "one"))?;

        for (j, i) in BitIterator::new(P::X).skip(1).enumerate() {
            let mut cs = cs.ns(|| format!("Iteration {}", j));
            f.square_in_place(cs.ns(|| "square"))?;

            for (k, &mut (p, ref mut coeffs)) in pairs.iter_mut().enumerate() {
                let cs = cs.ns(|| format!("Double input {}", k));
                Self::ell(cs, &mut f, coeffs.next().unwrap(), &p.0)?;
            }

            if i {
                for (k, &mut (p, ref mut coeffs)) in pairs.iter_mut().enumerate() {
                    let cs = cs.ns(|| format!("Addition input {}", k));
                    Self::ell(cs, &mut f, &coeffs.next().unwrap(), &p.0)?;
                }
            }
        }

        if P::X_IS_NEGATIVE {
            f.conjugate_in_place(cs.ns(|| "f conjugate"))?;
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
            r.frobenius_map_in_place(cs.ns(|| "frobenius map 2"), 2)?;

            // r = f^((p^6 - 1)(p^2) + (p^6 - 1))
            // r = f^((p^6 - 1)(p^2 + 1))
            r.mul_in_place(cs.ns(|| "mul 0"), &f2)?;

            // Hard part of the final exponentation is below:
            // From https://eprint.iacr.org/2016/130.pdf, Table 1
            let mut y0 = r.cyclotomic_square(cs.ns(|| "cyclotomic_sq 1"))?;
            y0.conjugate_in_place(&mut cs.ns(|| "conjugate 2"))?;

            let mut y5 = Self::exp_by_x(&mut cs.ns(|| "exp_by_x 1"), &r)?;

            let mut y1 = y5.cyclotomic_square(&mut cs.ns(|| "square 1"))?;
            let mut y3 = y0.mul(&mut cs.ns(|| "mul 1"), &y5)?;
            y0 = Self::exp_by_x(cs.ns(|| "exp_by_x 2"), &y3)?;
            let y2 = Self::exp_by_x(cs.ns(|| "exp_by_x 3"), &y0)?;
            let mut y4 = Self::exp_by_x(cs.ns(|| "exp_by_x 4"), &y2)?;
            y4.mul_in_place(cs.ns(|| "mul 2"), &y1)?;
            y1 = Self::exp_by_x(cs.ns(|| "exp_by_x 5"), &y4)?;
            y3.conjugate_in_place(cs.ns(|| "conjugate 3"))?;
            y1.mul_in_place(cs.ns(|| "mul 3"), &y3)?;
            y1.mul_in_place(cs.ns(|| "mul 4"), &r)?;
            y3 = r.clone();
            y3.conjugate_in_place(cs.ns(|| "conjugate 4"))?;
            y0.mul_in_place(cs.ns(|| "mul 5"), &r)?;
            y0.frobenius_map_in_place(cs.ns(|| "frobmap 3"), 3)?;
            y4.mul_in_place(cs.ns(|| "mul 6"), &y3)?;
            y4.frobenius_map_in_place(cs.ns(|| "frobmap 4"), 1)?;
            y5.mul_in_place(cs.ns(|| "mul 7"), &y2)?;
            y5.frobenius_map_in_place(cs.ns(|| "frobmap 5"), 2)?;
            y5.mul_in_place(cs.ns(|| "mul 8"), &y0)?;
            y5.mul_in_place(cs.ns(|| "mul 9"), &y4)?;
            y5.mul_in_place(cs.ns(|| "mul 10"), &y1)?;
            Ok(y5)
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
