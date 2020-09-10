use r1cs_core::SynthesisError;

use super::PairingVar as PG;

use crate::{
    fields::{fp::FpVar, fp12::Fp12Var, fp2::Fp2Var, FieldVar},
    groups::bls12::{G1AffineVar, G1PreparedVar, G1Var, G2PreparedVar, G2Var},
};
use algebra::{
    curves::bls12::{Bls12, Bls12Parameters, TwistType},
    fields::BitIteratorBE,
};
use core::marker::PhantomData;

/// Specifies the constraints for computing a pairing in a BLS12 bilinear group.
pub struct PairingVar<P: Bls12Parameters>(PhantomData<P>);

type Fp2V<P> = Fp2Var<<P as Bls12Parameters>::Fp2Params>;

impl<P: Bls12Parameters> PairingVar<P> {
    // Evaluate the line function at point p.
    #[tracing::instrument(target = "r1cs")]
    fn ell(
        f: &mut Fp12Var<P::Fp12Params>,
        coeffs: &(Fp2V<P>, Fp2V<P>),
        p: &G1AffineVar<P>,
    ) -> Result<(), SynthesisError> {
        let zero = FpVar::<P::Fp>::zero();

        match P::TWIST_TYPE {
            TwistType::M => {
                let c0 = coeffs.0.clone();
                let mut c1 = coeffs.1.clone();
                let c2 = Fp2V::<P>::new(p.y.clone(), zero);

                c1.c0 = c1.c0 * &p.x;
                c1.c1 = c1.c1 * &p.x;
                *f = f.mul_by_014(&c0, &c1, &c2)?;
                Ok(())
            }
            TwistType::D => {
                let c0 = Fp2V::<P>::new(p.y.clone(), zero);
                let mut c1 = coeffs.0.clone();
                let c2 = coeffs.1.clone();

                c1.c0 = c1.c0 * &p.x;
                c1.c1 = c1.c1 * &p.x;
                *f = f.mul_by_034(&c0, &c1, &c2)?;
                Ok(())
            }
        }
    }

    #[tracing::instrument(target = "r1cs")]
    fn exp_by_x(f: &Fp12Var<P::Fp12Params>) -> Result<Fp12Var<P::Fp12Params>, SynthesisError> {
        let mut result = f.optimized_cyclotomic_exp(P::X)?;
        if P::X_IS_NEGATIVE {
            result = result.unitary_inverse()?;
        }
        Ok(result)
    }
}

impl<P: Bls12Parameters> PG<Bls12<P>, P::Fp> for PairingVar<P> {
    type G1Var = G1Var<P>;
    type G2Var = G2Var<P>;
    type G1PreparedVar = G1PreparedVar<P>;
    type G2PreparedVar = G2PreparedVar<P>;
    type GTVar = Fp12Var<P::Fp12Params>;

    #[tracing::instrument(target = "r1cs")]
    fn miller_loop(
        ps: &[Self::G1PreparedVar],
        qs: &[Self::G2PreparedVar],
    ) -> Result<Self::GTVar, SynthesisError> {
        let mut pairs = vec![];
        for (p, q) in ps.iter().zip(qs.iter()) {
            pairs.push((p, q.ell_coeffs.iter()));
        }
        let mut f = Self::GTVar::one();

        for i in BitIteratorBE::new(P::X).skip(1) {
            f.square_in_place()?;

            for &mut (p, ref mut coeffs) in pairs.iter_mut() {
                Self::ell(&mut f, coeffs.next().unwrap(), &p.0)?;
            }

            if i {
                for &mut (p, ref mut coeffs) in pairs.iter_mut() {
                    Self::ell(&mut f, &coeffs.next().unwrap(), &p.0)?;
                }
            }
        }

        if P::X_IS_NEGATIVE {
            f = f.unitary_inverse()?;
        }

        Ok(f)
    }

    #[tracing::instrument(target = "r1cs")]
    fn final_exponentiation(f: &Self::GTVar) -> Result<Self::GTVar, SynthesisError> {
        // Computing the final exponentation following
        // https://eprint.iacr.org/2016/130.pdf.
        // We don't use their "faster" formula because it is difficult to make
        // it work for curves with odd `P::X`.
        // Hence we implement the slower algorithm from Table 1 below.

        let f1 = f.frobenius_map(6)?;

        f.inverse().and_then(|mut f2| {
            // f2 = f^(-1);
            // r = f^(p^6 - 1)
            let mut r = f1;
            r *= &f2;

            // f2 = f^(p^6 - 1)
            f2 = r.clone();
            // r = f^((p^6 - 1)(p^2))
            r.frobenius_map_in_place(2)?;

            // r = f^((p^6 - 1)(p^2) + (p^6 - 1))
            // r = f^((p^6 - 1)(p^2 + 1))
            r *= &f2;

            // Hard part of the final exponentation is below:
            // From https://eprint.iacr.org/2016/130.pdf, Table 1
            let mut y0 = r.cyclotomic_square()?;
            y0 = y0.unitary_inverse()?;

            let mut y5 = Self::exp_by_x(&r)?;

            let mut y1 = y5.cyclotomic_square()?;
            let mut y3 = y0 * &y5;
            y0 = Self::exp_by_x(&y3)?;
            let y2 = Self::exp_by_x(&y0)?;
            let mut y4 = Self::exp_by_x(&y2)?;
            y4 *= &y1;
            y1 = Self::exp_by_x(&y4)?;
            y3 = y3.unitary_inverse()?;
            y1 *= &y3;
            y1 *= &r;
            y3 = r.clone();
            y3 = y3.unitary_inverse()?;
            y0 *= &r;
            y0.frobenius_map_in_place(3)?;
            y4 *= &y3;
            y4.frobenius_map_in_place(1)?;
            y5 *= &y2;
            y5.frobenius_map_in_place(2)?;
            y5 *= &y0;
            y5 *= &y4;
            y5 *= &y1;
            Ok(y5)
        })
    }

    #[tracing::instrument(target = "r1cs")]
    fn prepare_g1(p: &Self::G1Var) -> Result<Self::G1PreparedVar, SynthesisError> {
        Self::G1PreparedVar::from_group_var(p)
    }

    #[tracing::instrument(target = "r1cs")]
    fn prepare_g2(q: &Self::G2Var) -> Result<Self::G2PreparedVar, SynthesisError> {
        Self::G2PreparedVar::from_group_var(q)
    }
}
