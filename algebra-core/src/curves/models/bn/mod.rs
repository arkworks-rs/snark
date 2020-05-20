use crate::{
    curves::{
        models::{ModelParameters, SWModelParameters},
        PairingEngine,
    },
    fields::{
        fp12_2over3over2::{Fp12, Fp12Parameters},
        fp2::Fp2Parameters,
        fp6_3over2::Fp6Parameters,
        Field, Fp2, PrimeField, SquareRootField,
    },
};
use derivative::Derivative;
use num_traits::One;

use std::{marker::PhantomData, ops::MulAssign};

pub mod g1;
pub mod g2;

pub trait BnParameters: 'static {
    const SIX_U_PLUS_2_NAF: &'static [i8];
    const U: &'static [u64];

    type Fp: PrimeField + SquareRootField + Into<<Self::Fp as PrimeField>::BigInt>;
    type Fp2Params: Fp2Parameters<Fp = Self::Fp>;
    type Fp6Params: Fp6Parameters<Fp2Params = Self::Fp2Params>;
    type Fp12Params: Fp12Parameters<Fp6Params = Self::Fp6Params>;
    type G1Parameters: SWModelParameters<BaseField = Self::Fp>;
    type G2Parameters: SWModelParameters<
        BaseField = Fp2<Self::Fp2Params>,
        ScalarField = <Self::G1Parameters as ModelParameters>::ScalarField,
    >;

    const CUBIC_NONRESIDUE_TO_Q_MINUS_1_OVER_2: Fp2<Self::Fp2Params>;
}

#[derive(Derivative)]
#[derivative(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct Bn<P: BnParameters>(PhantomData<fn() -> P>);

pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};

impl<P: BnParameters> Bn<P> {
    // Final steps of the line function on prepared coefficients
    fn ell(
        f: &mut Fp12<P::Fp12Params>,
        coeffs: &(Fp2<P::Fp2Params>, Fp2<P::Fp2Params>, Fp2<P::Fp2Params>),

        p: &G1Affine<P>,
    ) {
        let mut c0 = coeffs.0;
        let mut c1 = coeffs.1;

        c0.c0.mul_assign(&p.y);
        c0.c1.mul_assign(&p.y);

        c1.c0.mul_assign(&p.x);
        c1.c1.mul_assign(&p.x);

        // Sparse multiplication in Fq12
        f.mul_by_034(&c0, &c1, &coeffs.2);
    }

    fn exp_by_x(f: &mut Fp12<P::Fp12Params>) {
        *f = f.pow(&P::U);
    }
}

impl<P: BnParameters> PairingEngine for Bn<P>
/*
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
    >, */
{
    type Fr = <P::G1Parameters as ModelParameters>::ScalarField;
    type G1Projective = G1Projective<P>;
    type G1Affine = G1Affine<P>;
    type G1Prepared = G1Prepared<P>;
    type G2Projective = G2Projective<P>;
    type G2Affine = G2Affine<P>;
    type G2Prepared = G2Prepared<P>;
    type Fq = P::Fp;
    type Fqe = Fp2<P::Fp2Params>;
    type Fqk = Fp12<P::Fp12Params>;

    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>,
    {
        let mut pairs = vec![];
        for (p, q) in i {
            if !p.is_zero() && !q.is_zero() {
                pairs.push((p, q.ell_coeffs.iter()));
            }
        }

        let mut f = Self::Fqk::one();

        for i in (1..P::SIX_U_PLUS_2_NAF.len()).rev() {
            if i != P::SIX_U_PLUS_2_NAF.len() - 1 {
                f.square_in_place();
            }
            for (p, ref mut coeffs) in &mut pairs {
                Self::ell(&mut f, coeffs.next().unwrap(), &p.0);
            }
            let x = P::SIX_U_PLUS_2_NAF[i - 1];
            match x {
                1 => {
                    for (p, ref mut coeffs) in &mut pairs {
                        Self::ell(&mut f, coeffs.next().unwrap(), &p.0);
                    }
                },
                -1 => {
                    for (p, ref mut coeffs) in &mut pairs {
                        Self::ell(&mut f, coeffs.next().unwrap(), &p.0);
                    }
                },
                _ => continue,
            }
        }

        // two additional steps: for q1 and minus q2

        for (p, ref mut coeffs) in &mut pairs {
            Self::ell(&mut f, coeffs.next().unwrap(), &p.0);
        }

        for (p, ref mut coeffs) in &mut pairs {
            Self::ell(&mut f, coeffs.next().unwrap(), &p.0);
        }

        for (_p, ref mut coeffs) in &mut pairs {
            assert_eq!(coeffs.next(), None);
        }

        f
    }

    fn final_exponentiation(r: &Self::Fqk) -> Option<Self::Fqk> {
        let mut f1 = *r;
        f1.conjugate();

        match r.inverse() {
            Some(mut f2) => {
                let mut r = f1;
                r.mul_assign(&f2);
                f2 = r;
                r.frobenius_map(2);
                r.mul_assign(&f2);

                let mut fp = r;
                fp.frobenius_map(1);

                let mut fp2 = r;
                fp2.frobenius_map(2);
                let mut fp3 = fp2;
                fp3.frobenius_map(1);

                let mut fu = r;
                Self::exp_by_x(&mut fu);

                let mut fu2 = fu;
                Self::exp_by_x(&mut fu2);

                let mut fu3 = fu2;
                Self::exp_by_x(&mut fu3);

                let mut y3 = fu;
                y3.frobenius_map(1);

                let mut fu2p = fu2;
                fu2p.frobenius_map(1);

                let mut fu3p = fu3;
                fu3p.frobenius_map(1);

                let mut y2 = fu2;
                y2.frobenius_map(2);

                let mut y0 = fp;
                y0.mul_assign(&fp2);
                y0.mul_assign(&fp3);

                let mut y1 = r;
                y1.conjugate();

                let mut y5 = fu2;
                y5.conjugate();

                y3.conjugate();

                let mut y4 = fu;
                y4.mul_assign(&fu2p);
                y4.conjugate();

                let mut y6 = fu3;
                y6.mul_assign(&fu3p);
                y6.conjugate();

                y6.square_in_place();
                y6.mul_assign(&y4);
                y6.mul_assign(&y5);

                let mut t1 = y3;
                t1.mul_assign(&y5);
                t1.mul_assign(&y6);

                y6.mul_assign(&y2);

                t1.square_in_place();
                t1.mul_assign(&y6);
                t1.square_in_place();

                let mut t0 = t1;
                t0.mul_assign(&y1);

                t1.mul_assign(&y0);

                t0.square_in_place();
                t0.mul_assign(&t1);

                Some(t0)
            },
            None => None,
        }
    }
}
