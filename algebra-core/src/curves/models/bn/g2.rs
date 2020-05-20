use crate::{
    bytes::ToBytes,
    curves::{
        models::bn::BnParameters,
        short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        AffineCurve,
    },
    fields::{fp6_3over2::Fp6Parameters, Field, Fp2},
};
use derivative::Derivative;
use std::{
    io::{Result as IoResult, Write},
    ops::{AddAssign, MulAssign, Neg, SubAssign},
};
use num_traits::One;

pub type G2Affine<P> = GroupAffine<<P as BnParameters>::G2Parameters>;
pub type G2Projective<P> = GroupProjective<<P as BnParameters>::G2Parameters>;

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: BnParameters"),
    Debug(bound = "P: BnParameters"),
    PartialEq(bound = "P: BnParameters"),
    Eq(bound = "P: BnParameters")
)]
pub struct G2Prepared<P: BnParameters> {
    // Stores the coefficients of the line evaluations as calculated in
    // https://eprint.iacr.org/2013/722.pdf
    pub ell_coeffs: Vec<(Fp2<P::Fp2Params>, Fp2<P::Fp2Params>, Fp2<P::Fp2Params>)>,
    pub infinity:   bool,
}

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: BnParameters"),
    Copy(bound = "P: BnParameters"),
    Debug(bound = "P: BnParameters")
)]
struct G2HomProjective<P: BnParameters> {
    x: Fp2<P::Fp2Params>,
    y: Fp2<P::Fp2Params>,
    z: Fp2<P::Fp2Params>,
}

impl<P: BnParameters> Default for G2Prepared<P> {
    fn default() -> Self {
        Self::from(G2Affine::<P>::prime_subgroup_generator())
    }
}

impl<P: BnParameters> ToBytes for G2Prepared<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        for coeff in &self.ell_coeffs {
            coeff.0.write(&mut writer)?;
            coeff.1.write(&mut writer)?;
            coeff.2.write(&mut writer)?;
        }
        self.infinity.write(writer)
    }
}

impl<P: BnParameters> From<G2Affine<P>> for G2Prepared<P> {
    fn from(q: G2Affine<P>) -> Self {
        let mut coeffs = vec![];
        let mut r: G2HomProjective<P> = G2HomProjective {
            x: q.x,
            y: q.y,
            z: Fp2::one(),
        };

        let negq = q.neg();

        for i in (1..P::SIX_U_PLUS_2_NAF.len()).rev() {
            coeffs.push(doubling_step(&mut r));
            let x = P::SIX_U_PLUS_2_NAF[i - 1];
            match x {
                1 => {
                    coeffs.push(addition_step(&mut r, &q));
                },
                -1 => {
                    coeffs.push(addition_step(&mut r, &negq));
                },
                _ => continue,
            }
        }

        let mut q1 = q;

        q1.x.c1 = q1.x.c1.neg();
        q1.x.mul_assign(&P::Fp6Params::FROBENIUS_COEFF_FP6_C1[1]);

        q1.y.c1 = q1.y.c1.neg();
        q1.y.mul_assign(&P::CUBIC_NONRESIDUE_TO_Q_MINUS_1_OVER_2);

        coeffs.push(addition_step(&mut r, &q1));

        let mut minusq2 = q;
        minusq2
            .x
            .mul_assign(&P::Fp6Params::FROBENIUS_COEFF_FP6_C1[2]);

        coeffs.push(addition_step(&mut r, &minusq2));

        Self {
            ell_coeffs: coeffs,
            infinity:   false,
        }
    }
}

impl<P: BnParameters> G2Prepared<P> {
    pub fn is_zero(&self) -> bool {
        self.infinity
    }
}

fn doubling_step<B: BnParameters>(
    r: &mut G2HomProjective<B>,
) -> (Fp2<B::Fp2Params>, Fp2<B::Fp2Params>, Fp2<B::Fp2Params>) {
    // Adaptation of Algorithm 26, https://eprint.iacr.org/2010/354.pdf
    let mut tmp0 = r.x;
    tmp0.square_in_place();

    let mut tmp1 = r.y;
    tmp1.square_in_place();

    let mut tmp2 = tmp1;
    tmp2.square_in_place();

    let mut tmp3 = tmp1;
    tmp3.add_assign(&r.x);
    tmp3.square_in_place();
    tmp3.sub_assign(&tmp0);
    tmp3.sub_assign(&tmp2);
    tmp3.double_in_place();

    let mut tmp4 = tmp0;
    tmp4.double_in_place();
    tmp4.add_assign(&tmp0);

    let mut tmp6 = r.x;
    tmp6.add_assign(&tmp4);

    let mut tmp5 = tmp4;
    tmp5.square_in_place();

    let zsquared = r.z.square();

    r.x = tmp5;
    r.x.sub_assign(&tmp3);
    r.x.sub_assign(&tmp3);

    r.z.add_assign(&r.y);
    r.z.square_in_place();
    r.z.sub_assign(&tmp1);
    r.z.sub_assign(&zsquared);

    r.y = tmp3;
    r.y.sub_assign(&r.x);
    r.y.mul_assign(&tmp4);

    tmp2.double_in_place();
    tmp2.double_in_place();
    tmp2.double_in_place();

    r.y.sub_assign(&tmp2);

    // up to here everything was by algorith, line 11
    // use R instead of new T

    // tmp3 is the first part of line 12
    tmp3 = tmp4;
    tmp3.mul_assign(&zsquared);
    tmp3.double_in_place();
    tmp3 = tmp3.neg();

    // tmp6 is from line 14
    tmp6.square_in_place();
    tmp6.sub_assign(&tmp0);
    tmp6.sub_assign(&tmp5);

    tmp1.double_in_place();
    tmp1.double_in_place();

    tmp6.sub_assign(&tmp1);

    // tmp0 is the first part of line 16
    tmp0 = r.z;
    tmp0.mul_assign(&zsquared);
    tmp0.double_in_place();

    (tmp0, tmp3, tmp6)
}

fn addition_step<B: BnParameters>(
    r: &mut G2HomProjective<B>,
    q: &G2Affine<B>,
) -> (Fp2<B::Fp2Params>, Fp2<B::Fp2Params>, Fp2<B::Fp2Params>) {
    // Adaptation of Algorithm 27, https://eprint.iacr.org/2010/354.pdf
    let zsquared = r.z.square();

    let ysquared = q.y.square();

    // t0 corresponds to line 1
    let mut t0 = zsquared;
    t0.mul_assign(&q.x);

    // t1 corresponds to lines 2 and 3
    let mut t1 = q.y;
    t1.add_assign(&r.z);
    t1.square_in_place();
    t1.sub_assign(&ysquared);
    t1.sub_assign(&zsquared);
    t1.mul_assign(&zsquared);

    // t2 corresponds to line 4
    let mut t2 = t0;
    t2.sub_assign(&r.x);

    // t3 corresponds to line 5
    let mut t3 = t2;
    t3.square_in_place();

    // t4 corresponds to line 6
    let mut t4 = t3;
    t4.double_in_place();
    t4.double_in_place();

    // t5 corresponds to line 7
    let mut t5 = t4;
    t5.mul_assign(&t2);

    // t6 corresponds to line 8
    let mut t6 = t1;
    t6.sub_assign(&r.y);
    t6.sub_assign(&r.y);

    // t9 corresponds to line 9
    let mut t9 = t6;
    t9.mul_assign(&q.x);

    // corresponds to line 10
    let mut t7 = t4;
    t7.mul_assign(&r.x);

    // corresponds to line 11, but assigns to r.x instead of T.x
    r.x = t6;
    r.x.square_in_place();
    r.x.sub_assign(&t5);
    r.x.sub_assign(&t7);
    r.x.sub_assign(&t7);

    // corresponds to line 12, but assigns to r.z instead of T.z
    r.z.add_assign(&t2);
    r.z.square_in_place();
    r.z.sub_assign(&zsquared);
    r.z.sub_assign(&t3);

    // corresponds to line 13
    let mut t10 = q.y;
    t10.add_assign(&r.z);

    // corresponds to line 14
    let mut t8 = t7;
    t8.sub_assign(&r.x);
    t8.mul_assign(&t6);

    // corresponds to line 15
    t0 = r.y;
    t0.mul_assign(&t5);
    t0.double_in_place();

    // corresponds to line 12, but assigns to r.y instead of T.y
    r.y = t8;
    r.y.sub_assign(&t0);

    // corresponds to line 17
    t10.square_in_place();
    t10.sub_assign(&ysquared);

    let ztsquared = r.z.square();

    t10.sub_assign(&ztsquared);

    // corresponds to line 18
    t9.double_in_place();
    t9.sub_assign(&t10);

    // t10 = 2*Zt from Algo 27, line 19
    t10 = r.z;
    t10.double_in_place();

    // t1 = first multiplicator of line 21
    t6 = t6.neg();

    t1 = t6;
    t1.double_in_place();

    // t9 corresponds to t9 from Algo 27
    (t10, t1, t9)
}
