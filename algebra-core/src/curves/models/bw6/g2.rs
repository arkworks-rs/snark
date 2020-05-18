use crate::{
    bytes::ToBytes,
    curves::{
        bw6::{BW6Parameters, TwistType},
        models::SWModelParameters,
        short_weierstrass_jacobian::{GroupAffine, GroupProjective},
        AffineCurve,
    },
    fields::{BitIterator, Field},
    io::{Result as IoResult, Write},
    Vec,
};
use num_traits::{One, Zero};

pub type G2Affine<P> = GroupAffine<<P as BW6Parameters>::G2Parameters>;
pub type G2Projective<P> = GroupProjective<<P as BW6Parameters>::G2Parameters>;

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: BW6Parameters"),
    Debug(bound = "P: BW6Parameters"),
    PartialEq(bound = "P: BW6Parameters"),
    Eq(bound = "P: BW6Parameters")
)]
pub struct G2Prepared<P: BW6Parameters> {
    // Stores the coefficients of the line evaluations as calculated in
    // https://eprint.iacr.org/2013/722.pdf
    pub ell_coeffs_1: Vec<(P::Fp, P::Fp, P::Fp)>,
    pub ell_coeffs_2: Vec<(P::Fp, P::Fp, P::Fp)>,
    pub infinity: bool,
}

#[derive(Derivative)]
#[derivative(
    Clone(bound = "P: BW6Parameters"),
    Copy(bound = "P: BW6Parameters"),
    Debug(bound = "P: BW6Parameters")
)]
struct G2HomProjective<P: BW6Parameters> {
    x: P::Fp,
    y: P::Fp,
    z: P::Fp,
}

impl<P: BW6Parameters> Default for G2Prepared<P> {
    fn default() -> Self {
        Self::from(G2Affine::<P>::prime_subgroup_generator())
    }
}

impl<P: BW6Parameters> ToBytes for G2Prepared<P> {
    fn write<W: Write>(&self, mut writer: W) -> IoResult<()> {
        for coeff_1 in &self.ell_coeffs_1 {
            coeff_1.0.write(&mut writer)?;
            coeff_1.1.write(&mut writer)?;
            coeff_1.2.write(&mut writer)?;
        }
        for coeff_2 in &self.ell_coeffs_2 {
            coeff_2.0.write(&mut writer)?;
            coeff_2.1.write(&mut writer)?;
            coeff_2.2.write(&mut writer)?;
        }
        self.infinity.write(writer)
    }
}

impl<P: BW6Parameters> From<G2Affine<P>> for G2Prepared<P> {
    fn from(q: G2Affine<P>) -> Self {
        let two_inv = P::Fp::one().double().inverse().unwrap();
        if q.is_zero() {
            return Self {
                ell_coeffs_1: vec![],
                ell_coeffs_2: vec![],
                infinity: true,
            };
        }

        // f_{u+1,Q}(P)
        let mut ell_coeffs_1 = vec![];
        let mut r_1 = G2HomProjective {
            x: q.x,
            y: q.y,
            z: P::Fp::one(),
        };

        for i in BitIterator::new(P::ATE_LOOP_COUNT_1).skip(1) {
            ell_coeffs_1.push(doubling_step::<P>(&mut r_1, &two_inv));

            if i {
                ell_coeffs_1.push(addition_step::<P>(&mut r_1, &q));
            }
        }

        // f_{u^3-u^2-u,Q}(P)
        let mut ell_coeffs_2 = vec![];
        let mut r_2 = G2HomProjective {
            x: q.x,
            y: q.y,
            z: P::Fp::one(),
        };

        for j in BitIterator::new(P::ATE_LOOP_COUNT_2).skip(1) {
            ell_coeffs_2.push(doubling_step::<P>(&mut r_2, &two_inv));

            if j {
                ell_coeffs_2.push(addition_step::<P>(&mut r_2, &q));
            }
        }

        Self {
            ell_coeffs_1,
            ell_coeffs_2,
            infinity: false,
        }
    }
}
impl<P: BW6Parameters> G2Prepared<P> {
    pub fn is_zero(&self) -> bool {
        self.infinity
    }
}

fn doubling_step<B: BW6Parameters>(
    r: &mut G2HomProjective<B>,
    two_inv: &B::Fp,
) -> (B::Fp, B::Fp, B::Fp) {
    // Formula for line function when working with
    // homogeneous projective coordinates.

    let mut a = r.x * &r.y;
    a *= &two_inv;
    let b = r.y.square();
    let c = r.z.square();
    let e = B::G2Parameters::COEFF_B * &(c.double() + &c);
    let f = e.double() + &e;
    let mut g = b + &f;
    g *= &two_inv;
    let h = (r.y + &r.z).square() - &(b + &c);
    let i = e - &b;
    let j = r.x.square();
    let e_square = e.square();

    r.x = a * &(b - &f);
    r.y = g.square() - &(e_square.double() + &e_square);
    r.z = b * &h;
    match B::TWIST_TYPE {
        TwistType::M => (i, -B::TWIST * &h, j.double() + &j),
        TwistType::D => (B::TWIST * &i, -h, j.double() + &j),
    }
}

fn addition_step<B: BW6Parameters>(
    r: &mut G2HomProjective<B>,
    q: &G2Affine<B>,
) -> (B::Fp, B::Fp, B::Fp) {
    // Formula for line function when working with
    // homogeneous projective coordinates.
    let theta = r.y - &(q.y * &r.z);
    let lambda = r.x - &(q.x * &r.z);
    let c = theta.square();
    let d = lambda.square();
    let e = lambda * &d;
    let f = r.z * &c;
    let g = r.x * &d;
    let h = e + &f - &g.double();
    r.x = lambda * &h;
    r.y = theta * &(g - &h) - &(e * &r.y);
    r.z *= &e;
    let j = theta * &q.x - &(lambda * &q.y);

    match B::TWIST_TYPE {
        TwistType::M => (j, -theta, B::TWIST * &lambda),
        TwistType::D => (B::TWIST * &j, -theta, lambda),
    }
}
