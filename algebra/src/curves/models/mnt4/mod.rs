use crate::{Fp2, BigInteger768 as BigInteger, PrimeField, SquareRootField, Fp2Parameters, Fp4Parameters, SWModelParameters, ModelParameters, PairingEngine, Fp4, PairingCurve, ProjectiveCurve, Field};
use std::marker::PhantomData;
use std::ops::{Add, Mul, Sub};

pub trait MNT4Parameters: 'static {
    const ATE_LOOP_COUNT: &'static [u64];
    const WNAF: &'static [i32];
    const ATE_IS_LOOP_COUNT_NEG: bool;
    const TWIST: Fp2<Self::Fp2Params>;
    const TWIST_COEFF_A: Fp2<Self::Fp2Params>;
    const FINAL_EXPONENT_LAST_CHUNK_1: BigInteger;
    const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger;
    const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool;
    type Fp: PrimeField + SquareRootField + Into<<Self::Fp as PrimeField>::BigInt>;
    type Fp2Params: Fp2Parameters<Fp = Self::Fp>;
    type Fp4Params: Fp4Parameters<Fp2Params = Self::Fp2Params>;
    type G1Parameters: SWModelParameters<BaseField = Self::Fp>;
    type G2Parameters: SWModelParameters<
        BaseField = Fp2<Self::Fp2Params>,
        ScalarField = <Self::G1Parameters as ModelParameters>::ScalarField,
    >;
}

pub mod g1;
pub mod g2;

pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};
use crate::curves::models::mnt4::g2::G2PreparedCoefficients;

#[derive(Derivative)]
#[derivative(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct MNT4p<P: MNT4Parameters>(PhantomData<fn() -> P>);

impl<P: MNT4Parameters> MNT4p<P> {
    /// Takes as input a point in G1 in projective coordinates, and outputs a
    /// precomputed version of it for pairing purposes.
    fn ate_precompute_g1(value: &G1Projective<P>) -> G1Prepared<P> {
        let g1 = value.into_affine();
        let mut py_twist_squared = P::TWIST.square();
        py_twist_squared.mul_by_fp(&g1.y);

        G1Prepared {p: g1, py_twist_squared}
    }

    /// Takes as input a point in `G2` in projective coordinates, and outputs a
    /// precomputed version of it for pairing purposes.
    fn ate_precompute_g2(value: &G2Projective<P>) -> G2Prepared<P> {

        let mut g2p = G2Prepared {
            q: value.into_affine(),
            coeffs: vec![],
        };

        let mut s = value.into_affine();

        for &n in P::WNAF.iter().rev() {

            //Doubling step
            let gamma = {
                let sx_squared = s.x.square();
                let three_sx_squared_plus_a = sx_squared.double().add(&sx_squared).add(&P::TWIST_COEFF_A);
                let two_sy_inv = s.y.double().inverse().unwrap();
                three_sx_squared_plus_a.mul(&two_sy_inv)
            };
            let gamma_x = gamma.mul(&s.x);
            let new_sx = {
                let two_sx = s.x.double();
                gamma.square().sub(&two_sx)
            };
            let new_sy = {
                let sx_minus_new_sx = s.x.sub(&new_sx);
                gamma.mul(&sx_minus_new_sx).sub(&s.y)
            };
            let c = G2PreparedCoefficients{r_y: s.y, gamma, gamma_x};
            g2p.coeffs.push(c);
            s.x = new_sx;
            s.y = new_sy;

            if n != 0 {
                //Addition step
                let sx_minus_x_inv = s.x.sub(&value.x).inverse().unwrap();
                let numerator = if n > 0  { s.y.sub(&value.y) } else { s.y.add(&value.y) };
                let gamma = numerator.mul(&sx_minus_x_inv);
                let gamma_x = gamma.mul(&value.x);
                let new_sx = {
                    let sx_plus_x = s.x.add(&value.x);
                    gamma.square().sub(&sx_plus_x)
                };
                let new_sy = {
                    let sx_minus_new_sx = s.x.sub(&new_sx);
                    gamma.mul(&sx_minus_new_sx).sub(&s.y)
                };
                let c = G2PreparedCoefficients{r_y: s.y, gamma, gamma_x};
                g2p.coeffs.push(c);
                s.x = new_sx;
                s.y = new_sy;
            }
        }
        g2p
    }


    pub fn ate_miller_loop(p: &G1Prepared<P>, q: &G2Prepared<P>) -> Fp4<P::Fp4Params> {

        let mut f = Fp4::<P::Fp4Params>::one();

        let mut idx: usize = 0;


        for &n in P::WNAF.iter().rev() {
            // code below gets executed for all bits (EXCEPT the MSB itself) of
            // mnt4_param_p (skipping leading zeros) in MSB to LSB order

            f = f.square();
            let c = &q.coeffs[idx];
            idx += 1;

            let mut gamma_twist_times_x = c.gamma.mul(&P::TWIST);
            gamma_twist_times_x.mul_by_fp(&p.p.x);
            let g_rr_at_p = Fp4::<P::Fp4Params>::new(
                p.py_twist_squared,
                c.gamma_x - &gamma_twist_times_x  -&c.r_y,
            );

            f = f.mul_by_023(&g_rr_at_p);

            if n != 0 {
                let c = &q.coeffs[idx];
                idx += 1;

                let mut gamma_twist_times_x = c.gamma.mul(&P::TWIST);
                gamma_twist_times_x.mul_by_fp(&p.p.x);
                let g_rq_at_p_c1 = if n > 0 {
                    c.gamma_x - &gamma_twist_times_x - &q.q.y
                } else {
                    c.gamma_x - &gamma_twist_times_x + &q.q.y
                };
                let g_rq_at_p = Fp4::<P::Fp4Params>::new(
                    p.py_twist_squared,
                    g_rq_at_p_c1,
                );
                f = f.mul_by_023(&g_rq_at_p);
            }
        }

        if P::ATE_IS_LOOP_COUNT_NEG {
            f = f.unitary_inverse();
        }

        f
    }

    pub fn final_exponentiation(value: &Fp4<P::Fp4Params>) -> Fp4<P::Fp4Params> {
        let value_inv = value.inverse().unwrap();
        let value_to_first_chunk = Self::final_exponentiation_first_chunk(value, &value_inv);
        let value_inv_to_first_chunk = Self::final_exponentiation_first_chunk(&value_inv, value);
        Self::final_exponentiation_last_chunk(&value_to_first_chunk, &value_inv_to_first_chunk)
    }

    fn final_exponentiation_first_chunk(elt: &Fp4<P::Fp4Params>, elt_inv: &Fp4<P::Fp4Params>) -> Fp4<P::Fp4Params> {

        let mut elt_q2 = elt.clone();
        elt_q2.frobenius_map(2);
        let elt_q2_over_elt = elt_q2 * &elt_inv;
        elt_q2_over_elt
    }

    //Checked
    fn final_exponentiation_last_chunk(elt: &Fp4<P::Fp4Params>, elt_inv: &Fp4<P::Fp4Params>) -> Fp4<P::Fp4Params> {
        let elt_clone = elt.clone();
        let elt_inv_clone = elt_inv.clone();

        let mut elt_q = elt.clone();
        elt_q.frobenius_map(1);

        let w1_part = elt_q.cyclotomic_exp(&P::FINAL_EXPONENT_LAST_CHUNK_1);
        let w0_part;
        if P::FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG {
            w0_part = elt_inv_clone.cyclotomic_exp(&P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0);
        } else {
            w0_part = elt_clone.cyclotomic_exp(&P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0);
        }

        w1_part * &w0_part
    }
}

impl<P: MNT4Parameters> PairingEngine for MNT4p<P>
    where
        G1Affine<P>: PairingCurve<
            BaseField = <P::G1Parameters as ModelParameters>::BaseField,
            ScalarField = <P::G1Parameters as ModelParameters>::ScalarField,
            Projective = G1Projective<P>,
            PairWith = G2Affine<P>,
            Prepared = G1Prepared<P>,
            PairingResult = Fp4<P::Fp4Params>,
        >,
        G2Affine<P>: PairingCurve<
            BaseField = <P::G2Parameters as ModelParameters>::BaseField,
            ScalarField = <P::G1Parameters as ModelParameters>::ScalarField,
            Projective = G2Projective<P>,
            PairWith = G1Affine<P>,
            Prepared = G2Prepared<P>,
            PairingResult = Fp4<P::Fp4Params>,
        >,

{
    type Fr = <P::G1Parameters as ModelParameters>::ScalarField;
    type G1Projective = G1Projective<P>;
    type G1Affine = G1Affine<P>;
    type G2Projective = G2Projective<P>;
    type G2Affine = G2Affine<P>;
    type Fq = P::Fp;
    type Fqe = Fp2<P::Fp2Params>;
    type Fqk = Fp4<P::Fp4Params>;

    fn miller_loop<'a, I>(i: I) -> Self::Fqk
        where
            I: IntoIterator<
                Item = &'a (
                    &'a <Self::G1Affine as PairingCurve>::Prepared,
                    &'a <Self::G2Affine as PairingCurve>::Prepared,
                ),
            >,
    {
        let mut result = Self::Fqk::one();
        for &(ref p, ref q) in i {
            result *= &Self::ate_miller_loop(p, q);
        }
        result
    }

    fn final_exponentiation(r: &Self::Fqk) -> Option<Self::Fqk> {
        Some(Self::final_exponentiation(r))
    }
}

