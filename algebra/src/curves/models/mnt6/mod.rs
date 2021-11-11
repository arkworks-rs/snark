use crate::{
    BigInteger768 as BigInteger, Error, Field, Fp3, Fp3Parameters, Fp6, Fp6Parameters,
    ModelParameters, PairingEngine, PrimeField, SWModelParameters, SquareRootField,
};
use std::marker::PhantomData;
use std::ops::{Add, Mul, Sub};

// Ate pairing e: G_1 x G_2 -> G_T for MNT6 curves over prime fields
//
//     E: y^2 = x^3 + a*x + b mod p.
//
// Its embedding field F6 is regarded as towered extension
//
//     F6 = F2[Y]/(Y^2-X),
//     F3 = Fp[X]/(X^3-alpha),
//
// using a "non-residue" alpha mod p such that (X^6-alpha) is irreducible over Fp.
// We apply standard efficiency measures (see, e.g. ): G_2 is represented by a subgroup
// of prime order r=ord(G_1) of the quadratic twist
//
//     E': y^2 = x^3 + (a*twist^2) x + b*twist^3
//
// over F3, with twist = X = Y^2, the Frobenius operator is applied to reduce the cost of the
// final exponentiation, and we do pre-computations of (essentially) the line coefficients
// of the Miller loop.
// The loop count allows signed bit representation, so this variant supports curves with Frobenius
// trace having low Hamming weight NAF..

pub trait MNT6Parameters: 'static {
    // the loop count for the Miller loop, equals the |Frobenius trace of E - 1|
    const ATE_LOOP_COUNT: &'static [u64];
    // the non-adjacent normal form of ATE_LOOP_COUNT trimmed of leading zeroes and
    // without MSB, starting with the least significant bit
    const WNAF: &'static [i32];
    // true/false depending whether the Frobenius trace is negative/positive
    const ATE_IS_LOOP_COUNT_NEG: bool;
    // The twist factor twist=Y^2  for
    // E': y'^2  = x'^3 + a*twist^2*x + twist^3 * b
    // as needed for the point evaluation of the Miller loop lines
    const TWIST: Fp3<Self::Fp3Params>;
    // Weierstrass coefficient a'=a*omega^4= a*alpha of the quadratic twist E'
    // as needed for the point evaluation of the Miller loop lines
    // translated via the twist map
    const TWIST_COEFF_A: Fp3<Self::Fp3Params>;
    // the final pairing exponent is decomposed as
    //      (p^6-1)/r = (p^3-1)(p+1) (p^2 - p + 1)/r,
    // wheras
    //      (p^2 - p + 1)/r = m_1*p + m_0,
    // with 0<= m_0 < p, m_0
    const FINAL_EXPONENT_LAST_CHUNK_1: BigInteger;
    // and m_0, |m_0| <= p/2, equal to
    const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger;
    // is set true/false depending on the sign of m_0
    const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool;

    // base field F of the curve
    type Fp: PrimeField + SquareRootField + Into<<Self::Fp as PrimeField>::BigInt>;
    // scalar field of the curve
    type Fr: PrimeField + SquareRootField + Into<<Self::Fr as PrimeField>::BigInt>;
    // parameters of the quadratic extension field F3
    type Fp3Params: Fp3Parameters<Fp = Self::Fp>;
    // paramters of the embedding field F6
    type Fp6Params: Fp6Parameters<Fp3Params = Self::Fp3Params>;
    // parameters for E with defining field F
    type G1Parameters: SWModelParameters<BaseField = Self::Fp, ScalarField = Self::Fr>;
    // parameters for the quadratic twist E' over F3
    type G2Parameters: SWModelParameters<
        BaseField = Fp3<Self::Fp3Params>,
        ScalarField = <Self::G1Parameters as ModelParameters>::ScalarField,
    >;
}

pub mod g1;
pub mod g2;

pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};
use crate::curves::models::mnt6::g2::G2PreparedCoefficients;

#[derive(Derivative)]
#[derivative(Copy, Clone, PartialEq, Eq, Debug, Hash)]
pub struct MNT6p<P: MNT6Parameters>(PhantomData<fn() -> P>);

impl<P: MNT6Parameters> MNT6p<P> {
    // Takes as input a point in G1 in affine coordinates, and outputs a
    // precomputed version of it for pairing purposes.
    fn ate_precompute_g1(value: &G1Affine<P>) -> G1Prepared<P> {
        let mut py_twist_squared = P::TWIST.square();
        py_twist_squared.mul_assign_by_fp(&value.y);

        G1Prepared {
            p: *value,
            py_twist_squared,
        }
    }

    // Takes as input a (non-zero) point Q from G2 in affine coordinates, and outputs the
    // (P-independent) pre-computable coefficients for all operations of the Miller loop.
    // These are comprised of the line coefficients in an optimized variant:
    //     s.y = the y-coordinate of internal state S,
    //     gamma = the F3-slope of the tangent/P-chord at S,
    //     gamma_x = the F3-slope times the x-coordinate s.x of S.
    fn ate_precompute_g2(value: &G2Affine<P>) -> Result<G2Prepared<P>, Error> {
        let mut g2p = G2Prepared {
            q: *value,
            coeffs: vec![],
        };

        let mut s = value.clone();

        // signed binary representation of the Ate loop count in big endian order
        for &n in P::WNAF.iter().rev() {
            //Doubling step
            let gamma = {
                let sx_squared = s.x.square();
                let three_sx_squared_plus_a =
                    sx_squared.double().add(&sx_squared).add(&P::TWIST_COEFF_A);
                if value.y.is_zero() {
                    Err(format!("Invalid Q-point value"))?
                }
                let two_sy_inv = s.y.double().inverse().unwrap();
                three_sx_squared_plus_a.mul(&two_sy_inv) // the F3-slope of the tangent at S=(s.x,s.y)
            };
            let gamma_x = gamma.mul(&s.x);
            let new_sx = {
                let two_sx = s.x.double();
                gamma.square().sub(&two_sx) //x-coordinate after doubling
            };
            let new_sy = {
                let sx_minus_new_sx = s.x.sub(&new_sx);
                gamma.mul(&sx_minus_new_sx).sub(&s.y) //y-coordinate after doubling
            };
            let c = G2PreparedCoefficients {
                r_y: s.y,
                gamma,
                gamma_x,
            };
            g2p.coeffs.push(c);
            s.x = new_sx;
            s.y = new_sy;

            if n != 0 {
                //Addition/substraction step depending on the sign of n
                if s.x == value.x {
                    Err(format!("Invalid Q-point value"))?
                }
                let sx_minus_x_inv = s.x.sub(&value.x).inverse().unwrap();
                let numerator = if n > 0 {
                    s.y.sub(&value.y)
                } else {
                    s.y.add(&value.y)
                };
                let gamma = numerator.mul(&sx_minus_x_inv); //the F3 slope of the chord Q'R'
                let gamma_x = gamma.mul(&value.x);
                let new_sx = {
                    let sx_plus_x = s.x.add(&value.x);
                    gamma.square().sub(&sx_plus_x)
                };
                let new_sy = {
                    let sx_minus_new_sx = s.x.sub(&new_sx);
                    gamma.mul(&sx_minus_new_sx).sub(&s.y)
                };
                let c = G2PreparedCoefficients {
                    r_y: s.y,
                    gamma,
                    gamma_x,
                };
                g2p.coeffs.push(c);
                s.x = new_sx;
                s.y = new_sy;
            }
        }

        Ok(g2p)
    }

    pub fn ate_miller_loop(p: &G1Prepared<P>, q: &G2Prepared<P>) -> Fp6<P::Fp6Params> {
        let mut f = Fp6::<P::Fp6Params>::one();

        let mut idx: usize = 0;

        for &n in P::WNAF.iter().rev() {
            // code below gets executed for all bits (EXCEPT the MSB itself) of
            // mnt4_param_p (skipping leading zeros) in MSB to LSB order

            // doubling step
            f = f.square();
            let c = &q.coeffs[idx];
            idx += 1;

            // evaluate the tangent line g_{R,R} at P in F6 (scaled by twist^2) using the
            // pre-computed data:
            //      g_{R,R}(P) = (y_P - lambda*x_p - d) * X^2,
            // where
            //      lambda = gamma * Y/twist,
            //      d = (y'-gamma * x')* Y/twist^2,
            // with (x',y') being the twist coordinates of R.
            // Thus
            //     g_{R,R}(P) = y_p*twist^2 + (gamma*x'- gamma*twist*x_p - y') *Y.
            // The scale factor twist^2 from F3 is cancelled out by the final exponentiation.

            let mut gamma_twist_times_x = c.gamma.mul(&P::TWIST);
            gamma_twist_times_x.mul_assign_by_fp(&p.p.x);
            let g_rr_at_p = Fp6::<P::Fp6Params>::new(
                p.py_twist_squared,
                c.gamma_x - &gamma_twist_times_x - &c.r_y,
            );
            //and cumulate it to f
            f = f.mul_by_2345(&g_rr_at_p);

            //addition/substraction step
            if n != 0 {
                let c = &q.coeffs[idx];
                idx += 1;

                //evaluate chord g_{RQ}(P) in F6 using pre-computed data as above
                //I suggest to write a separate function for the point evaluation
                //as done in the implementation of the sw6 Miller loop
                let mut gamma_twist_times_x = c.gamma.mul(&P::TWIST);
                gamma_twist_times_x.mul_assign_by_fp(&p.p.x);
                let g_rq_at_p_c1 = if n > 0 {
                    c.gamma_x - &gamma_twist_times_x - &q.q.y
                } else {
                    c.gamma_x - &gamma_twist_times_x + &q.q.y
                };
                let g_rq_at_p = Fp6::<P::Fp6Params>::new(p.py_twist_squared, g_rq_at_p_c1);
                //and cumulate it to f
                f = f.mul_by_2345(&g_rq_at_p);
            }
        }

        if P::ATE_IS_LOOP_COUNT_NEG {
            f = f.unitary_inverse();
        }

        f
    }

    pub fn final_exponentiation(value: &Fp6<P::Fp6Params>) -> Result<Fp6<P::Fp6Params>, Error> {
        if value.is_zero() {
            Err(format!("Invalid exponentiation value: 0"))?
        }
        let value_inv = value.inverse().unwrap();
        // "easy part" of the exponentiation
        let value_to_first_chunk = Self::final_exponentiation_first_chunk(value, &value_inv);
        let value_inv_to_first_chunk = Self::final_exponentiation_first_chunk(&value_inv, value);
        // "hard part"
        Ok(Self::final_exponentiation_last_chunk(
            &value_to_first_chunk,
            &value_inv_to_first_chunk,
        ))
    }

    fn final_exponentiation_first_chunk(
        elt: &Fp6<P::Fp6Params>,
        elt_inv: &Fp6<P::Fp6Params>,
    ) -> Fp6<P::Fp6Params> {
        // use the Frobenius map and elt^{-1} to compute the "easy part"
        // elt^{(q^3-1)*(q+1)}

        let mut elt_q3 = elt.clone();
        // elt^{q^3}
        elt_q3.conjugate();
        // elt^{q^3-1}
        let mut elt_q3_over_elt = elt_q3 * elt_inv;
        let elt_q3_over_elt_clone = elt_q3_over_elt.clone();
        // elt^{(q^3-1)q}
        elt_q3_over_elt.frobenius_map(1);
        // elt^{(q^3-1)*(q+1)}
        elt_q3_over_elt *= &elt_q3_over_elt_clone;

        elt_q3_over_elt
    }

    fn final_exponentiation_last_chunk(
        elt: &Fp6<P::Fp6Params>,
        elt_inv: &Fp6<P::Fp6Params>,
    ) -> Fp6<P::Fp6Params> {
        // remaining exponentiaton by m_1*q + m_0, m_0 can be signed.
        let elt_clone = elt.clone();
        let elt_inv_clone = elt_inv.clone();

        let mut elt_q = elt.clone();
        //elt^{q}
        elt_q.frobenius_map(1);

        // exponentiation by m_1 and m_0 using optimized exponentiation for r-th roots of unity
        //elt^{q*m_1}
        let w1_part = elt_q.cyclotomic_exp(&P::FINAL_EXPONENT_LAST_CHUNK_1);
        //elt^{m_0}
        let w0_part;
        if P::FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG {
            w0_part = elt_inv_clone.cyclotomic_exp(&P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0);
        } else {
            w0_part = elt_clone.cyclotomic_exp(&P::FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0);
        }
        //elt^{q*m_1+m_0}
        w1_part * &w0_part
    }
}

impl<P: MNT6Parameters> PairingEngine for MNT6p<P> {
    type Fr = <P::G1Parameters as ModelParameters>::ScalarField;
    type G1Projective = G1Projective<P>;
    type G1Affine = G1Affine<P>;
    type G1Prepared = G1Prepared<P>;
    type G2Projective = G2Projective<P>;
    type G2Affine = G2Affine<P>;
    type G2Prepared = G2Prepared<P>;
    type Fq = P::Fp;
    type Fqe = Fp3<P::Fp3Params>;
    type Fqk = Fp6<P::Fp6Params>;

    fn miller_loop<'a, I>(i: I) -> Result<Self::Fqk, Error>
    where
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>,
    {
        let mut result = Self::Fqk::one();
        for &(ref p, ref q) in i {
            result *= &Self::ate_miller_loop(p, q);
        }
        Ok(result)
    }

    fn final_exponentiation(r: &Self::Fqk) -> Result<Self::Fqk, Error> {
        Self::final_exponentiation(r)
    }
}
