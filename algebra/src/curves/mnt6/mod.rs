use crate::field_new;
use crate::{
    biginteger::BigInteger320,
    curves::{PairingCurve, PairingEngine, ProjectiveCurve},
    fields::{
        mnt6::{
            fq::{Fq, FqParameters},
            Fq3, Fq6, Fr,
        },
        BitIterator, Field, FpParameters,
    },
};

pub mod g1;
pub mod g2;
#[cfg(test)]
mod tests;

use self::g2::{AteAdditionCoefficients, AteDoubleCoefficients, G2ProjectiveExtended};
pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};

pub type GT = Fq6;

#[derive(Copy, Clone, Debug)]
pub struct MNT6;

impl PairingEngine for MNT6 {
    type Fr = Fr;
    type G1Projective = G1Projective;
    type G1Affine = G1Affine;
    type G2Projective = G2Projective;
    type G2Affine = G2Affine;
    type Fq = Fq;
    type Fqe = Fq3;
    type Fqk = Fq6;

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
            result *= &MNT6::ate_miller_loop(p, q);
        }
        result
    }

    fn final_exponentiation(r: &Self::Fqk) -> Option<Self::Fqk> {
        Some(MNT6::final_exponentiation(r))
    }
}

impl MNT6 {
    /// Takes as input a point in G1 in projective coordinates, and outputs a
    /// precomputed version of it for pairing purposes.
    fn ate_precompute_g1(value: &G1Projective) -> G1Prepared {
        let g1 = value.into_affine();

        let mut x_twist = TWIST.clone();
        x_twist.mul_assign_by_fp(&g1.x);

        let mut y_twist = TWIST.clone();
        y_twist.mul_assign_by_fp(&g1.y);

        G1Prepared {
            x: g1.x,
            y: g1.y,
            x_twist,
            y_twist,
        }
    }

    /// Takes as input a point in `G2` in projective coordinates, and outputs a
    /// precomputed version of it for pairing purposes.
    fn ate_precompute_g2(value: &G2Projective) -> G2Prepared {
        let g2 = value.into_affine();

        let twist_inv = TWIST.inverse().unwrap();

        let mut g2p = G2Prepared {
            x:                     g2.x,
            y:                     g2.y,
            x_over_twist:          g2.x * &twist_inv,
            y_over_twist:          g2.y * &twist_inv,
            double_coefficients:   vec![],
            addition_coefficients: vec![],
        };

        let mut r = G2ProjectiveExtended {
            x: g2.x,
            y: g2.y,
            z: Fq3::one(),
            t: Fq3::one(),
        };

        for (idx, value) in ATE_LOOP_COUNT.iter().rev().enumerate() {
            let mut tmp = *value;
            let skip_extraneous_bits = 64 - value.leading_zeros();
            let mut v = Vec::with_capacity(16);
            for i in 0..64 {
                if idx == 0 && (i == 0 || i >= skip_extraneous_bits) {
                    continue;
                }
                v.push(tmp & 1 == 1);
                tmp >>= 1;
            }

            for bit in v.iter().rev() {
                let (r2, coeff) = MNT6::doubling_step_for_flipped_miller_loop(&r);
                g2p.double_coefficients.push(coeff);
                r = r2;

                if *bit {
                    let (r2, coeff) =
                        MNT6::mixed_addition_step_for_flipped_miller_loop(&g2.x, &g2.y, &r);
                    g2p.addition_coefficients.push(coeff);
                    r = r2;
                }

                tmp >>= 1;
            }
        }

        if ATE_IS_LOOP_COUNT_NEG {
            let rz_inv = r.z.inverse().unwrap();
            let rz2_inv = rz_inv.square();
            let rz3_inv = rz_inv * &rz2_inv;

            let minus_r_affine_x = r.x * &rz2_inv;
            let minus_r_affine_y = -r.y * &rz3_inv;

            let add_result = MNT6::mixed_addition_step_for_flipped_miller_loop(
                &minus_r_affine_x,
                &minus_r_affine_y,
                &r,
            );
            g2p.addition_coefficients.push(add_result.1);
        }

        g2p
    }

    fn doubling_step_for_flipped_miller_loop(
        r: &G2ProjectiveExtended,
    ) -> (G2ProjectiveExtended, AteDoubleCoefficients) {
        let a = r.t.square();
        let b = r.x.square();
        let c = r.y.square();
        let d = c.square();
        let e = (r.x + &c).square() - &b - &d;
        let f = (b + &b + &b) + &(TWIST_COEFF_A * &a);
        let g = f.square();

        let d_eight = d.double().double().double();

        let x = -(e + &e + &e + &e) + &g;
        let y = -d_eight + &(f * &(e + &e - &x));
        let z = (r.y + &r.z).square() - &c - &r.z.square();
        let t = z.square();

        let r2 = G2ProjectiveExtended { x, y, z, t };
        let coeff = AteDoubleCoefficients {
            c_h:  (r2.z + &r.t).square() - &r2.t - &a,
            c_4c: c + &c + &c + &c,
            c_j:  (f + &r.t).square() - &g - &a,
            c_l:  (f + &r.x).square() - &g - &b,
        };

        (r2, coeff)
    }

    fn mixed_addition_step_for_flipped_miller_loop(
        x: &Fq3,
        y: &Fq3,
        r: &G2ProjectiveExtended,
    ) -> (G2ProjectiveExtended, AteAdditionCoefficients) {
        let a = y.square();
        let b = r.t * x;
        let d = ((r.z + y).square() - &a - &r.t) * &r.t;
        let h = b - &r.x;
        let i = h.square();
        let e = i + &i + &i + &i;
        let j = h * &e;
        let v = r.x * &e;
        let l1 = d - &(r.y + &r.y);

        let x = l1.square() - &j - &(v + &v);
        let y = l1 * &(v - &x) - &(j * &(r.y + &r.y));
        let z = (r.z + &h).square() - &r.t - &i;
        let t = z.square();

        let r2 = G2ProjectiveExtended { x, y, z, t };
        let coeff = AteAdditionCoefficients { c_l1: l1, c_rz: z };

        (r2, coeff)
    }

    pub fn ate_miller_loop(p: &G1Prepared, q: &G2Prepared) -> Fq6 {
        let l1_coeff = field_new!(Fq3, p.x, Fq::zero(), Fq::zero()) - &q.x_over_twist;

        let mut f = Fq6::one();

        let mut dbl_idx: usize = 0;
        let mut add_idx: usize = 0;

        let mut found_one = false;

        for bit in BitIterator::new(ATE_LOOP_COUNT) {
            // code below gets executed for all bits (EXCEPT the MSB itself) of
            // mnt6_param_p (skipping leading zeros) in MSB to LSB order
            if !found_one && bit {
                found_one = true;
                continue;
            } else if !found_one {
                continue;
            }

            let dc = &q.double_coefficients[dbl_idx];
            dbl_idx += 1;

            let g_rr_at_p = Fq6::new(
                -dc.c_4c - &(dc.c_j * &p.x_twist) + &dc.c_l,
                dc.c_h * &p.y_twist,
            );

            f = f.square() * &g_rr_at_p;

            if bit {
                let ac = &q.addition_coefficients[add_idx];
                add_idx += 1;

                let g_rq_at_p = Fq6::new(
                    ac.c_rz * &p.y_twist,
                    -(q.y_over_twist * &ac.c_rz + &(l1_coeff * &ac.c_l1)),
                );
                f = f * &g_rq_at_p;
            }
        }

        if ATE_IS_LOOP_COUNT_NEG {
            let ac = &q.addition_coefficients[add_idx];

            let g_rnegr_at_p = Fq6::new(
                ac.c_rz * &p.y_twist,
                -(q.y_over_twist * &ac.c_rz + &(l1_coeff * &ac.c_l1)),
            );
            f = (f * &g_rnegr_at_p).inverse().unwrap();
        }

        f
    }

    pub fn final_exponentiation(value: &Fq6) -> GT {
        let value_inv = value.inverse().unwrap();
        let value_to_first_chunk = MNT6::final_exponentiation_first_chunk(value, &value_inv);
        let value_inv_to_first_chunk = MNT6::final_exponentiation_first_chunk(&value_inv, value);
        MNT6::final_exponentiation_last_chunk(&value_to_first_chunk, &value_inv_to_first_chunk)
    }

    fn final_exponentiation_first_chunk(elt: &Fq6, elt_inv: &Fq6) -> Fq6 {
        // (q^3-1)*(q+1)

        // elt_q3 = elt^(q^3)
        let mut elt_q3 = elt.clone();
        elt_q3.frobenius_map(3);
        // elt_q3_over_elt = elt^(q^3-1)
        let elt_q3_over_elt = elt_q3 * &elt_inv;
        // alpha = elt^((q^3-1) * q)
        let mut alpha = elt_q3_over_elt.clone();
        alpha.frobenius_map(1);
        // beta = elt^((q^3-1)*(q+1)
        alpha * &elt_q3_over_elt
    }

    fn final_exponentiation_last_chunk(elt: &Fq6, elt_inv: &Fq6) -> Fq6 {
        let elt_clone = elt.clone();
        let elt_inv_clone = elt_inv.clone();

        let mut elt_q = elt.clone();
        elt_q.frobenius_map(1);

        let w1_part = elt_q.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_1);
        let w0_part;
        if FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG {
            w0_part = elt_inv_clone.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0);
        } else {
            w0_part = elt_clone.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0);
        }

        w1_part * &w0_part
    }
}

pub const TWIST: Fq3 = field_new!(Fq3, FQ_ZERO, FQ_ONE, FQ_ZERO);
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger320([0, 0, 0, 0, 0]));
pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
pub const TWIST_COEFF_A: Fq3 = field_new!(Fq3, 
    FQ_ZERO,
    FQ_ZERO,
    field_new!(Fq, BigInteger320([
        0xb9b2411bfd0eafef,
        0xc61a10fadd9fecbd,
        0x89f128e59811f3fb,
        0x980c0f780adadabb,
        0x9ba1f11320,
    ])),
);

pub const ATE_LOOP_COUNT: [u64; 3] = [0xdc9a1b671660000, 0x46609756bec2a33f, 0x1eef55];

pub const ATE_IS_LOOP_COUNT_NEG: bool = true;

pub const FINAL_EXPONENT_LAST_CHUNK_1: BigInteger320 = BigInteger320([0x1, 0x0, 0x0, 0x0, 0x0]);

pub const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = true;

pub const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger320 =
    BigInteger320([0xdc9a1b671660000, 0x46609756bec2a33f, 0x1eef55, 0x0, 0x0]);
