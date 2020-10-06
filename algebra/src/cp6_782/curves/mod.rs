use crate::{
    biginteger::BigInteger832,
    curves::{models::SWModelParameters, PairingEngine},
    field_new,
    fields::{BitIteratorBE, Field, FpParameters},
    One,
};

use crate::cp6_782::{Fq, Fq3, Fq6, FqParameters, Fr};

pub mod g1;
pub use self::g1::{G1Affine, G1Projective};

pub mod g2;
pub use self::g2::{G2Affine, G2Projective};

#[cfg(test)]
mod tests;

pub type GT = Fq6;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct CP6_782;

impl PairingEngine for CP6_782 {
    type Fr = Fr;
    type G1Projective = G1Projective;
    type G1Affine = G1Affine;
    type G1Prepared = G1Affine;
    type G2Projective = G2Projective;
    type G2Affine = G2Affine;
    type G2Prepared = G2Affine;
    type Fq = Fq;
    type Fqe = Fq3;
    type Fqk = Fq6;

    fn miller_loop<'a, I>(i: I) -> Self::Fqk
    where
        I: IntoIterator<Item = &'a (Self::G1Prepared, Self::G2Prepared)>,
    {
        let mut result = Self::Fqk::one();
        for &(ref p, ref q) in i {
            result *= &CP6_782::ate_miller_loop(p, q);
        }
        result
    }

    fn final_exponentiation(r: &Self::Fqk) -> Option<Self::Fqk> {
        Some(CP6_782::final_exponentiation(r))
    }
}

impl CP6_782 {
    pub fn ate_pairing(p: &G1Affine, q: &G2Affine) -> GT {
        CP6_782::final_exponentiation(&CP6_782::ate_miller_loop(p, q))
    }

    fn ate_miller_loop(p: &G1Affine, q: &G2Affine) -> Fq6 {
        let px = p.x;
        let py = p.y;
        let qx = q.x;
        let qy = q.y;
        let mut py_twist_squared = TWIST.square();
        py_twist_squared.mul_assign_by_fp(&py);

        let mut old_rx;
        let mut old_ry;
        let mut rx = qx;
        let mut ry = qy;
        let mut f = Fq6::one();

        // The for loop is executed for all bits (EXCEPT the MSB itself) of
        // cp6_782_param_p (skipping leading zeros) in MSB to LSB order
        for bit in BitIteratorBE::without_leading_zeros(ATE_LOOP_COUNT).skip(1) {
            old_rx = rx;
            old_ry = ry;

            let old_rx_square = old_rx.square();
            let old_rx_square_3 = old_rx_square.double() + &old_rx_square;
            let old_rx_square_3_a = old_rx_square_3 + &g2::Parameters::COEFF_A;
            let old_ry_double_inverse = old_ry.double().inverse().unwrap();

            let gamma = old_rx_square_3_a * &old_ry_double_inverse;
            let gamma_twist = gamma * &TWIST;
            let gamma_old_rx = gamma * &old_rx;
            let mut gamma_twist_px = gamma_twist;
            gamma_twist_px.mul_assign_by_fp(&px);

            let x = py_twist_squared;
            let y = gamma_old_rx - &old_ry - &gamma_twist_px;
            let ell_rr_at_p = Fq6::new(x, y);

            rx = gamma.square() - &old_rx.double();
            ry = gamma * &(old_rx - &rx) - &old_ry;
            f = f.square() * &ell_rr_at_p;

            if bit {
                old_rx = rx;
                old_ry = ry;

                let gamma = (old_ry - &qy) * &((old_rx - &qx).inverse().unwrap());
                let gamma_twist = gamma * &TWIST;
                let gamma_qx = gamma * &qx;
                let mut gamma_twist_px = gamma_twist;
                gamma_twist_px.mul_assign_by_fp(&px);

                let x = py_twist_squared;
                let y = gamma_qx - &qy - &gamma_twist_px;
                let ell_rq_at_p = Fq6::new(x, y);

                rx = gamma.square() - &old_rx - &qx;
                ry = gamma * &(old_rx - &rx) - &old_ry;
                f = f * &ell_rq_at_p;
            }
        }
        f
    }

    fn final_exponentiation(value: &Fq6) -> GT {
        let value_inv = value.inverse().unwrap();
        let value_to_first_chunk = CP6_782::final_exponentiation_first(value, &value_inv);
        let value_inv_to_first_chunk = CP6_782::final_exponentiation_first(&value_inv, value);
        CP6_782::final_exponentiation_last(&value_to_first_chunk, &value_inv_to_first_chunk)
    }

    fn final_exponentiation_first(elt: &Fq6, elt_inv: &Fq6) -> Fq6 {
        // (q^3-1)*(q+1)

        // elt_q3 = elt^(q^3)
        let mut elt_q3 = elt.clone();
        elt_q3.frobenius_map(3);
        // elt_q3_over_elt = elt^(q^3-1)
        let elt_q3_over_elt = elt_q3 * elt_inv;
        // alpha = elt^((q^3-1) * q)
        let mut alpha = elt_q3_over_elt.clone();
        alpha.frobenius_map(1);
        // beta = elt^((q^3-1)*(q+1)
        alpha * &elt_q3_over_elt
    }

    fn final_exponentiation_last(elt: &Fq6, elt_inv: &Fq6) -> Fq6 {
        let mut elt_q = elt.clone();
        elt_q.frobenius_map(1);

        let w1_part = elt_q.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_W1);
        let w0_part = if FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG {
            elt_inv.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0)
        } else {
            elt.cyclotomic_exp(&FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0)
        };

        w1_part * &w0_part
    }
}

/// FQ_ZERO = 0
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger832([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));

/// FQ_ONE = 1
pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);

/// TWIST = (0, 1, 0)
pub const TWIST: Fq3 = field_new!(Fq3, FQ_ZERO, FQ_ONE, FQ_ZERO);

/// ATE_IS_LOOP_COUNT_NEG = false
pub const ATE_IS_LOOP_COUNT_NEG: bool = false;

/// ATE_LOOP_COUNT =
/// 506464946133393486072777102926336625944849939610982267859828541006717966526573193706126370441346337661774335955699621
pub const ATE_LOOP_COUNT: [u64; 13] = [
    0x55c5b9b57b942ae8,
    0x3d52287d3dfd424a,
    0xcf1ff9d6a543deb7,
    0x820c9c5711ceeebc,
    0x549a2d44305d20fe,
    0x50f5c131afd70235,
    0xab3596c8617c5792,
    0x830c728d80f9d78b,
    0x6a7223ee72023d07,
    0xbc5d176b746af026,
    0xe959283d8f526663,
    0xc4d2263babf8941f,
    0x3848,
];

/// FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG = true
pub const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = true;

/// FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0 =
/// 7000705447348627246181409558336018323010329260726930841638672011287206690002601216854775649561085256265269640040570922609783227469279331691880282815325569032149343779036142830666859805506518426649197067288711084398033
pub const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger832 = BigInteger832([
    0xb62ef36af72855d1,
    0x676b5cef49d290fa,
    0xd17fcf3c60947427,
    0x5b93d992bc1b2849,
    0x2171887cecd072cb,
    0x879a2873f1516f4a,
    0x8cc6856bd2cdf24e,
    0xbff4fb6644d01993,
    0x5dcbeea3e31ea667,
    0x5f256f47681649f3,
    0x2355a2b0839967fe,
    0x144ed,
    0x0,
]);

/// FINAL_EXPONENT_LAST_CHUNK_W1 =
/// 86482221941698704497288378992285180119495364068003923046442785886272123124361700722982503222189455144364945735564951562986
pub const FINAL_EXPONENT_LAST_CHUNK_W1: BigInteger832 = BigInteger832([
    0x5657b9b57b942aea,
    0x84f9a65f3bd54eaf,
    0x5ea4214e35cd127,
    0xe3cbcbc14ec1501d,
    0xf196cb845a3092ab,
    0x7e14627ad0e19017,
    0x217db4,
    0x0,
    0x0,
    0x0,
    0x0,
    0x0,
    0x0,
]);
