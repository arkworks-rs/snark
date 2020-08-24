use algebra_core::{
    biginteger::BigInteger320,
    curves::models::mnt4::{MNT4Parameters, MNT4},
    field_new,
    fields::FpParameters,
    Fp2,
};

use crate::mnt4_298::{Fq, Fq2, Fq2Parameters, Fq4Parameters, FqParameters, Fr, FrParameters};

pub mod g1;
pub mod g2;

#[cfg(test)]
mod tests;

pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};

pub type MNT4_298 = MNT4<Parameters>;

pub struct Parameters;

impl MNT4Parameters for Parameters {
    const TWIST: Fp2<Self::Fp2Params> = field_new!(Fq2, FQ_ZERO, FQ_ONE);
    // A coefficient of MNT4-298 G2 =
    // ```
    // mnt4298_twist_coeff_a = mnt4298_Fq2(mnt4298_G1::coeff_a * non_residue, mnt6298_Fq::zero());
    //  = (A_COEFF * NONRESIDUE, ZERO)
    //  = (34, ZERO)
    // ```
    #[rustfmt::skip]
    const TWIST_COEFF_A: Fp2<Self::Fp2Params> = field_new!(Fq2,
        G1_COEFF_A_NON_RESIDUE,
        FQ_ZERO,
    );

    const ATE_LOOP_COUNT: &'static [u64] = &[993502997770534912, 5071219579242586943, 2027349];
    const ATE_IS_LOOP_COUNT_NEG: bool = false;
    const FINAL_EXPONENT_LAST_CHUNK_1: BigInteger320 = BigInteger320([0x1, 0x0, 0x0, 0x0, 0x0]);
    const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = false;
    const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger320 =
        BigInteger320([993502997770534913, 5071219579242586943, 2027349, 0, 0]);
    type Fp = Fq;
    type Fr = Fr;
    type Fp2Params = Fq2Parameters;
    type Fp4Params = Fq4Parameters;
    type G1Parameters = self::g1::Parameters;
    type G2Parameters = self::g2::Parameters;
}

// 34
pub const G1_COEFF_A_NON_RESIDUE: Fq = field_new!(
    Fq,
    BigInteger320([
        9379015694948865065,
        3933863906897692531,
        7183785805598089445,
        17382890709766103498,
        3934325337380,
    ])
);
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger320([0, 0, 0, 0, 0]));
pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
pub const FR_ZERO: Fr = field_new!(Fr, BigInteger320([0, 0, 0, 0, 0]));
pub const FR_ONE: Fr = field_new!(Fr, FrParameters::R);
