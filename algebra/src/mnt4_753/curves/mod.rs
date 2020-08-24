use algebra_core::{
    biginteger::BigInteger768,
    curves::models::mnt4::{MNT4Parameters, MNT4},
    field_new,
    fields::FpParameters,
    Fp2,
};

use crate::mnt4_753::{Fq, Fq2, Fq2Parameters, Fq4Parameters, FqParameters, Fr, FrParameters};

pub mod g1;
pub mod g2;

#[cfg(test)]
mod tests;

pub use self::{
    g1::{G1Affine, G1Prepared, G1Projective},
    g2::{G2Affine, G2Prepared, G2Projective},
};

pub type MNT4_753 = MNT4<Parameters>;

pub struct Parameters;

impl MNT4Parameters for Parameters {
    const TWIST: Fp2<Self::Fp2Params> = field_new!(Fq2, FQ_ZERO, FQ_ONE);
    // A coefficient of MNT4-753 G2 =
    // ```
    // mnt4753_twist_coeff_a = mnt4753_Fq2(mnt4753_G1::coeff_a * non_residue, mnt6753_Fq::zero());
    //  = (A_COEFF * NONRESIDUE, ZERO)
    //  = (26, ZERO)
    // ```
    #[rustfmt::skip]
    const TWIST_COEFF_A: Fp2<Self::Fp2Params> = field_new!(Fq2,
        G1_COEFF_A_NON_RESIDUE,
        FQ_ZERO,
    );
    // https://github.com/o1-labs/snarky/blob/9c21ab2bb23874604640740d646a932e813432c3/snarkette/mnt4753.ml
    const ATE_LOOP_COUNT: &'static [u64] = &[
        8824542903220142080,
        7711082599397206192,
        8303354903384568230,
        5874150271971943936,
        9717849827920685054,
        95829799234282493,
    ];
    const ATE_IS_LOOP_COUNT_NEG: bool = true;
    const FINAL_EXPONENT_LAST_CHUNK_1: BigInteger768 =
        BigInteger768([0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]);
    const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = true;
    const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger768 = BigInteger768([
        8824542903220142079,
        7711082599397206192,
        8303354903384568230,
        5874150271971943936,
        9717849827920685054,
        95829799234282493,
        0,
        0,
        0,
        0,
        0,
        0,
    ]);
    type Fp = Fq;
    type Fr = Fr;
    type Fp2Params = Fq2Parameters;
    type Fp4Params = Fq4Parameters;
    type G1Parameters = self::g1::Parameters;
    type G2Parameters = self::g2::Parameters;
}

// 26
pub const G1_COEFF_A_NON_RESIDUE: Fq = field_new!(
    Fq,
    BigInteger768([
        16948538951764659373,
        10775354577659735631,
        12766795894854242596,
        8684022258823474090,
        973489465296612807,
        3883945490221946200,
        16178634811223492029,
        16155746945640075033,
        17642042187059426365,
        10295720303844380352,
        13265853240981244259,
        39422991244875,
    ])
);
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger768([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
pub const FR_ZERO: Fr = field_new!(Fr, BigInteger768([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
pub const FR_ONE: Fr = field_new!(Fr, FrParameters::R);
