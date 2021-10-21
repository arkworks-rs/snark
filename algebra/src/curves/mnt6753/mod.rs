use self::{g1::MNT6G1Parameters, g2::MNT6G2Parameters};
use crate::curves::models::mnt6::{
    G1Affine as MNT6G1Affine, G1Projective as MNT6G1Projective, G2Affine as MNT6G2Affine,
    G2Projective as MNT6G2Projective, MNT6Parameters, MNT6p,
};
use crate::field_new;
use crate::{
    fields::{
        mnt6753::{
            fq::{Fq, FqParameters},
            fq3::Fq3Parameters,
            fq6::Fq6Parameters,
            Fq3, Fr,
        },
        FpParameters,
    },
    BigInteger768 as BigInteger,
};

pub mod g1;
pub mod g2;
#[cfg(test)]
mod tests;

pub struct MNT6_753Parameters;

impl MNT6Parameters for MNT6_753Parameters {
    /// The Frobenius trace of the MNT6 curve is
    /// t = 204691208819330962009469868104636132783269696790011977400223898462431810102935615891307667367766898917669754470401
    /// Our Ate pairing Miller loop count is the absolute value of the Frobenius trace minus 1
    const ATE_LOOP_COUNT: &'static [u64] = &[
        0x7a7713041ba18000,
        0x6b0344c4e2c428b0,
        0x733b714aa43c31a6,
        0x51852c8cbe26e600,
        0x86dcbcee5dcda7fe,
        0x15474b1d641a3fd,
    ];

    //Output of find_wnaf(ate_loop_count), already trimmed of leading zeros and MSB
    //starting with least significant bit
    const WNAF: &'static [i32] = &[
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, -1, 0, 0,
        1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 1, 0, 0, 0, -1, 0, 0, -1, 0, 0, 0, 1, 0, 1,
        0, -1, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, -1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, -1,
        0, -1, 0, 1, 0, 0, -1, 0, 0, 1, 0, 1, 0, 0, 0, -1, 0, 1, 0, 1, 0, 0, 0, 1, 0, -1, 0, 1, 0,
        0, 0, 0, 0, -1, 0, -1, 0, -1, 0, 0, 1, 0, -1, 0, 1, 0, 1, 0, -1, 0, 1, 0, 0, -1, 0, 1, 0,
        0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, -1, 0,
        0, -1, 0, 0, -1, 0, 0, 0, 1, 0, -1, 0, 1, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0,
        1, 0, -1, 0, 0, -1, 0, 0, 1, 0, 1, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 1, 0, -1, 0, 1, 0, 0, 1,
        0, 0, -1, 0, -1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, -1, 0, 1, 0, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 1, 0, -1, 0, -1, 0, 0, 1, 0, -1, 0, 0, -1, 0, 0, 0, -1, 0, 1, 0,
        -1, 0, 0, -1, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0, 1, 0, -1, 0, 0, -1, 0, 0, -1, 0, 0, 1, 0,
        0, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, -1, 0, 1, 0, 0, 0, 0, 1, 0, 0, -1,
        0, -1, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, -1, 0, 1, 0, 1, 0, -1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1,
        0,
    ];
    // Frobenius trace of this curve is non-negative
    const ATE_IS_LOOP_COUNT_NEG: bool = false;

    const TWIST: Fq3 = field_new!(Fq3, FQ_ZERO, FQ_ONE, FQ_ZERO);

    // I would do the hard coded definition inside G2, and just refer to from here.
    const TWIST_COEFF_A: Fq3 = field_new!(
        Fq3,
        FQ_ZERO,
        FQ_ZERO,
        field_new!(
            Fq,
            BigInteger([
                // = COEFF_A
                0x4768931cfff9c7d4,
                0xc45e46d6ada96ca0,
                0x479b0bdb0b3c0107,
                0x362a089610f8d41b,
                0xdbafcec2c8a91aaf,
                0x78428b0ff9d96a06,
                0xf2e4472a9080c353,
                0xc9006ed33f0e971c,
                0x0794d9d10bdb7288,
                0x3c1e44cab5419e2c,
                0x49b5fc6c81f4560c,
                0x1c287777c30ba,
            ])
        ),
    );

    // m_1 = 1
    const FINAL_EXPONENT_LAST_CHUNK_1: BigInteger =
        BigInteger([0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]);
    // m_0 = 691208819330962009469868104636132783269696790011977400223898462431810102935615891307667367766898917669754470400
    const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger = BigInteger([
        0x7a7713041ba18000,
        0x6b0344c4e2c428b0,
        0x733b714aa43c31a6,
        0x51852c8cbe26e600,
        0x86dcbcee5dcda7fe,
        0x15474b1d641a3fd,
        0x0,
        0x0,
        0x0,
        0x0,
        0x0,
        0x0,
    ]);
    // m_0 is positive
    const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = false;

    type Fp = Fq;
    type Fr = Fr;
    type Fp3Params = Fq3Parameters;
    type Fp6Params = Fq6Parameters;
    type G1Parameters = MNT6G1Parameters;
    type G2Parameters = MNT6G2Parameters;
}

pub type MNT6 = MNT6p<MNT6_753Parameters>;
pub type G1Affine = MNT6G1Affine<MNT6_753Parameters>;
pub type G1Projective = MNT6G1Projective<MNT6_753Parameters>;
pub type G2Affine = MNT6G2Affine<MNT6_753Parameters>;
pub type G2Projective = MNT6G2Projective<MNT6_753Parameters>;

pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
