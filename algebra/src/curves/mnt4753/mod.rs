use self::{g1::MNT4G1Parameters, g2::MNT4G2Parameters};
use crate::curves::models::mnt4::{
    G1Affine as MNT4G1Affine, G1Projective as MNT4G1Projective, G2Affine as MNT4G2Affine,
    G2Projective as MNT4G2Projective, MNT4Parameters, MNT4p,
};
use crate::field_new;
use crate::{
    fields::{
        mnt4753::{
            fq::{Fq, FqParameters},
            fq2::Fq2Parameters,
            fq4::Fq4Parameters,
            Fq2, Fr,
        },
        FpParameters,
    },
    BigInteger768 as BigInteger,
};

pub mod g1;
pub mod g2;
#[cfg(test)]
mod tests;

pub struct MNT4_753Parameters;

impl MNT4Parameters for MNT4_753Parameters {
    /// The Frobenius trace of the MNT4 curve is
    /// t = -204691208819330962009469868104636132783269696790011977400223898462431810102935615891307667367766898917669754470399
    /// Our Ate pairing Miller loop count is the absolute value of the Frobenius trace minus 1
    const ATE_LOOP_COUNT: &'static [u64] = &[
        0x7a7713041ba18000,
        0x6b0344c4e2c428b0,
        0x733b714aa43c31a6,
        0x51852c8cbe26e600,
        0x86dcbcee5dcda7fe,
        0x15474b1d641a3fd,
    ];

    /// Output of find_wnaf(ate_loop_count), already trimmed of leading zeros and MSB,
    /// starting with least significant bit
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

    /// Frobenius trace of this curve is negative
    const ATE_IS_LOOP_COUNT_NEG: bool = true;

    const TWIST: Fq2 = field_new!(Fq2, FQ_ZERO, FQ_ONE);

    // I would do the hard coded definition inside G2, and just refer to from here.
    const TWIST_COEFF_A: Fq2 = field_new!(
        Fq2,
        field_new!(
            Fq,
            BigInteger([
                // = COEFF_A
                0xeb354e6121cdccad,
                0x9589bfe5ea49ae4f,
                0xb12cc53998b3d124,
                0x7883d83c06c22baa,
                0xd828782cb96edc7,
                0x35e68bd867a8d558,
                0xe0860ea489bec5bd,
                0xe034be400ffa8f19,
                0xf4d51fe5c821f43d,
                0x8ee1bf11396bcec0,
                0xb819c73cb726c963,
                0x23dae1639e4b,
            ])
        ),
        FQ_ZERO,
    );

    // m_1 = 1
    const FINAL_EXPONENT_LAST_CHUNK_1: BigInteger =
        BigInteger([0x1, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0]);
    // |m_0| =  204691208819330962009469868104636132783269696790011977400223898462431810102935615891307667367766898917669754470399
    const FINAL_EXPONENT_LAST_CHUNK_ABS_OF_W0: BigInteger = BigInteger([
        0x7a7713041ba17fff,
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
    //sign of m_0 is negative
    const FINAL_EXPONENT_LAST_CHUNK_W0_IS_NEG: bool = true;

    type Fp = Fq;
    type Fr = Fr;
    type Fp2Params = Fq2Parameters;
    type Fp4Params = Fq4Parameters;
    type G1Parameters = MNT4G1Parameters;
    type G2Parameters = MNT4G2Parameters;
}

pub type MNT4 = MNT4p<MNT4_753Parameters>;
pub type G1Affine = MNT4G1Affine<MNT4_753Parameters>;
pub type G1Projective = MNT4G1Projective<MNT4_753Parameters>;
pub type G2Affine = MNT4G2Affine<MNT4_753Parameters>;
pub type G2Projective = MNT4G2Projective<MNT4_753Parameters>;

// field element 0 in Montgomery representation
pub const FQ_ZERO: Fq = field_new!(Fq, BigInteger([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]));
// field element 1 in Montgomery representation
pub const FQ_ONE: Fq = field_new!(Fq, FqParameters::R);
