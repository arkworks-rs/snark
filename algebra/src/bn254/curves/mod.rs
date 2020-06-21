use crate::bn254::*;
use algebra_core::{
    biginteger::BigInteger256,
    curves::{
        bn,
        bn::{Bn, BnParameters, TwistType},
    },
    field_new,
};
pub mod g1;
pub mod g2;

#[cfg(test)]
mod tests;

pub struct Parameters;

impl BnParameters for Parameters {
    const X: &'static [u64] = &[4965661367192848881];
    /// `x` is positive.
    const X_IS_NEGATIVE: bool = false;
    const ATE_LOOP_COUNT: &'static [i8] = &[
        0, 0, 0, 1, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 0, 0, 1, 1, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0,
        0, 1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 1, 1, 0,
        -1, 0, 0, 1, 0, 1, 1,
    ];
    /// `ate_loop_count` is positive.
    const ATE_LOOP_COUNT_IS_NEGATIVE: bool = false;
    const TWIST_MUL_BY_Q_X: Fq2 = field_new!(
        Fq2,
        field_new!(
            Fq,
            BigInteger256([
                0xb5773b104563ab30,
                0x347f91c8a9aa6454,
                0x7a007127242e0991,
                0x1956bcd8118214ec,
            ])
        ),
        field_new!(
            Fq,
            BigInteger256([
                0x6e849f1ea0aa4757,
                0xaa1c7b6d89f89141,
                0xb6e713cdfae0ca3a,
                0x26694fbb4e82ebc3,
            ])
        ),
    );
    const TWIST_MUL_BY_Q_Y: Fq2 = field_new!(
        Fq2,
        field_new!(
            Fq,
            BigInteger256([
                0xe4bbdd0c2936b629,
                0xbb30f162e133bacb,
                0x31a9d1b6f9645366,
                0x253570bea500f8dd,
            ])
        ),
        field_new!(
            Fq,
            BigInteger256([
                0xa1d77ce45ffe77c7,
                0x07affd117826d1db,
                0x6d16bd27bb7edc6b,
                0x2c87200285defecc,
            ])
        ),
    );
    const TWIST_TYPE: TwistType = TwistType::D;
    type Fp = Fq;
    type Fp2Params = Fq2Parameters;
    type Fp6Params = Fq6Parameters;
    type Fp12Params = Fq12Parameters;
    type G1Parameters = g1::Parameters;
    type G2Parameters = g2::Parameters;
}

pub type Bn254 = Bn<Parameters>;

pub type G1Affine = bn::G1Affine<Parameters>;
pub type G1Projective = bn::G1Projective<Parameters>;
pub type G2Affine = bn::G2Affine<Parameters>;
pub type G2Projective = bn::G2Projective<Parameters>;
